#!/usr/bin/env python

import pandas as pd
import pyranges as pr
import sys
from pathlib import Path
import file_handler as fh


def get_cmd_line_args(params):
    """
    #default params.

    params = {
    "bed_path": "",
    "matrix_name": "matrix_from_regions",
    "region_type": "REGION",
    }
    """

    cmd = sys.argv
    script_path = cmd.pop(0)

    help_text = """
 -b:\tPath to bed file of regions [REQUIRED].
 -m:\tMatrix name. Default = matrix_from_regions [OPTIONAL].
 -r:\tRegion type. Default = REGION [OPTIONAL].
	"""

    if len(cmd) == 0:
        print("\n NO COMMAND LINE INPUT GIVEN. ALLOWED OPTIONS:\n")
        print(help_text)
        exit()

    # cycle through cmd line, get keys and check vals
    bad_options = []

    while len(cmd) > 0:
        _opt = cmd.pop(0).lower()

        if _opt == "-b":
            params["bed_path"] = cmd.pop(0)
            continue

        if _opt == "-m":
            params["matrix_name"] = cmd.pop(0)
            continue

        if _opt == "-r":
            params["region_type"] = cmd.pop(0)
            continue

        if _opt in ["-h", "--h", "-help", "--help"]:
            print(help_text)
            exit()

        bad_options.append(_opt)

    if len(bad_options) == 0:
        return params

    else:
        print("\n FOLLOWING TERMS NOT FOUND.")
        print(bad_options)

        print("\n ALLOWED OPTIONS:\n")
        print(help_text)
        exit()


######################################################################################################


def read_narrow_peaks_bed(file_path):
    bed_table = pr.read_bed(str(file_path))
    bed_table.columns = [
        "Chromosome",
        "Start",
        "End",
        "Name",
        "Score",
        "Strand",
        "SignalValue",
        "PValue",
        "QValue",
        "Peak",
    ]
    return bed_table


######################################################################################################

# default params.
script_path = sys.argv[0]
_dir = Path(script_path).parent

params = {
    "bed_path": "",
    "matrix_name": "matrix_from_regions",
    "region_type": "REGION",
}

# get params from command line
params = get_cmd_line_args(params)

bed_input_path = Path(params["bed_path"])
matrix_output_path = Path(f"{_dir}/Matrices/{params['matrix_name']}.pkl.gz")
encode_metadata_path = Path(f"{_dir}/Encode_Data/metadata.tsv")
encode_data_folder = Path(f"{_dir}/Encode_Data")
region_colname = params["region_type"]

# Read in metadata for Encode data
encode_metadata = pd.read_table(encode_metadata_path)
encode_metadata["transcription_factor"] = (
    encode_metadata["Experiment target"].str.split("-").str[0]
)
encode_metadata["file_id"] = encode_metadata["File accession"]

# read in bed file with regions of interest
background_regions = pr.read_bed(str(bed_input_path), as_df=True)
# the other columns are not used and may have an incorrect column name from read_bed
background_regions = background_regions[["Chromosome", "Start", "End"]]
background_regions = background_regions.sort_values(
    ["Chromosome", "Start", "End"]
).reset_index()
background_regions["Name"] = (
    background_regions.Chromosome.astype(str)
    + ":"
    + background_regions.Start.astype(str)
    + "-"
    + background_regions.End.astype(str)
)
# Store back in PyRanges object after sorting
background_regions = pr.PyRanges(background_regions)

# scrape encode data for transcription factor ChIP-seq data
magic_matrix = pd.DataFrame(data={region_colname: background_regions.Name})
filenames = fh.list_filenames(encode_data_folder, ["bed.gz", "bed"])
# TODO: just for testing; remove this
filenames = [
    filename
    for filename in filenames
    if "ENCFF002DAF" in filename or "ENCFF003AMQ" in filename
]

for idx, filename in enumerate(filenames):
    print(f"Parsing chip seq file [{idx+1}/{len(filenames)}]: {filename}")
    file_path = Path(encode_data_folder / filename)
    encode_bed = read_narrow_peaks_bed(file_path)

    exp_id = filename.split(".")[0]
    exp_metadata = encode_metadata[encode_metadata.file_id == exp_id]
    transcription_factor = exp_metadata.transcription_factor.values[0]
    exp_name = ":".join([str(exp_id), str(transcription_factor)])

    exp_regions_matrix = pd.DataFrame(
        data={region_colname: magic_matrix[region_colname], exp_name: 0}
    )

    # loop over regions of interest to find highest signal value for each region
    # TODO: replace with dataframe.apply for speed
    for index, row in exp_regions_matrix.iterrows():
        region_id = row[region_colname]
        region = background_regions[background_regions.Name == region_id]

        if index % 1000 == 0:
            print(
                f"   Parsing region [{index+1}/{len(exp_regions_matrix)}]: {region_id}"
            )

        # Select the highest chip-seq motif signal value that overlaps with the enhancers
        region_overlap = encode_bed.overlap(region)

        if region_overlap.empty:
            highest_signal_value = 0
        else:
            highest_signal_value = region_overlap.SignalValue.max()

        exp_regions_matrix.loc[index, exp_name] = highest_signal_value

    # Add column to magic matrix if signal from experiment was found in at least one region of interest
    if not (exp_regions_matrix[exp_name] == 0).all():
        magic_matrix = pd.concat([magic_matrix, exp_regions_matrix[exp_name]], axis=1)

# Pickle the magic matrix
magic_matrix.to_pickle(matrix_output_path)
