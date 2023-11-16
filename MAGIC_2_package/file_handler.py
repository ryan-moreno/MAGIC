#!/usr/bin/env python

from pathlib import Path
import os

###################################################################################


def list_filenames(folder, allowed_extensions):
    # get all files in folder
    files = os.listdir(folder)
    # subset files with correct extension
    subset_files = []
    for ext in allowed_extensions:
        subset_files += [f for f in files if f.endswith(ext)]
    # remove duplicates
    subset_files = list(set(subset_files))

    return files


###################################################################################


def get_dir(_path):
    # get parent directory
    parent_dir = Path(_path).parent
    return parent_dir


###################################################################################


def make_magic_folder(_output_folder):
    # make results folder name

    results_dir = "MAGIC"
    counter = 2

    while Path(_output_folder / results_dir).exists():
        results_dir = f"MAGIC_{counter}"
        counter += 1

    Path.mkdir(_output_folder / results_dir)

    MAGIC_output_dir = str(_output_folder / results_dir)

    return MAGIC_output_dir


###################################################################################


def make_subfolders(output_folder, analysis_name):
    # make results folder for gene list
    # using Path from pathlib to generate paths without having to worry about OS
    MAGIC_folder = Path(output_folder)
    results_folder = Path(MAGIC_folder / analysis_name)
    Path.mkdir(results_folder)

    # make auxilary folder for CDF folder and targets folder
    auxilary_folder = Path(results_folder / "Auxiliary_Files")
    Path.mkdir(auxilary_folder)

    # make target gene folder for gene list
    targets_folder = Path(auxilary_folder / "Target_Data")
    Path.mkdir(targets_folder)

    # make figures folder for gene list
    dist_folder = Path(auxilary_folder / "Distributions")
    Path.mkdir(dist_folder)

    # make dict to hold paths
    paths = dict()
    paths["results_folder"] = str(results_folder)

    # make file paths for gene list
    paths["target_fol"] = targets_folder
    paths["distributions"] = dist_folder
    paths["summary"] = f"{results_folder}/{analysis_name}_summary.xlsx"
    paths["details"] = f"{results_folder}/{analysis_name}_details.xlsx"
    paths["Factor_gmx"] = f"{results_folder}/{analysis_name}_targets.gmx"
    paths["summary_figure"] = f"{results_folder}/{analysis_name}_summary.html"

    return paths

    return str(results_folder)


if __name__ == "__main__":
    paths = file_handler()
    print(paths)
