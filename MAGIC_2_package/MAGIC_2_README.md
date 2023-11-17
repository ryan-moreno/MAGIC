# MAGIC

The manuscript accompanying this script is:
MAGIC: A tool for predicting transcription factors and cofactors driving gene sets using ENCODE data.
Roopra, A. 2020
PLoS Computational Biology



MAGIC_2_ requires Matrix files that can be downloaded from:
https://go.wisc.edu/magic

Download Matrices.zip and unpack.
Place unzipped folder in the same directory/folder as _magic.py prior to running MAGIC_2_.py

Implementation of MAGIC_2.
MAGIC2 accepts command line inputs.  If none are provided, a file path to the lists file will be requested by the script upon running.
Command line options are:

-f	File path to lists file [REQUIRED].
-p:	Maximum acceptable adjusted p value for figure generation. Default=0.25  [OPTIONAL].
-t:	Minimum number of target genes for Factor to be called. Default=10 [OPTIONAL].
-m:	Path to matrix file. Default = /Matrices/1Kb_gene.pkl.gz [OPTIONAL].
-g:	Keep General Transcription Factors and RNApol in analysis (y/n).  Default = y [OPTIONAL].  
-z:	Keep zero ChiP values (y/n). Default = n [OPTIONAL]

NOTE on z option:  Keeping zero Chip signals will recreate MAGIC outputs very similar to MAGIC outlined in Roopra,2020.  Eliminating zero chip values will place more weight on Factors with higher ChIP values and de-emphasize those Factors that may regulate genes via moderate/lower binding affinities.  This tends to produce outputs more like - but not the same as - tests based on Fisher Exact Tests with gene lists and gene sets such as EnrichR.

A tab delimited text file is requested by MAGIC (lists file).  The first column contains a list of all genes expressed in the experiment (background list).  Any number of other columns are then added containing query lists.  For example, a query list may be all genes that go up under some criterion and another may be all genes that go down.  The first row is the header and must have unique names for each column.  These names will be used to generate file and folder names so do not use special characters, spaces, period etc in column names.

MAGIC analyzes each query list and generate a series of output files and sub-directories. A series of directories are generated named after each query list in the lists file.  Each directory contains sub-directories and files for the stand-alone analysis of that query list. A Query_List_Details.xlsx file contains statistical information for all non-triaged ENCODE experiments. Data reported in Query_List_Details.xls are:

Experiment: The ENCODE experiment.  format = cell-line_experiment-ID:Factor e.g K562_ENCFF558VPP:REST
D:					Kolmogorov Smirnov statistic
argD:				Argument (ChipValue) of D
p:					Kolmogorov Smirnov p value
Score:			-log(p)xD
TF:					Factor name 

A summary file contains the best scoring Experiment for each Factor.  Each Factor appears once in this file.  The format is:
TF:					Factor Name
Experiment: The best scoring ENCODE experiment for a Factor.  format = cell-line_experiment-ID:Factor e.g K562_ENCFF558VPP:REST
D:					Kolmogorov Smirnov statistic
argD:				Argument (ChipValue) of D
p:					Kolmogorov Smirnov p value
Score:			-log(p)xD
padj:				Benjamini-Hochberg corrected p value 

Query_list_Drivers.gmx is a tab delimited file in the GMX format utilized by GSEA.  The first column contains the background list of genes.  Subsequent columns contain all target genes of each Factor.  Target genes are defined as those genes in the query list whose ChIP signal is greater than the argD i.e. the ChIP at which there is the maximal difference between the population and query cumulative.

An Auxiliary_Files directory contains 2 sub directories:
'Distributions' contains html files showing the PDFs and CDFs for ChIP values in the query and master list for each gene in the ENCODE experiment
'Target_Data_' contains comma separated text files for each Factor with padj less than user defined or default cutoff.  Each file contains a list of target genes for that Factor and associated MAGIC Matrix ChIP value.

## Magic with Generic Regions

### Creating the magic matrix

To run Magic on generic regions (rather than genes), you will need to create the magic matrix relating your regions to the Encode datasets. For this purpose, `create_magic_matrix.py` is included. `create_magic_matrix.py` will scrape Encode datasets (which you must download, instructions below) to create a matrix relating your regions to these datasets. The matrix will be saved as a python Pickle file at `Matrices/[matrix_name].pkl.gz`. The encode data is expected to be saved to `Encode_Data/`
Command line options are:

 -b:\tPath to bed file of regions [REQUIRED]. This must contain all background regions.
 -m:\tMatrix name. Default = matrix_from_regions [OPTIONAL]. The matrix will be stored at `Matrices/matrix_from_regions.pkl.gz`
 -r:\tRegion type. Default = REGION [OPTIONAL]. This will be the column name of the region IDs in the matrix

To speed up the process, you can run `create_magic_matrix.py` on subsets of the encode data and then merge the resulting matrices together.

### Downloading encode datasets

The encode datasets for all chip-seq assays can be downloaded using the following query:
`https://www.encodeproject.org/matrix/?searchTerm=chip-seq&type=Experiment&target.investigated_as=transcription+factor&files.file_type=bed+narrowPeak&award.project=ENCODE`

After clicking downloaded, a single file named `files.txt` will be donwloaded. This file contains a url to each dataset. Move `files.txt` into the `Encode_Data` folder and run the command `xargs -L 1 curl -O -J -L < files.txt` to download the datasets into the `Encode_Data` folder.

### Running Magic

To run Magic on generic regions, you will need to use the `-m` flag to specify the matrix you created. You will also need to use the `-r` flag to set the expected column name for the regions of the matrix. The list file will have the regions listed as `chr[#]:[start]-[end]` rather than `[genename]`.