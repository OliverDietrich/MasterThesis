# Computational Workflow

## Bash scripts for the Raw Data Preprocessing
The preprocessing starts from raw sequencing data which requires both large storage space and computing power. Therefore, all preprocessing steps are performed on the HZI bioinf cluster which requires submission of processes to the Sun Grid Engine ([SGE](https://en.wikipedia.org/wiki/Oracle_Grid_Engine)). 

### dropEst Pipeline

### Cell Ranger
The [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) pipeline is used for preprocessing of 10x Genomics Chromium droplet sequencing data. 

## R scripts for the Biological Analysis

### 1. Setting up the R dataset

### 2. Quality Control

### 3. Feature Selection

### 4. Dimensional Reduction

### 5. Unsupervised Clustering

### 6. Expression Plots

### 7. Labelling

### 8. Gene Set Enrichment Analysis

### 9. Gene Ontology

### Using bash scripts to submit R scripts as background processes
R scripts are run in sequence and submitted as a background process using a bash script (dropseq).

The order is specified above ranging from 2--9 starting with the Quality Control. 

The dropseq bash script is used to call an R script by specifying the step (e.g. QC) and the dataset (e.g. G7). This will access the R script (e.g. dropseQC.R) and submit it as a background process while producing an output log in the working directory.

dropSetup.R is run from the dataset directory (/home/user/Data/Project_00000/datasets) and will create a directory based on the dataset appended with the current date (e.g. G7_2019-05-13). dropseq (a bash script) will be run from this directory (/home/user/Data/Project_00000/analysis/G7_2019-05-13). 

The R dataset will be accessed from the RDS directory. Different steps (e.g. Quality Control) will create subdirectories and fill them with visualizations. The data will be stored in the same R dataset that has been loaded, thus there will always be just one R dataset per analysis. 

## Patching
Patches are used to update the commonly used version of a script (e.g. dropseQC.R) with changes from the latest draft (e.g. dropseQC-1.7) without replacing the file.

The code used for patching such a script is:

> diff -c dropseQC.R dropseQC-1.x.R | patch

Alternatively, the patch (file containing only the changes) can be distributed and applied to update the script, for example to multiple servers.

> diff -c dropseQC.R dropseQC-1.x.R > patch-1.x

Now you only have to download the patch and run

> patch -i patch-1.x

## Conda environment
The bash scripts activate a conda environment that includes all necessary packages (dependencies) in the correct version for the execution of the script. 

The conda environment file (.yml) is available in the file collection and can be installed by the simple command

> conda env create -f Seurat-2.3.4.yml
