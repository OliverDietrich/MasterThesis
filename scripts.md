# Computational Workflow

## Bash scripts for the Raw Data Preprocessing
The preprocessing starts from raw sequencing data which requires both large storage space and computing power. Therefore, all preprocessing steps are performed on the HZI bioinf cluster which requires submission of processes to the Sun Grid Engine ([SGE](https://en.wikipedia.org/wiki/Oracle_Grid_Engine)). 

### Cell Ranger
The [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) pipeline is used for preprocessing of 10x Genomics Chromium droplet sequencing data. 

### dropEst Pipeline
[dropEst](https://dropest.readthedocs.io/en/latest/index.html) is a "pipeline for estimating molecular count matrices for droplet-based single-cell RNA-seq measurements". It serves as a more flexible alternative to Cell Ranger and is not specific for 10x Chromium libraries but can also handle inDrop, iCLIP, SPLiT-seq, Seq-Well and Drop-seq. Additionally, it includes UMI count correction and cell quality classification.

## R scripts for the Biological Analysis

### The filesystem (directory tree)
Since all scripts are executable from the command line ([bash](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) terminal) they rely on a fixed architecture of one branch of the filesystem. 

The general path to all data is:

> /home/user/Data/project/datasets

This can either be the directory where Cell Ranger (or dropEst) has been deposited or a clone of this directory containing the raw and/or filtered count matrices and optionally the web summary, metrics summary and loupe file. 

The datasets from 10x Genomics are labelled by indices ranging from A1 to H12. In all scripts the dataset is abbreviated by DS representing such indices (e.g. G3). The minimal directory conforms to:

> DS/outs/
> - filtered_gene_bc_matrix or raw_gene_bc_matrix
>   - matrix.mtx
>   - barcodes.tsv
>   - features.tsv
> - web_summary.html, metrics_summary.csv, cloupe.cloupe

Setting up the Seurat object ([dropSetup](https://github.com/OliverDietrich/MasterThesis/blob/master/dropSetup.sh)) will create the directory

> /home/user/Data/project/analysis/DS_YYYY-MM-DD

However, to files are needed to start the analysis. A tsv file called [datasets](https://github.com/OliverDietrich/MasterThesis/blob/master/datasets) in the datasets-folder that specifies the metadata which is specific for the datasets. And a tsv file called [projects](https://github.com/OliverDietrich/MasterThesis/blob/master/projects) in the Data-folder that specifies all the datasets that belong to a project and allows for access from root. Additionally, a directory "markers" can be put into the projects-folder that contains the [markergenes.csv](https://github.com/OliverDietrich/MasterThesis/blob/master/markergenes.csv) which can be used to only print plots showing the expression of interesting markergenes instead of highly variable genes (HVG) or all genes.

### 1. [Setting up](https://github.com/OliverDietrich/MasterThesis/blob/master/dropSetup.R) the R dataset
The Cell Ranger output can be imported as a matrix (.mtx) or in the hierarchical data format ([HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format), .h5). 

For the matrix import two files containing the feature (gene) names and barcodes are required. The format needs to be checked as it varies between Cell Ranger versions. The matrix rownames are specified as Ensembl IDs (no redundancy) and can be converted using the features.tsv file. For ease of use a csv-file containing visualization parameters is produced which is passed on to all further scripts.

> HDF5 import does not work in the current version!

### 2. Quality Control

### 3. Feature Selection

### 4. Dimensional Reduction

### 5. Unsupervised Clustering

### 6. Expression Plots

### 7. Labelling

### 8. Gene Set Enrichment Analysis

### 9. Gene Ontology

### Using bash scripts to submit R scripts as background processes
R scripts are run in sequence and submitted as a background process using a bash script ([dropseq](https://github.com/OliverDietrich/MasterThesis/blob/master/dropseq)).

The order is specified above ranging from 2-9 starting with the Quality Control. 

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

> conda env create -f seurat-2.3.4.yml

Sometimes a package causes trouble during the installation. It can be removed from the .yml-file manually. If it was not essential this will have fixed the issue, otherwise it may need re-installation. 

If packages are missing and the scripts abort with error warnings printed to output.log the missing packages have to be installed using either conda install r-"package" or from within R install.packages("package").
