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

> /home/$USER/Data/project/datasets

This can either be the directory where Cell Ranger (or dropEst) has been deposited or a clone of this directory containing the raw and/or filtered count matrices and optionally the web summary, metrics summary and loupe file. 

The datasets from 10x Genomics are labelled by indices ranging from A1 to H12. In all scripts the dataset is abbreviated by DS representing such indices (e.g. G3). The minimal directory conforms to:

> DS/outs/
> - filtered_gene_bc_matrix or raw_gene_bc_matrix
>   - matrix.mtx
>   - barcodes.tsv
>   - features.tsv
> - web_summary.html, metrics_summary.csv, cloupe.cloupe

### 1. [Setting up](https://github.com/OliverDietrich/MasterThesis/blob/master/R/dropSetup.R) the R dataset
The Cell Ranger output can be imported as a matrix (.mtx) or in the hierarchical data format ([HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format), .h5). 

> HDF5 import does not work in the current version!

For the matrix import two files containing the feature (gene) names and barcodes are required. The format needs to be checked as it varies between Cell Ranger versions. The matrix rownames are specified as Ensembl IDs (no redundancy) and can be converted using the features.tsv file. For ease of use a csv-file containing visualization parameters is produced which is passed on to all further scripts.

> Make sure the files are unzipped, otherwise the script will not work
>
> gunzip DS/outs/filtered_feature_bc_matrix/*

Setting up the Seurat object will create the directory

> /home/$USER/Data/project/analysis/DS_YYYY-MM-DD

However, to files are needed to start the analysis. A tsv file called [datasets](https://github.com/OliverDietrich/MasterThesis/blob/master/docs/datasets) in the datasets-folder that specifies the metadata which is specific for the datasets. And a tsv file called [projects](https://github.com/OliverDietrich/MasterThesis/blob/master/docs/projects) in the Data-folder that specifies all the datasets that belong to a project and allows for access from root. Additionally, a directory "markers" can be put into the projects-folder that contains the [markergenes.csv](https://github.com/OliverDietrich/MasterThesis/blob/master/docs/markergenes.csv) which can be used to only print plots showing the expression of interesting markergenes instead of highly variable genes (HVG) or all genes.

### 2. [Quality Control](https://github.com/OliverDietrich/MasterThesis/blob/master/R/dropseQC.R)
Once the data is imported into R (and stored as an R dataset) the quality of all barcodes (cells) is assessed. This can be done based on either the raw or filtered matrix. However, since our filtering criteria are generally more stringent than the thresholds set by Cell Ranger it is often preferable to import the filtered matrix to reduce memory consumption and runtime.

The data matrix is filtered based on three main criteria: library size (number of unique molecular identifiers, nUMI), number of genes (nGene) and the percentage of mitochondrial genes (pMito) per barcode. The most common visualizations are violin plots (one each) and scatterplots (nGene v. nUMI, nGene v. pMito). The Scatterplots are usually the most informative. Obviously, there should be no correlation between the number of expressed genes and the percentage of mitochondrial genes whereas the number of UMI should strictly correlate to the former.

> A interactive "gating" strategy would be useful but so far filtering is done by choosing upper and lower thresholds for each of the quality metrics

Currently, the mitochondrial genes are inferred from the gene names (e.g. MT-ND1 vs mt-Nd1). Alternatively, the chromosome location could be used to select the genes used to calculate this percentage.

### 3. [Feature Selection](https://github.com/OliverDietrich/MasterThesis/blob/master/R/dropseqFS.R)
The total number of expressed genes in a population of cells often amounts to roughly 20,000. To infer cell types from transcriptomic data it is necessary to cluster the cells based on some distance metric that is derived from these features. The number of features is too high to properly analyze (see [Curse of Dimensionality](https://medium.freecodecamp.org/the-curse-of-dimensionality-how-we-can-save-big-data-from-itself-d9fa0f872335)) so we need to reduce the number of features by selecting only those which can actually explain the biological differences in question. 

The first step of selection is the calculation of so-called highly variable genes (HVG) which show a high variance of expression in the population.

The second step is a mathematical procedure called principal component analysis ([PCA](https://en.wikipedia.org/wiki/Principal_component_analysis), alternatives are [CCA](https://en.wikipedia.org/wiki/Canonical_correlation) or [SVD](https://en.wikipedia.org/wiki/Singular_value_decomposition)) which entails matrix rotations to identify principal components that explain decreasing amounts of variation in the dataset. 

### 4. [Dimensional Reduction](https://github.com/OliverDietrich/MasterThesis/blob/master/R/dropseqDR.R)
Visualization of multi-dimensional data requires the reduction from multiple (e.g. 20 principle components) to two dimensions. There are multiple algorithms that can perform this task among the most popular are t-distributed stochastic neighbor embedding ([tSNE](https://lvdmaaten.github.io/tsne/)) and uniform manifold approximation and projection ([UMAP](https://github.com/lmcinnes/umap)).

### 5. [Unsupervised Clustering](https://github.com/OliverDietrich/MasterThesis/blob/master/R/dropseqUC.R)

### 6. [Expression Plots](https://github.com/OliverDietrich/MasterThesis/blob/master/R/dropseqExpressionPlots.R)

### 7. [Labelling](https://github.com/OliverDietrich/MasterThesis/blob/master/R/dropseqLabel.R)

### 8. [Gene Set Enrichment Analysis](https://github.com/OliverDietrich/MasterThesis/blob/master/R/dropseqGSEA.R)

### 9. [Gene Ontology](https://github.com/OliverDietrich/MasterThesis/blob/master/R/dropseqGO.R)

### Using bash scripts to submit R scripts as background processes
R scripts are run in sequence and submitted as a background process using a bash script ([dropseq](https://github.com/OliverDietrich/MasterThesis/blob/master/bash/dropseq)).

The order is specified above ranging from 2-9 starting with the Quality Control. 

The dropseq bash script is used to call an R script by specifying the step (e.g. QC) and the dataset (e.g. G7). This will access the R script (e.g. dropseQC.R) and submit it as a background process while producing an output log in the working directory.

dropSetup.R is run from the dataset directory (/home/$USER/Data/project/datasets) and will create a directory based on the dataset appended with the current date (e.g. G7_2019-05-13). dropseq (a bash script) will be run from this directory (/home/$USER/Data/project/analysis/G7_2019-05-13). 

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

## Sun Grid Engine ([SGE](http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html))
On HZI bioinformatics cluster a scheduling system is used to direct jobs to the CPU. Therefore, tasks have to be submitted in a specialized manner using the command qsub and a script that specifies the parameters. 

The status of scheduled jobs can be retrieved using the command 

> qstat

The availability of CPU cores and distribution of scheduled tasks can be showed by 

> qstat -g c

## Shiny Apps
Visualization is key when it comes to data analysis. Changing plots in R is very flexible and easy but requires coding skills. Exporting image files to provide browsable data is therefore necessary when working with computationally less trained collaborators. This requires fixed export criteria which severely decrease the flexibility as well as adding a large footprint of the analysis in terms of disk usage. 

Currently, flexibility is achieved by using an external file that specifies image metrics and can be freely manipulated. This requires, however, the repeated execution of scripts and optimally the separation of scripts into computation and visualization. 

In the future, this rather tedious workflow should be replaced by using [Shiny](https://shiny.rstudio.com/) applications that produce reactive plots that can be easily manipulated. Moreover, only the R dataset is required which substantially simplifies the transfer of data while reducing the amount of blocked disk space. 

> The preliminary code for the Shiny app can be found [here](https://github.com/OliverDietrich/MasterThesis/blob/master/shinyApp.R)
