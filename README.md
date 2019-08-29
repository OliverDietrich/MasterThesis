# MasterThesis

Analysis of single-cell RNA-seq data for the Single Cell Analysis ([SIGA](https://www.helmholtz-hiri.de/en/research/organisation/teams/team/single-cell-analysis/)) group of the Helmholtz Centre for Infection Research ([HZI](https://www.helmholtz-hzi.de/en/)) Institute for RNA-based Infection Research ([HIRI](https://www.helmholtz-hiri.de/)). 

The data analysis workflow developed here is based on [Seurat](https://satijalab.org/seurat/) and guided by [Scanpy](https://scanpy.readthedocs.io/en/stable/), [Scran](https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html), the [single-cell RNA-seq course](https://scrnaseq-course.cog.sanger.ac.uk/website/introduction-to-single-cell-rna-seq.html) from the University of Cambridge and [RNA-seq Data Analysis - A Practical Approach](https://doi.org/10.1201/b17457). 

## Preprocessing 

> From Raw-Sequencing Data to Constructing the Expression Matrix

The library preparation differs greatly between conventional and single-cell RNA-seq and the choice of technology depends on the biological question at hand. The different single-cell technologies have been reviewed by [Ziegenhain, et al. 2017](https://doi.org/10.1016/j.molcel.2017.01.023) and [Svensson et al. 2017](https://doi.org/10.1038/nmeth.4220). 

The processing of next-generation sequencing data is similar between single-cell and traditional RNA-seq. The general workflow for the analysis of Illumina sequencing data can be broken down into the following steps:

1. Demultiplexing and conversion to [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format)-format

The primary output for data analysis are per-cycle BCL basecall files from Illumina Sequencers. First, these files need to be de-multiplexed and converted to FASTQ-format by the Illumina software [bcl2fastq](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf). 

2. Quality control & preprocessing

The quality of reads depends on the library preparation that can induce technical artifacts including low-confidence bases, sequence-specific bias, 3′/5′ positional bias, polymerase chain reaction (PCR) artifacts, untrimmed adapters, and sequence contamination ([Korpelainen, E. RNA-seq a Practical Approach](https://doi.org/10.1201/b17457)) which are inherent to RNA-seq as well as gene dropouts (Kharchenko, Silberstein, and Scadden [2014](https://doi.org/10.1038/nmeth.2967)) which are specific for the single-cell approach. 

Tools for checking read quality include [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [PRINSEQ](https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtr026) which provide a visual report of several quality metrics. Further processing of reads includes filtering and trimming which can be performed by a range of other tools including PRINSEQ, [Trimmomatic](https://doi.org/10.1093/bioinformatics/btu170), [Cutadapt](https://doi.org/10.14806/ej.17.1.200), and [FastX](http://hannonlab.cshl.edu/fastx_toolkit/index.html).

Specific for single-cell RNA-seq is the removal of barcodes derived from empty droplets prior to alignment. Even though no cell was captured the presence of ambient RNA can lead to read counts for those barcodes and this background of small libraries has to be removed. This can be achieved by setting a certain threshold, usually associated with a significant drop in library size, or by using newly developed tools for this purpose (e.g. [dropletUtils](https://bioconductor.org/packages/devel/bioc/html/DropletUtils.html)).

3. Aligning reads to a reference

The remaining high quality reads must be aligned to a reference in order to generate gene counts. There is a large variety of alignment programs available which offer distinct qualities. Some aligners ([bowtie2](https://doi.org/10.1038/nmeth.1923), [BWA](https://doi.org/10.1093/bioinformatics/btp324)) only align contiguously and are not able to handle genomes containing introns. Examples for spliced aligners are [TopHat2](https://doi.org/10.1186/gb-2013-14-4-r36), [STAR](https://doi.org/10.1093/bioinformatics/bts635), or [GSNAP](https://doi.org/10.1007/978-1-4939-3578-9_15) and the list can be extended by pseudo-aligners like [Kallisto](https://doi.org/10.1038/nbt.3519).

4. Alignment quality control

The output generated from the alignment are SAM files. These have to be converted to BAM files by [SAMtools](https://doi.org/10.1093/bioinformatics/btp352). The quality of the alignment can be assessed by [RseQC](https://doi.org/10.1093/bioinformatics/bts356).

5. Quantitation of gene expression

To quantify the gene expression level of each gene for each cell the mapped reads can be counted per genomic feature based on the location information. Several counting tools are available, however, many single-cell protocols employ unique molecular identifiers ([UMI](https://doi.org/10.1038/s41598-018-31064-7)) to count the absolute number of molecules (and get rid of amplification bias) which has to be supported by the program. A popular tool used for both conventional and single-cell RNA-seq is [HTseq](https://doi.org/10.1093/bioinformatics/btu638).

6. Correcting for Errors

Collapsing reads based on a shared UMI-sequence is necessary for UMI-based scRNAseq protocols, however, it is still an active area of research how to process and use the UMIs correctly. Popular tools available at the moment are [UMI-tools](https://genome.cshlp.org/content/27/3/491), [zUMIs](https://doi.org/10.1093/gigascience/giy059) and [dropEst](https://doi.org/10.1186/s13059-018-1449-6).

> The pre-processing of 10x Chromium Drop-Seq data has been streamlined by the [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) software.

> [dropEst](https://dropest.readthedocs.io/en/latest/index.html) is a "pipeline for estimating molecular count matrices for droplet-based single-cell RNA-seq measurements". It serves as a more flexible alternative to Cell Ranger and is not specific for 10x Chromium libraries but can also handle inDrop, iCLIP, SPLiT-seq, Seq-Well and Drop-seq. Additionally, it includes UMI count correction and cell quality classification.

## Biological Analysis
The Biological Analysis of single-cell RNA-seq data starts with an expression matrix where the columns represent barcodes (cells) and the rows represent features (genes). The challenges in the general workflow have been described well by [Kiselev, Andrews and Hemberg](https://www.nature.com/articles/s41576-018-0088-9). A recent article by [Luecken and Theis] (https://doi.org/10.15252/msb.20188746) presents the current best practices in single-cell RNA-seq analysis.

The script for the analysis of Drop-Seq data is currently based on the [Seurat](https://satijalab.org/seurat/) workflow and executed in [R](https://www.r-project.org/). 

The steps for the analysis include:

- Quality Control ([QC](https://doi.org/10.1093/bioinformatics/btw777))
    - Adding Metadata
    - Filtering
- Normalization
- [Feature Selection](https://doi.org/10.1101/574574)
    - Variance based selection of genes
    - Principle component analysis ([PCA](https://doi.org/10.1038/nmeth.4346))
- Dimensional Reduction
    - t-distributed Stochastic Neighbor Embedding ([tSNE](https://lvdmaaten.github.io/tsne/))
    - Fast interpolation-based t-SNE ([FI-tSNE](https://doi.org/10.1038/s41592-018-0308-4))
    - Uniform Manifold Approximation and Projection ([UMAP](https://umap-learn.readthedocs.io/en/latest/))
    - Diffusion Maps ([DM](https://doi.org/10.1093/bioinformatics/btv325))
    - [SPRING](https://doi.org/10.1093/bioinformatics/btx792)
- Clustering
    - [Louvain](https://perso.uclouvain.be/vincent.blondel/research/louvain.html) method for community detection in large networks
    - Smart local moving algorithm for large-scale modularity-based community detection ([SLM](http://www.ludowaltman.nl/slm/))
    - Hierarchical Density-Based Spatial Clustering of Applications with Noise ([HDBSCAN](https://hdbscan.readthedocs.io/en/latest/index.html)).
- Differential Gene Expression (DE)
based on annotated clusters
    - [Wilcoxon rank sum test](https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test)
    - [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
    - [MAST](https://bioconductor.org/packages/3.9/bioc/html/MAST.html)

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
The smart local moving algorithm for large-scale modularity-based community detection ([SLM](https://link.springer.com/article/10.1140%2Fepjb%2Fe2013-40829-0)) is employed to find density-based clusters in the data. UMAP/tSNE plots colored by cluster are printed as output. 

### 6. [Expression Plots](https://github.com/OliverDietrich/MasterThesis/blob/master/R/dropseqExpressionPlots.R)
Depending on the settings (HVG, all, markergenes file) a number of UMAP/tSNE plots colored by gene expression are printed.

### 7. [Labelling](https://github.com/OliverDietrich/MasterThesis/blob/master/R/dropseqLabel.R)
Labels for the clusters in a defined resolution are added to the cell metadata.

### 8. [Differential Expression](https://github.com/OliverDietrich/MasterThesis/blob/master/R/dropseqDE.R)
The Wilcoxon rank sum test is employed to find differentially expressed genes across the labelled clusters.

### 9. [Gene Ontology](https://github.com/OliverDietrich/MasterThesis/blob/master/R/dropseqGO.R)
From the list of differentially expressed genes (only upregulated genes) the enriched gene ontology categories are calculated and visualized as a dotplot.

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

### Conda environment
The bash scripts activate a conda environment that includes all necessary packages (dependencies) in the correct version for the execution of the script. 

The conda environment file (.yml) is available in the file collection and can be installed by the simple command

> conda env create -f seurat-2.3.4.yml

Sometimes a package causes trouble during the installation. It can be removed from the .yml-file manually. If it was not essential this will have fixed the issue, otherwise it may need re-installation. 

If packages are missing and the scripts abort with error warnings printed to output.log the missing packages have to be installed using either conda install r-"package" or from within R install.packages("package").

### Using bash scripts to submit R scripts as background processes
R scripts are run in sequence and submitted as a background process using a bash script ([dropseq](https://github.com/OliverDietrich/MasterThesis/blob/master/bash/dropseq)).

The order is specified above ranging from 2-9 starting with the Quality Control. 

The dropseq bash script is used to call an R script by specifying the step (e.g. QC) and the dataset (e.g. G7). This will access the R script (e.g. dropseQC.R) and submit it as a background process while producing an output log in the working directory.

dropSetup.R is run from the dataset directory (/home/$USER/Data/project/datasets) and will create a directory based on the dataset appended with the current date (e.g. G7_2019-05-13). dropseq (a bash script) will be run from this directory (/home/$USER/Data/project/analysis/G7_2019-05-13). 

The R dataset will be accessed from the RDS directory. Different steps (e.g. Quality Control) will create subdirectories and fill them with visualizations. The data will be stored in the same R dataset that has been loaded, thus there will always be just one R dataset per analysis. 

### Shiny App
Visualization is key when it comes to data analysis. Changing plots in R is very flexible and easy but requires coding skills. Exporting image files to provide browsable data is therefore necessary when working with computationally less trained collaborators. This requires fixed export criteria which severely decrease the flexibility as well as adding a large footprint of the analysis in terms of disk usage. 

Currently, flexibility is achieved by using an external file that specifies image metrics and can be freely manipulated. This requires, however, the repeated execution of scripts and optimally the separation of scripts into computation and visualization. 

In the future, this rather tedious workflow should be replaced by using [Shiny](https://shiny.rstudio.com/) applications that produce reactive plots that can be easily manipulated. Moreover, only the R dataset is required which substantially simplifies the transfer of data while reducing the amount of blocked disk space. 

> The preliminary code for the Shiny app can be found [here](https://github.com/OliverDietrich/MasterThesis/blob/master/shinyApp.R)

## Future development
1. Inference Regularory Networks ([SCENIC](https://github.com/aertslab/SCENIC))
2. Rare Cell Identification ([RaceID](https://doi.org/10.1038/nature14966))
3. Scoring of Batch Effect ([k-Bet](https://github.com/theislab/kBET/blob/master/README.md))
4. Removal of Batch Effect ([CCA](https://doi.org/10.1038/nbt.4096), [MNN](https://doi.org/10.1038/nbt.4091), [scMerge](https://doi.org/10.1073/pnas.1820006116))
5. Automating Feature Selection and Clustering ([GAP statistic](https://statweb.stanford.edu/~gwalther/gap))
6. Python-based workflow ([Scanpy](https://scanpy.readthedocs.io/en/latest/))
7. Bioconductor workflow ([simpleSingleCell](https://bioconductor.org/packages/release/workflows/html/simpleSingleCell.html))
8. [SCONE](https://www.bioconductor.org/packages/devel/bioc/vignettes/scone/inst/doc/sconeTutorial.html#the-scone-workflow)
9. [clusterExperiment](https://bioconductor.org/packages/release/bioc/vignettes/clusterExperiment/inst/doc/clusterExperimentTutorial.html)
10. Gene Ontology ([clusterProfiler](https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html))
11. [SCDE](http://hms-dbmi.github.io/scde/)
12. Pseudotime analysis ([Monocle](http://cole-trapnell-lab.github.io/monocle-release/), [Slingshot](https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/slingshot.html), [PAGA](https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html))
