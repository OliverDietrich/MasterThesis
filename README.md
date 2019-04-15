# MasterThesis

Analysis of single-cell RNA-seq data for the Single Cell Analysis (SIGA) group of the Helmholtz Centre for Infection Research ([HZI](https://www.helmholtz-hzi.de/en/)) Institute for RNA-based Infection Research ([HIRI](https://www.helmholtz-hiri.de/)). 

This public repository contains the source code that has been used to analyze data from the following projects.  

## Projects:
1. Transcriptional Profiling of the Immune Response to Influenza Vaccination
2. Dissecting the Cellular Composition of Pediatric Kidney Tumors
3. Exploring the Transcriptional Regulation in *Trypanosoma brucei*
4. Dissecting the Microenvironment of the Squamo-Columnar Junction of the Cervix
5. Exploring the Cellular Repertoire of Myeloid-Derived Suppressor Cells

## Workflows
The data analysis workflow developed here is based on [Seurat](https://satijalab.org/seurat/), [Scanpy](https://scanpy.readthedocs.io/en/stable/), [Scran](https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html), the [single-cell RNA-seq course](https://hemberg-lab.github.io/scRNA.seq.course/index.html) from the University of Cambridge and [RNA-seq Data Analysis - A Practical Approach](https://doi.org/10.1201/b17457). 

### Preprocessing Raw-Sequencing Data & Constructing the Expression Matrix
The processing of next-generation sequencing data is similar between single-cell and traditional RNA-seq and consists of the following steps:

1. Demultiplexing and conversion to FASTQ-format

The primary output for data analysis are per-cycle BCL basecall files from Illumina Sequencers. First, these files need to be de-multiplexed and converted to FASTQ-format by the program [bcl2fastq](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf). 

2. Quality control & preprocessing

The quality of reads has to be assessed by programs like [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Low quality reads can be either removed or trimmed by programs like [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/), [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) or [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). 

3. Aligning reads to a reference

The remaining high quality reads must be aligned to a reference in order to generate gene counts. Different aligners ([Bowtie](http://bowtie-bio.sourceforge.net/index.shtml), [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml), [STAR](https://github.com/alexdobin/STAR), [Kallisto](https://pachterlab.github.io/kallisto/about)) are available BUT MORE INFORMATION WILL BE ADDED.

4. Annotation based quality control

5. Quantitation of gene expression

### Biological Analysis
The Biological Analysis of single-cell RNA-seq data starts with an expression matrix where the columns represent barcodes (cells) and the rows represent features (genes). The challenges in the general workflow have been described well by [Kiselev, Andrews and Hemberg](https://www.nature.com/articles/s41576-018-0088-9). 

The script for the analysis of Drop-Seq data is currently based on the [Seurat](https://satijalab.org/seurat/) workflow and executed in [R](https://www.r-project.org/). 

The steps for the analysis include:

1. Quality Control (QC)
    - Adding Metadata
    - Filtering
2. Feature Selection
    - Normalization
    - Variance based selection of genes
    - Principle component analysis ([PCA](https://doi.org/10.1038/nmeth.4346))
3. Dimensional Reduction
    - t-distributed Stochastic Neighbor Embedding ([tSNE](https://lvdmaaten.github.io/tsne/))
    - Fast interpolation-based t-SNE ([FI-tSNE](https://doi.org/10.1038/s41592-018-0308-4))
    - Uniform Manifold Approximation and Projection ([UMAP](https://umap-learn.readthedocs.io/en/latest/))
    - Diffusion Map
4. Clustering
    - [Louvain](https://perso.uclouvain.be/vincent.blondel/research/louvain.html) method for community detection in large networks
    - Smart local moving algorithm for large-scale modularity-based community detection ([SLM](http://www.ludowaltman.nl/slm/))
    - Hierarchical Density-Based Spatial Clustering of Applications with Noise ([HDBSCAN](https://hdbscan.readthedocs.io/en/latest/index.html)).
5. Gene Set Enrichment Analysis
based on annotated clusters
    - Wilcoxon rank sum test  
