# MasterThesis

Analysis of single-cell RNA-seq data for the Single Cell Analysis (SIGA) group of the Helmholtz Centre for Infection Research ([HZI](https://www.helmholtz-hzi.de/en/)) Institute for RNA-based Infection Research ([HIRI](https://www.helmholtz-hiri.de/)). 

This public repository contains the source code that has been used to analyze the data. 

## Projects:
1. Transcriptional Profiling of the Immune Response to Influenza Vaccination
2. Dissecting the Cellular Composition of Pediatric Kidney Tumors
3. Exploring the Transcriptional Regulation in *Trypanosoma brucei*
4. Dissecting the Microenvironment of the Squamo-Columnar Junction of the Cervix
5. Exploring the Cellular Repertoire of Myeloid-Derived Suppressor Cells

## Workflows
The data analysis frameworks, workflows and tools that were used to develop this script are [Seurat](https://satijalab.org/seurat/), [Scanpy](https://scanpy.readthedocs.io/en/stable/), [Scran](https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html) and the [single-cell RNA-seq course](https://hemberg-lab.github.io/scRNA.seq.course/index.html) from the University of Cambridge. 

### Preprocessing Raw-Sequencing Data & Constructing the Expression Matrix
1. The primary output for data analysis are per-cycle BCL basecall files from Illumina Sequencers. First, these files need to be de-multiplexed and converted to FASTQ-format by the program [bcl2fastq](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf). 
2. The quality of reads has to be assessed by programs like [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). 
3. Low quality reads can be either removed or trimmed by programs like [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/), [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) or [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). 
4. The remaining high quality reads must be aligned to a reference in order to generate gene counts. Different aligners ([Bowtie](http://bowtie-bio.sourceforge.net/index.shtml), [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml), [STAR](https://github.com/alexdobin/STAR), [Kallisto](https://pachterlab.github.io/kallisto/about)) are available BUT MORE INFORMATION WILL BE ADDED.

### Biological Analysis
The Biological Analysis of single-cell RNA-seq data starts with an expression matrix where the columns represent barcodes (cells) and the rows represent features (genes). The challenges in the general workflow have been described well by [Kiselev, Andrews and Hemberg](https://www.nature.com/articles/s41576-018-0088-9). 

#### Chromium droplet-sequencing data analysis
The script for the analysis of Drop-Seq data is currently based on the [Seurat](https://satijalab.org/seurat/) workflow and executed in [R](https://www.r-project.org/). 
The steps for the analysis include:
1. Adding Metadata
2. Filtering
3. Normalization
4. Feature Selection
5. Dimensional Reduction
6. Clustering
7. Gene Enrichment Analysis

#### Smart-Seq2 data analysis
No Smart-Seq2 data has been processed yet.
