#!/bin/bash
source activate seurat_v3
Rscript /home/oliver/Data/bin/R/dropSetup.R $1 $2 $3 &> output.log &
