#!/bin/bash
source activate seurat-2.3.4
Rscript /home/$USER/Data/bin/R/dropSetup.R $1 $2 $3 &> output.log &
