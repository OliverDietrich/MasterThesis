#!/bin/bash

project=$1
samples=$2
password=$3

for i in $samples
do
mkdir $i
mkdir $i/outs
#mkdir $i/outs/filtered_feature_bc_matrix
#mkdir $i/outs/raw_feature_bc_matrix

scp -r odietric@bioinf001:/vol/projects/odietric/$project/datasets/$i/outs/\{web_summary.html,metrics_summary.csv,cloupe.cloupe,raw_feature_bc_matrix,filtered_feature_bc_matrix\} $i/outs/

#scp -r odietric@bioinf001:/vol/projects/odietric/$project/datasets/$i/outs/raw_feature_bc_matrix \
#$i/outs/

#scp -r odietric@bioinf001:/vol/projects/odietric/$project/datasets/$i/outs/filtered_feature_bc_matrix \
#$i/outs/

#scp odietric@bioinf001:/vol/projects/odietric/$project/datasets/$i/outs/metrics_summary.csv \
#$i/outs/metrics_summary.html

#scp odietric@bioinf001:/vol/projects/odietric/$project/datasets/$i/outs/cloupe.cloupe \
#$i/outs/cloupe.cloupe
done
