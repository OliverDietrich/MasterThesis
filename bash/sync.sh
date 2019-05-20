#!/bin/bash

project=$1
samples=$2

for i in $samples
do
mkdir $i
mkdir $i/outs

scp -r odietric@bioinf001:/vol/projects/odietric/$project/datasets/$i/outs/\{web_summary.html,metrics_summary.csv,cloupe.cloupe,raw_feature_bc_matrix,filtered_feature_bc_matrix\} $i/outs/
done
