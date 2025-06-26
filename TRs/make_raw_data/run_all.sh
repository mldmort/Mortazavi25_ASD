#!/bin/bash

set -e 
date

mkdir -p logs

./00_make_vep_table.sh 2>&1 | tee logs/out_vep.txt

python ./01_make_score_table.py 2>&1 | tee logs/out_score.txt

python ./02_make_combined_table.py 2>&1 | tee logs/out_combine.txt

date
