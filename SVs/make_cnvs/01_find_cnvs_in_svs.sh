#!/bin/bash

### If this table does not exist you should make it with code in folder: SVs
TAB_IN=../outputs/table_SVs_SQ20_SAMPLES_freq1.tsv

BED_IN=outputs/Shanta_CNVs_shanta.bed
OUT=outputs/out_shanta_all.txt
bedtools intersect -a $BED_IN -b <(awk 'BEGIN{FS="\t";OFS="\t"}$15<=1' $TAB_IN) -wa -wb -r -f 0.5 > $OUT

BED_IN=outputs/Shanta_CNVs_clinical.bed
OUT=outputs/out_clinical_all.txt
bedtools intersect -a $BED_IN -b <(awk 'BEGIN{FS="\t";OFS="\t"}$15<=1' $TAB_IN) -wa -wb -r -f 0.5 > $OUT
