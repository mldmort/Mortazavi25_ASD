#!/bin/bash

### omitted one of the copies of 22,19049819,21112437
FILE_IN=Shanta_CNVs_mod.csv

#group = ""
#if ($1~"Shanta") {group = "Shanta"} else {group = "clinical"};\

BED_OUT=outputs/Shanta_CNVs_shanta.bed
awk 'BEGIN{FS=",";OFS="\t"}\
$1 ~ "Shanta" \
{\
	print "chr"$2, $3, $4, "Shanta";\
}' $FILE_IN | sort -k1,1V -k2,2n -k3,3n > $BED_OUT

BED_OUT=outputs/Shanta_CNVs_clinical.bed
awk 'BEGIN{FS=",";OFS="\t"}\
$1 ~ "Clinical" \
{\
	print "chr"$2, $3, $4, "Clinical";\
}' $FILE_IN | sort -k1,1V -k2,2n -k3,3n > $BED_OUT
