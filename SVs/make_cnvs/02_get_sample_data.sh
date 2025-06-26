#!/bin/bash

TAB_IN=/expanse/projects/sebat1/miladm/UCSD/LONG_READ_COHORT/variant_prior/SV_process/plots_svs/table_SVs_SQ20_SAMPLES_freq1.tsv
COL_OFFSET=4

FILE_IN=outputs/out_shanta_all.txt
OUT=outputs/cnvs_in_svs_shanta_all.tsv
echo "+++++++++++++++++++++++++ $FILE_IN +++++++++++++++++++++++++++"
head -n 1 $TAB_IN |
awk 'BEGIN{FS="\t";OFS="\t"}\
{\
print $1,$2,$3,$4,$5,$6,$7,$8,$37,$38,$40,$41;\
}' > $OUT

awk -v offset=$COL_OFFSET \
'BEGIN{FS="\t";OFS="\t"}\
{\
print $(1+offset),$(2+offset),$(3+offset),$(4+offset),$(5+offset),$(6+offset),$(7+offset),$(8+offset),$(37+offset),$(38+offset),$(40+offset),$(41+offset);\
}' $FILE_IN >> $OUT

FILE_IN=outputs/out_clinical_all.txt
OUT=outputs/cnvs_in_svs_clinical_all.tsv
echo "+++++++++++++++++++++++++ $FILE_IN +++++++++++++++++++++++++++"
head -n 1 $TAB_IN |
awk 'BEGIN{FS="\t";OFS="\t"}\
{\
print $1,$2,$3,$4,$5,$6,$7,$8,$37,$38,$40,$41;\
}' > $OUT

awk -v offset=$COL_OFFSET \
'BEGIN{FS="\t";OFS="\t"}\
{\
print $(1+offset),$(2+offset),$(3+offset),$(4+offset),$(5+offset),$(6+offset),$(7+offset),$(8+offset),$(37+offset),$(38+offset),$(40+offset),$(41+offset);\
}' $FILE_IN >> $OUT
