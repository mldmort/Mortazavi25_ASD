#!/bin/bash

DIR_IN=<path to processed vcf directory>
VCF_IN=$DIR_IN/longtr_merge.sorted.vep.addInfo.vcf.gz

bcftools query -l $VCF_IN > sample_order.tsv
bcftools query -f '%ID[\t%GB]' $VCF_IN > genotypes.tsv

awk 'BEGIN{FS="\t";OFS="\t"}\
{\
	if (NR==FNR) {\
		samples[NR] = $1;\
		sample_num[NR] = 0;\
	}\
	else {\
		for (i=2; i<=NF; i++) {\
			if ($i != ".") {\
				sample_num[i-1] += 1;\
			}\
		}\
	}\
}\
END {\
	for (i=1;i<=length(sample_num);i++) {\
		print samples[i], sample_num[i];\
	}\
}' sample_order.tsv genotypes.tsv | sort -k2,2n > gt_counts.tsv
