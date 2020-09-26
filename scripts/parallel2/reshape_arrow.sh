#!/bin/bash

if [[ $1 =~ \.gz$ ]]; then

	output=${1%%.*}

else

	output=${1%.*}

fi


bcftools annotate -h $merfin/extra_header.vcf ${1} > ${output}.temp.reshaped.vcf

bcftools view -h ${output}.temp.reshaped.vcf | sed 's/\tINFO/\tINFO\tFORMAT\tIND/g' > ${output}.reshaped.vcf

rm ${output}.temp.reshaped.vcf 

bcftools view -H ${1} | awk -F"\t" -v OFS="\t" '{gsub(/DP=/,".\tGT:DP\t1/1:",$8);print $0}' >> ${output}.reshaped.vcf

bcftools view ${output}.reshaped.vcf -Oz > ${output}.reshaped.vcf.gz

rm ${output}.reshaped.vcf 
