#!/bin/bash

if [[ $1 =~ \.gz$ ]]; then

	output=${1%%.*}
	gunzip -c ${1} > ${1%.*}

else

	output=${1%.*}

fi

grep -v "#" ${output}.vcf | sed 's/,/;/g' > ${output}.temp.reshaped.vcf

bcftools view -h ${1} > ${output}.temp.reshaped.header.vcf

cat ${output}.temp.reshaped.header.vcf ${output}.temp.reshaped.vcf > ${output}.temp.reshaped.combined.vcf

rm ${output}.temp.reshaped.header.vcf ${output}.temp.reshaped.vcf

bcftools annotate -h $merfin/extra_header.vcf ${output}.temp.reshaped.combined.vcf > ${output}.temp.reshaped.vcf

bcftools view -h ${output}.temp.reshaped.vcf | sed 's/\tINFO/\tINFO\tFORMAT\tIND/g' > ${output}.reshaped.vcf

rm ${output}.temp.reshaped.vcf 

bcftools view -H ${output}.temp.reshaped.combined.vcf | awk -F"\t" -v OFS="\t" '{gsub(/DP=/,".\tGT:DP\t1/1:",$8);print $0}' >> ${output}.reshaped.vcf

bcftools view ${output}.reshaped.vcf -Oz > ${output}.reshaped.vcf.gz

rm ${output}.reshaped.vcf 
rm ${output}.temp.reshaped.combined.vcf
