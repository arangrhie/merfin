#!/bin/bash

bcftools concat -Oz $(cd vcfs && ls *vcf | sort -n | awk '{printf "vcfs/"$1" "}') -o polished_${1%.*}.vcf.gz

bcftools index polished_${1%.*}.vcf.gz

bcftools norm -Oz polished_${1%.*}.vcf.gz -f ${1} -o polished_${1%.*}.norm.vcf.gz -c w

rm polished_${1%.*}.vcf.gz polished_${1%.*}.vcf.gz.csi

bcftools index polished_${1%.*}.norm.vcf.gz

bcftools consensus polished_${1%.*}.norm.vcf.gz -f ${1} -Hla > polished_${1}

