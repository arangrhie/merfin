#!/bin/bash

bcftools concat -Oz $(cd vcfs && ls *vcf | sort -n | awk '{printf "vcfs/"$1" "}') -o polished_${1%.*}.vcf.gz

bcftools index polished_${1%.*}.vcf.gz

bcftools consensus polished_${1%.*}.vcf.gz -f ${1} > polished_${1}

