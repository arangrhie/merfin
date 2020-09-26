#!/bin/sh

if [[ "$#" -lt 3 ]]; then
  echo "Usage: ./consensus.sh asm.fa sample.list out"
  exit -1
fi

asm=$1
list=$2
out=$3

cpus=$SLURM_CPUS_PER_TASK
if [[ -z cpus ]]; then
  cpus=4
fi

module load samtools

cat vcf.*.list > vcf.list
echo "Concat
bcftools concat -f vcf.list --threads $cpus -Oz -o $out.merfin.vcf.gz"
bcftools concat -f vcf.list --threads $cpus -Oz -o $out.merfin.vcf.gz

echo "Index
bcftools index $out.merfin.vcf.gz"
bcftools index $out.merfin.vcf.gz

echo "Consensus
bcftools consensus -H 1 $out.merfin.vcf.gz -f $asm -c $out.merfin.chain > $out.merfin.fa"
bcftools consensus -H 1 $out.merfin.vcf.gz -f $asm -c $out.merfin.chain > $out.merfin.fa

