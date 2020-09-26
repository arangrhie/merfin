#!/bin/sh

if [[ "$#" -lt 6 ]]; then
  echo "Usage: ./merfin.sh asm.fa in.vcf asm.fai peak mem1(in GB) mem2(in GB) [array_id_i]"
  echo -e "\tasm.fa:  full assembly to polish"
  echo -e "\tin.vcf:  full vcf to polish. This step will filter 0/0 RefCalls or QUAL<=1."
  echo -e "\tasm.fai: from samtools faidx. Or len file with seqName <tab> seqLen"
  echo -e "\tpeak:    peak to use to normalize coverage to ploidy"
  echo -e "\tmem1:    maximum memory allowed for loading asm.meryl"
  echo -e "\tmem2:    maximum memory allowed for loading read.meryl"
  echo -e "\tarray_id_i: pass array job id or 1 ~ 100. DEFAULT=\$SLURM_ARRAY_TASK_ID"
  echo -e "\t\tEvery seq in %100 == \$ith line from asm.fai will be processed."
  exit -1
fi

asm=$1
vcf=$2
fai=$3
p=$4
mem1=$5
mem2=$6

i=$SLURM_ARRAY_TASK_ID
if [[ -z $i ]]; then
# Try to get i from $7
  i=$7
fi

if [[ -z $i ]]; then
  echo "No job array id provided. Give an array id ranging from 1 to 100."
  exit -1
fi

cpus=$SLURM_CPUS_PER_TASK
if [ -z $cpus ]; then
  cpus=8
fi

set -e

module load samtools

if [[ -s vcf.$i.list ]]; then
  rm vcf.$i.list
fi

LEN=`wc -l $fai | awk '{print $1}'`

mkdir -p vcf
mkdir -p fa
mkdir -p merfin

filt_vcf="noRR.QUALgt1.$i.vcf"
tmp_vcf="vcf/tmp.$filt_vcf"
filt_vcf="vcf/$filt_vcf"

extract_vcf=0;
if [[ ! -s "$tmp_vcf" ]]; then

  echo "
  # Get header first
  bcftools view -h $vcf > $tmp_vcf
  "
  bcftools view -h $vcf > $tmp_vcf
  extract_vcf=1;
fi

# Collect files for every %100 == $i th line
for j in $(seq $i 100 $LEN )
do
  region=`sed -n ${j}p $fai | awk '{print $1":0-"$2}'`

  if [[ "$extract_vcf" -eq 1 ]]; then
    # Collect vcf files
    echo "
    # Extract $region
    bcftools view --no-header -r $region --threads=$cpus -Ov $vcf >> $tmp_vcf
    "
    bcftools view --no-header -r $region --threads=$cpus -Ov $vcf >> $tmp_vcf
  fi

  # Collect fasta entries
  fa_list="fa/fa.$i.list"
  sed -n ${j}p $fai | awk '{print $1}' >> $fa_list
done

if [[ ! -s "$filt_vcf" ]]; then
  echo "
  bcftools filter -e 'GT=="RR" || QUAL<=1' --threads=$cpus -Ov $tmp_vcf >> $filt_vcf
  "
  bcftools filter -e 'GT=="RR" || QUAL<=1' --threads=$cpus -Ov $tmp_vcf >> $filt_vcf
fi
rm $tmp_vcf

target_fa="fa/to_polish.$i.fa"
if [[ ! -s $target_fa ]]; then
  echo "
  # Collect target asm to $target_fa
  samtools faidx -r $fa_list $asm > $target_fa
  "
  samtools faidx -r $fa_list $asm > $target_fa
fi

out="merfin/merfin.$i"
echo "
/data/rhiea/merfin/build/bin/merfin -vmer -memory1 $mem1 -memory2 $mem2 -sequence $target_fa -seqmers asm.meryl -readmers read.meryl -vcf $filt_vcf -peak $p -output $out"
/data/rhiea/merfin/build/bin/merfin -vmer -memory1 $mem1 -memory2 $mem2 -sequence $target_fa -seqmers asm.meryl -readmers read.meryl -vcf $filt_vcf -peak $p -output $out

if [[ -s $out.polish.vcf ]]; then
  bcftools view --threads=$cpus -Oz $out.polish.vcf > $out.polish.vcf.gz
  bcftools index $out.polish.vcf.gz

  echo "$out.polish.vcf.gz" >> vcf.$i.list

  echo "
  Clean up"
  rm $filt_vcf $target_fa $out.polish.vcf $fa_list
fi

