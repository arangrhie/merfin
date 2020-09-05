#!/bin/sh

if [[ "$#" -lt 4 ]]; then
  echo "Usage: ./_submit_merfin.sh full-asm.fasta variant-call.bcf sample peak memAsm memRead"
  echo -e "Submit merfin to job arrays on slurm"
	echo
  echo -e "  Requires asm.meryl and read.meryl in this path"
  echo -e "  full-asm.fasta    target assembly to polish"
	echo -e "  variant-call.bcf  variants to filter. will be split to 100 .vcf and distributed."
	echo -e "  sample            job name prefix and output prefix"
	echo -e "  peak              peak to use for normalizing copy numbers in read.meryl"
	echo -e "  memAsm            maximum memory allowed for loading asm.meryl"
  echo -e "  memRead           maximum memory allowed for loading read.meryl"
	echo
	echo -e "ver. 2020-09-05"
  exit -1
fi

if [[ ! -e asm.meryl ]] || [[ ! -e read.meryl ]]; then
  echo "Cannot find asm.meryl and read.meryl. Exit."
  exit -1
fi

asm=$1
vcf=$2
sample=$3
peak=$4
mem1=$5
mem2=$6
total_mem=$(($mem1 + $mem2 + 50))
module load samtools

if [[ ! -e $asm.fai ]]; then
  samtools faidx $asm
fi

if [[ ! -s $vcf.csi ]]; then
  bcftools index $vcf
fi

cut -f1 $asm.fai > $sample.list
LEN=`wc -l $sample.list | awk '{print $1}'`
if [[ $LEN -lt 100 ]]; then
  echo "$LEN < 100"
  MAX=$LEN
else
  MAX=100
fi

mkdir -p vcf
mkdir -p fa
mkdir -p merfin
mkdir -p logs

cpus=24
mem=${total_mem}g
name=$sample.merfin
script=merfin.sh
args="$asm $vcf $asm.fai $peak $mem1 $mem2"
walltime=8:00:00
log=logs/$name.%A_%a.log

echo "\
sbatch --partition=norm --array=1-$MAX -D $PWD --cpus-per-task=$cpus -J $name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=norm --array=1-$MAX -D $PWD --cpus-per-task=$cpus -J $name --mem=$mem --time=$walltime --error=$log --output=$log $script $args > merfin.jid
wait_for="--dependency=afterok:`cat merfin.jid`"

cpus=4
mem=8g
name=$sample.consensus
script=consensus.sh
args="$asm $sample.list $sample"
partition="quick"
walltime=4:00:00
log=logs/$name.%A.log

echo "\
sbatch --partition=$partition -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=$partition -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args > consensus.jid


