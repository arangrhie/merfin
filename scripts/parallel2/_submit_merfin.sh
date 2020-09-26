#!/bin/bash

set -e -o pipefail

#++++                  This script is part of:                    ++++
#++++                      the T2T project                        ++++
#++++     Credit: Giulio Formenti gformenti@rockefeller.edu       ++++
#++++             Arang  Rhie     arrhie@gmail.com                ++++

if [ -z $1 ]; then

	echo "use $0 -h for help"
	exit 0
elif [ $1 == "-h" ]; then

	cat << EOF
	This script evaluates variant calls based on the K* metric.
	
	Required inputs are:
	-a	assembly file
	-v 	variant file
	-d	assmebly meryl db
	-r	read meryl db
	-p	haploid peak
	-q	QUAL variant filter
	-o	output prefix
	-n	number of nodes
	-t	number of threads
EOF

exit 0

fi

printf "\n\n++++ running: _submit_merfin.sh ++++\n\n"

#set options

while getopts ":a:v:d:r:p:q:o:n:t:" opt; do

	case $opt in
		a)
			ASM="$OPTARG"
			export ASM
			echo "Assembly: -a $OPTARG"
			;;
		v)
			VAR="$OPTARG"
			export VAR
        		echo "Variants: -v $OPTARG"
			;;
		d)
			ADB="$OPTARG"
			echo "Assembly meryl db: -d $OPTARG"
			;;
		r)
			RDB="$OPTARG"
			echo "Read meryl db: -r $OPTARG"
			;;
		p)
			PEAK="$OPTARG"
			echo "Haploid peak: -p $OPTARG"
			;;
		q)
			QUAL="$OPTARG"
			echo "QUAL filter: -q $OPTARG"
			;;
		o)
			OUT="$OPTARG"
			echo "Output prefix: -o $OPTARG"
			;;
		n)
			NODE="$OPTARG"
			echo "Number of nodes: -n $OPTARG"
			;;  
		t)
			CPU="$OPTARG"
			echo "Number of threads: -t $OPTARG"
			;;    
		\?)
			echo "ERROR - Invalid option: -$OPTARG" >&2
			exit 1
			;;		
	esac
	printf "\n"
done

RAM=$(mktemp -dt "k_eval.XXXXXXXX" --tmpdir=/run/user/$(id -u))
export RAM

printf "Temporary files written to: "
echo $RAM

printf "\nScript directory: "
echo ${BASH_SOURCE%/*}

bcftools index -f $VAR

bcftools index -s $VAR > ${RAM}/${VAR%.*}_chrs_stats.txt

n_vars=$(awk '{sum+=$3}END{print sum}' ${RAM}/${VAR%.*}_chrs_stats.txt)

awk -v n_vars=${n_vars} -v nodes=${NODE} 'BEGIN{pass=int(n_vars/nodes)}{if(sum>=pass){sum=0;printf "%s\n",chr;chr=""};sum+=$3;chr=chr","$1}END{printf "%s\n",chr}' ${RAM}/${VAR%.*}_chrs_stats.txt > ${RAM}/${VAR%.*}_chrs.txt

mkdir -p vcfs

cat ${RAM}/${VAR%.*}_chrs.txt | parallel -j ${CPU} "bcftools view -Ov ${VAR} -r {} > vcfs/{#}.vcf"

mkdir -p merfin
mkdir -p logs

for fl in $(ls vcfs);
do

	filename=$(basename -- "${fl}")
	out="${filename%.*}"
	
	log=logs/${VAR%.*}.%A.log
	sbatch --partition=vgl --nice=10000 --cpus-per-task=${CPU} --error=$log --output=$log $merfin/merfin.sh $ASM $ADB $RDB $PEAK "vcfs/${fl}" $out

done
