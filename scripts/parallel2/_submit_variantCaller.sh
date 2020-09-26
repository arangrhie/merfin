#!/bin/bash

set -e -o pipefail

#++++                  This script is part of:                    ++++
#++++                      the VGP project                        ++++
#++++     Credit: Giulio Formenti gformenti@rockefeller.edu       ++++

#!/bin/bash

set -e

if [ -z $1 ]; then

	echo "use $0 -h for help"
	exit 0
elif [ $1 == "-h" ]; then

	cat << EOF

	Usage: '$0 -i species_ID -a assembly -b bam_folder -t threads'

	_submit_variantCaller.sh is used to polish a genome assembly generated using long reads.

	Required arguments are:
	-i the output prefix
	-a the genome assembly .fasta file
	-b the folder where raw data bam files are
	-n the number of nodes
	-t the number of threads

EOF

exit 0

fi

#set options

printf "\n"

while getopts ":i:a:b:n:t:" opt; do

	case $opt in
		i)
			ID=$OPTARG
			echo "Species ID: -i $OPTARG"
			;;
		a)
			ASM=$OPTARG
			echo "Genome assembly: -a $OPTARG"
			;;
        b)
        	BAM=$OPTARG
        	echo "Bam folder: -b $OPTARG"
            ;;
        n)
        	N_NODE=$OPTARG
        	echo "Nodes: -n $OPTARG"
            ;;
        t)
        	N_PROC=$OPTARG
        	echo "Threads: -t $OPTARG"
            ;;
		\?)
			echo "ERROR - Invalid option: -$OPTARG" >&2
			exit 1
			;;
	esac

printf "\n"

done

printf "\n"

RAM=$(mktemp -dt "k_eval.XXXXXXXX" --tmpdir=/run/user/$(id -u))
export RAM

printf "Temporary files written to: "
echo $RAM

printf "\nScript directory: "
echo ${BASH_SOURCE%/*}

if ! [[ -e "${ASM}.mmi" ]]; then

	echo "Generating index."

	pbmm2 index ${ASM} ${ASM}.mmi

	echo "Index ${ASM}.mmi generated."

fi

if ! [[ -e "${ID}_read_set.xml" ]]; then

	echo "Gathering bam files."

	dataset create --type SubreadSet --name ${ID} ${ID}_read_set.xml ${BAM}/*.subreads.bam

	echo "Bam files in ${ID}_read_set.xml"

fi

if ! [[ -e "aligned_reads.bam" ]]; then

	echo "Aligning..."

	pbmm2 align ${ASM}.mmi ${ID}_read_set.xml aligned_reads.bam -j ${N_PROC} --sort -m 10G
	
	echo "Generated aligned_reads.bam file."

fi

if ! [[ -e "aligned_reads.bam.pbi" ]]; then

	echo "Indexing aligned_reads.bam file."

	pbindex aligned_reads.bam

	echo "aligned_reads.bam file indexed."	

fi


if ! [[ -e "${ASM}.fai" ]]; then

	echo "Generating .fai index."	

	samtools faidx ${ASM}

	echo "${ASM}.fai file generated."	

fi

if ! [[ -e "polished_${ASM}" ]]; then
	
	if [[ $N_NODE == 1 ]]; then
	
		echo "Generating consensus..."
	
		variantCaller aligned_reads.bam -r ${ASM} -o polished_${ASM%.*}.fasta -o polished_${ASM%.*}.vcf --algorithm=arrow --numWorkers ${N_PROC}
		
		echo "Polished consensus written to polished_${ASM}."
		
	else

		n_bases=$(awk '{sum+=$2}END{print sum}' ${ASM}.fai)

		mkdir -p slices

		awk -v n_bases=${n_bases} -v nodes=${N_NODE} 'BEGIN{pass=int(n_bases/nodes);counter=1}{if(sum>=pass){sum=0;printf "%s\n",chr > "slices/"counter".txt";chr=""; counter+=1};sum+=$2;if(chr!=""){chr=chr"\n"$1":1-"$2}else{chr=$1":1-"$2}}END{printf "%s\n",chr > "slices/"counter".txt"}' ${ASM}.fai

		mkdir -p logs
		mkdir -p vcfs
		
		rm -f freebayes_jid

		for fl in $(ls slices)
		do

			filename=$(basename -- "${fl}")
			out="${filename%.*}"
	
			log=logs/${ASM%.*}.${fl%.*}.%A.log
			sbatch --partition=vgl --nice=10000 --cpus-per-task=${N_PROC} --error=$log --output=$log $merfin/variantCaller.sh aligned_reads.bam ${ASM%.*} ${fl} ${N_PROC} | awk '{print $4}' >> freebayes_jid

		done
		
		wait_for="--dependency=afterok:`awk '{printf $1","}' freebayes_jid | sed 's/.$//'`"
		
		cpus=2
		name=$1.consensus
		script=$merfin/variantCaller_cns.sh
		args=${ASM}
		walltime=2-0
		log=logs/$name.%A_%a.log

		echo "\
		sbatch --partition=vgl -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --time=$walltime --error=$log --output=$log $script $args"
		sbatch --partition=vgl -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --time=$walltime --error=$log --output=$log $script $args
	
	fi	

fi

