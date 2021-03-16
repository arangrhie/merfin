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
	-d	meryl db
	-k	kmer length
	-p	haploid peak
	-q	QUAL variant filter
	-c	max number of combinatory calls (default 15)
	-o	output prefix
	-t	number of threads
EOF

exit 0

fi

printf "\n\n++++ running: k_eval.sh ++++\n\n"

#set options

while getopts ":a:v:d:k:p:q:c:o:t:" opt; do

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
			DB="$OPTARG"
			echo "Meryl db: -d $OPTARG"
			;;
		k)
			K="$OPTARG"
			export K
			echo "Kmer length: -k $OPTARG"
			;;
		p)
			PEAK="$OPTARG"
			echo "Haploid peak: -p $OPTARG"
			;;
		q)
			QUAL="$OPTARG"
			echo "QUAL filter: -q $OPTARG"
			;;
		c)
			COMB="$OPTARG"
			echo "Combinations: -c $OPTARG"
			;;
		o)
			OUT="$OPTARG"
			echo "Output prefix: -o $OPTARG"
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

rm -f ${OUT}.lookup.table

if [[ -z ${COMB} ]]; then

	COMB=15

fi

filename=$(basename -- "${ASM}")
filename="${filename%.*}"

samtools faidx ${ASM}

if [[ ! -e ${filename}.meryl ]]; then

	meryl count k=${K} ${ASM} output ${filename}.meryl

fi

if [[ ! -e ${OUT}.asm.qv ]]; then

	meryl-lookup -dump -memory 200 -sequence ${ASM} -mers ${filename}.meryl | parallel -k --pipe -j ${CPU} -q awk '{print $1"\t"($NF+$(NF-2))}' > ${RAM}/${OUT}.asm
	meryl-lookup -dump -memory 200 -sequence ${ASM} -mers ${DB} | parallel -k --pipe -j ${CPU} -q awk '{COUNT=$NF+$(NF-2); print $1"\t"(COUNT) }' > ${RAM}/${OUT}.asm.read

	cut -f2 ${RAM}/${OUT}.asm.read | paste ${RAM}/${OUT}.asm - | parallel -k --pipe -j ${CPU} -q awk -v hap=${PEAK} '
	function round(x,   ival, aval, fraction)
	{
	   ival = int(x)    # integer part, int() truncates

	   # see if fractional part
	   if (ival == x)   # no fraction
		  return ival   # ensure no decimals

	   if (x < 0) {
		  aval = -x     # absolute value
		  ival = int(aval)
		  fraction = aval - ival
		  if (fraction >= .5)
			 return int(x) - 1   # -2.5 --> -3
		  else
			 return int(x)       # -2.3 --> -2
	   } else {
		  fraction = x - ival
		  if (fraction >= .5)
			 return ival + 1
		  else
			 return ival
	   }
	}	
	{
	as=$2;rd=$3/hap;
	if(rd==0/hap){rd_int=0}else if(rd>0/hap && rd<=1){rd_int=1}else{rd_int=rd};
	printf "%s\t%s\t%s\t",$0,rd,round(rd_int);
	if($3==0) {print "\tNA\tNA"}
	else if(rd>=as && as!=0){print (-1*(rd/as-1))}
	else if(rd>=as && as==0){print "NA\t"}
	else if(rd<as && as!=0){print (as/rd-1)}

	}' > ${RAM}/${OUT}.asm.combined
	
	Kasm=$(cat ${RAM}/${OUT}.asm.combined | parallel --pipe -j ${CPU} -q awk '{if($2>$5) sum+=1-$5/$2}END{print sum}' | awk '{sum+=$1}END{print sum}')
	Ktot=$(wc -l ${RAM}/${OUT}.asm.combined | awk '{print $1}')
	
	cp ${RAM}/${OUT}.asm.combined ${OUT}.asm.combined
	
	printf "%s\t%s\t%s\t" ${filename} $Kasm $Ktot > ${OUT}.asm.qv
	
	awk -v Kasm=${Kasm} -v Ktot=${Ktot} -v K=${K} 'BEGIN {printf -10*log(1-(1-Kasm/Ktot)^(1/K))/log(10)"\n"}' >> ${OUT}.asm.qv
	
fi

if [[ ! -e ${OUT}.lookup.table ]]; then

	bcftools view --threads $((${CPU}-1)) $VAR -i "QUAL>${QUAL} && (GT=\"0/1\" || GT=\"1/1\")" -Oz > ${RAM}/filtered.vcf.gz
	bcftools index -f ${RAM}/filtered.vcf.gz
	tabix -fp vcf ${RAM}/filtered.vcf.gz

	bcftools view --threads $((${CPU}-1)) -H ${RAM}/filtered.vcf.gz > ${RAM}/vars.vcf
	
	n_vars=$(wc -l ${RAM}/vars.vcf | awk '{print $1}')
	
	printf "\nNumber of variants considered: $n_vars\n"
	
	printf "\nComputing combinations:\n"

	a=1
	o=1
	u=1
	END=$(head -1 ${RAM}/vars.vcf | awk '{print $2}')
	CHR=$(head -1 ${RAM}/vars.vcf | awk '{print $1}')
	
	while read -r line
	
	do
	
	var=($line)
	
	STR=$(( ${var[1]} - ${K} + 1))
	
	if [[ $STR -le 1 ]]; then
	
		STR=1
		
	fi
		
		if ! [[ $END -gt $STR ]] || [[ $o -gt ${COMB} ]] || ! [[ ${var[0]} == $CHR ]]; then
	
			readarray -t var_ls < ${RAM}/var.vcf
	
			n=${#var_ls[@]}
			
		    chr_end=$(bcftools index -s ${RAM}/filtered.vcf.gz | grep "${CHR}" | awk '{print $2}')
	
			for (( e = 1; e < 2**n; e++ )); do

				combination=()
	
				idx=$e
	
				for (( j = 0; j < $n; j++ )); do
	
					if (( idx % 2 )); then combination=("${combination[@]}" "${var_ls[$j]}"); fi
	
					idx=$((idx>>1))
	
				done
		
				FIRST=$(echo ${combination[0]} | awk '{print $2}')
				LAST=$(echo ${combination[-1]} | awk '{print $2}')
				LAST_LN=$(echo ${combination[-1]} | awk '{print $4}')
		
				STR=$(( ${FIRST} - ${K} + 1))
				END=$(( ${LAST} + ${K} - 2 + ${#LAST_LN}))
				CHR=$(echo ${combination[0]} | awk '{print $1}')
				
				if [[ $END -gt $chr_end ]]; then
	
					END=$chr_end
		
				fi				
			
				printf "%s\t%s\t%s\t%s\t%s\t%s\n" ${u} ${a} ${CHR} ${STR} ${END} $(printf "%s\n"  "${combination[@]}" | awk '{printf "%s,",$1":"$2}' | sed 's/.$//' ) >> ${RAM}/${OUT}.lookup.table
				
				u=$((u+1))

			done
			
			o=1
			truncate -s 0 ${RAM}/var.vcf
	
		fi
		
		CHR=${var[0]}
		END=$(( ${var[1]} + ${K} - 2 + ${#var[3]}))

		printf "%s\n" "${line}" >> ${RAM}/var.vcf
		
		if [[ "$(( ${a} % 1000 ))" -eq 0 ]] || [[ $a == $n_vars ]]; then
	
			printf "$a processed (${u} combinations).\n"
		
		fi
		
		a=$((a+1))
		o=$((o+1))

	done<${RAM}/vars.vcf
	
	printf "\n"
	
	cp ${RAM}/${OUT}.lookup.table ${OUT}.lookup.table
	
else

	cp ${OUT}.lookup.table ${RAM}/${OUT}.lookup.table

fi

if [[ ! -e ${OUT}.results.vars ]]; then
	
	printf "\nGenerating REF/ALT kmers\n"
	
	parallel --record-env
	parallel -a ${RAM}/${OUT}.lookup.table -j ${CPU} -k --env _ --colsep '\t' "sh ../k_eval_kmers.sh {1} {3} {4} {5} {6}" >> ${RAM}/var.fa
	
	printf "\n"

	meryl-lookup -dump -memory 200 -sequence ${RAM}/var.fa -mers ${filename}.meryl | awk '{print $1"\t"($NF+$(NF-2))}' > ${RAM}/${OUT}.var
	meryl-lookup -dump -memory 200 -sequence ${RAM}/var.fa -mers ${DB} | awk '{COUNT=$NF+$(NF-2); print $1"\t"(COUNT) }' > ${RAM}/${OUT}.var.read
	
	cut -f2 ${RAM}/${OUT}.var.read | paste ${RAM}/${OUT}.var - | parallel -k --pipe -j ${CPU} -q awk -v hap=${PEAK} '{
	as=$2;rd=$3/hap;
	printf "%s\t%s\t%0.f\t",$0,rd,rd;
	if (index($1, "ref") != 0) {diff=-1} else {diff=1}; 
	if(rd==0) {print "\tNA\tNA"}
	else if(rd>=as && as!=0 && as+diff!=0){print (-1*(rd/as-1))"\t"(-1*(rd/(as+diff)-1))}
	else if(rd>=as && as!=0 && as+diff==0){print (-1*(rd/as-1))"\tNA"}
	else if(rd>=as && as==0 && as+diff!=0){print "NA\t"(-1*(rd/(as+diff)-1))}
	else if(rd<as && as!=0 && as+diff!=0){print (as/rd-1)"\t"((as+diff)/rd-1)}
	else if(rd<as && as!=0 && as+diff==0){print (as/rd-1)"\tNA"}
	else {printf "whooops"}

	}' > ${RAM}/${OUT}.var.combined

	awk '{if(NR%2==1) {split($0,var,"_"); if (header!=var[1]) {header=var[1];split(var[1],id,":");split(id[2],coords,"-");printf "%s\t%s\t%s\t%s\t%s\n",substr(header,2),substr(id[1],2),coords[1],coords[2],kmers; kmers=""}}else{kmers=kmers" "$1}}' ${RAM}/var.fa > ${OUT}.kmers.txt

	awk 'function abs(x) {return x<0 ? -x : x}BEGIN{type="ref"}
			{
				split($1,a,"#");split(a[2],b,":");split(b[2],c,"_");
				if(c[3]<=f)
				{
					if(type=="ref")
					{
						printf a[1]"\t"b[1]":"c[1]"\t"
						if(pre_NA==0){if(pre_sign>=0){var=1}else{var=-1};
						printf var*pre_abs_avg/f"\t0"}else{printf "NA\t"pre_NA};
						if(pst_NA==0){if(pst_sign>=0){var=1}else{var=-1};
						printf "\t"var*pst_abs_avg/f"\t0\t"}else{printf "\tNA\t"pst_NA"\t"}
					}
					 else
					{
						if(pre_NA==0){if(pre_sign>=0){var=1}else{var=-1};
						printf var*pre_abs_avg/f"\t0"}else{printf "NA\t"pre_NA};
						if(pst_NA==0){if(pst_sign>=0){var=1}else{var=-1};
						printf "\t"var*pst_abs_avg/f"\t0\n"}else{printf "\tNA\t"pst_NA"\n"}
					}
					pre_sign=0;
					pst_sign=0;
					pre_abs_avg=0;
					pst_abs_avg=0;
					pre_NA=0;
					pst_NA=0;				
				}
				if($6!="NA"){pre_sign+=$6;pre_abs_avg+=abs($6)}else{pre_NA+=1};
				if($7!="NA"){pst_sign+=$7;pst_abs_avg+=abs($7)}else{pst_NA+=1};
				type=c[2];
				f=c[3]
			}END{
					if(pre_NA==0){if(pre_sign>=0){var=1}else{var=-1};
					printf var*pre_abs_avg/f"\t0"}else{printf "NA\t"pre_NA};
					if(pst_NA==0){if(pst_sign>=0){var=1}else{var=-1};
					printf "\t"var*pst_abs_avg/f"\t0\n"}else{printf "\tNA\t"pst_NA"\n\n"}			
			}' ${RAM}/${OUT}.var.combined > ${RAM}/${OUT}.calls
			
			#awk 'function abs(x) {return x<0 ? -x : x} {if (abs((1-abs($8)))-abs(1-abs($9))>0 && abs((1-abs($10)))-abs(1-abs($11))>0) printf "%s\t%.5f\t%.5f\n", $0,abs((1-abs($8)))-abs(1-abs($9)),abs((1-abs($10)))-abs(1-abs($11))}' ram.calls > ram.diff

	cut -f1,2,6 ${RAM}/${OUT}.lookup.table | paste ${RAM}/${OUT}.calls - > ${OUT}.results
	
	awk 'function abs(x) {return x<0 ? -x : x}
			{
			 if(id==$12)
			 {
			 if(false<$10){next}
			 else if(false==$10 && abs($9)>abs(trend)){next}		 
			 }
			 else{print var}
			 trend=$9;false=$10;id=$12;var=$13
		}' ${OUT}.results > ${OUT}.results.vars

fi