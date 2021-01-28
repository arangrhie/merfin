#!/bin/bash

ulimit -s unlimited

mkdir -p results

if [ ! -z ${8} ]; then

arg="-lookup ${8}"

fi

printf "merfin \
-sequence ${1} \
-seqmers ${2} \
-readmers ${3} \
-peak ${4} \
-vcf ${5} \
-output results/${6} \
-memory1 30 \
-memory2 150 \
-vmer \
-threads ${7} \
${arg} \
\n\n"

time merfin \
-sequence ${1} \
-seqmers ${2} \
-readmers ${3} \
-peak ${4} \
-vcf ${5} \
-output results/${6} \
-memory1 30 \
-memory2 150 \
-vmer \
-threads ${7} \
${arg}
