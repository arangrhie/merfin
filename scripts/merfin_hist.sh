#!/bin/bash

if [ ! -z ${8} ]; then

arg="-lookup ${8}"

fi

echo "time merfin -hist \
-sequence ${1} \
-seqmers ${2} \
-readmers ${3} \
-peak ${4} \
-output ${5} \
-memory1 30 -memory2 ${6} \
-threads ${7} \
${arg}"

time merfin -hist \
-sequence ${1} \
-seqmers ${2} \
-readmers ${3} \
-peak ${4} \
-output ${5} \
-memory1 30 -memory2 ${6} \
-threads ${7} \
${arg}
