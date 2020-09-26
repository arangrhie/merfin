#!/bin/bash

variantCaller ${1} -r ${2}.fasta -o vcfs/${3%.*}.vcf --referenceWindowsFile slices/${3} --algorithm=arrow --numWorkers ${4}