#!/bin/bash

#cut -f3-5 ${2} | paste ${1} - | awk '{if($3==0){a="NA"}else{a=$5};if($6==0){b="NA"}else{b=$8}; print a"\t"b}' | sort -nk1 -nk2 | uniq -c | sed "s/^[ ]*//" | tr -s " " | tr [:blank:] \\t

cut -f3-5 ${2} | paste ${1} - | awk '{if($3==0){a="NA"}else{a=$5};if($6==0){b="NA"}else{b=$8}; dups[a"\t"b]++} END{for (num in dups) {printf dups[num]"\t"num"\n"}}' | sort -nrk1 -nk2
