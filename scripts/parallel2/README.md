# varMer, variant evaluation

Make sure your VCF is sorted. If not already sorted, you can sort with jvarkit:

```
cat $input.vcf | java -jar sortvcfonref2.jar -R $reference.fasta > $sorted.vcf
```
