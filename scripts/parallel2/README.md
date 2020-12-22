# varMer, variant evaluation

Make sure your VCF is sorted. If not already sorted, you can sort with jvarkit:

```
cat $input.vcf | java -jar sortvcfonref2.jar -R $reference.fasta > $sorted.vcf
```

# General purpose oneliners
Sometimes you need to append something to the scaffold names in a fasta or vcf (e.g. to distinguish parental haplotypes in trio assemblies). This is a useful oneliner:
```
sed 's/tig[[:digit:]]*/&_mat/'
```