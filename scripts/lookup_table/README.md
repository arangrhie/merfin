# Generate lookup table

The script `lookup.R` is based on [Genomescope 2.0](http://qb.cshl.edu/genomescope/genomescope2.0/).

In addition to the canonical Genomescope output it generates fitted read multiplicity values and probabilities (`lookup_table.txt`).
These can be given to merfin with the `-lookup` option to generate more accurate K* estimates.

```
Rscript lookup.R <kmer_histogram> <k_size> <output_folder> <haploid=TRUE/FALSE>
```



