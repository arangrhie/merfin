# Merfin

Improved variant filtering and polishing via k-mer validation

## Installation

* Required: git >=v.2.12, gcc >=v7.1 with OMP (only for parallelization)

```
git clone https://github.com/arangrhie/merfin.git
cd merfin/src
make -j 12
```

During installation, `meryl` and `meryl-utility` are cloned as submodules, and `meryl`, `merfin` will be installed under `merfin/build/bin/`.

* Recommended: [GenomeScope 2.0](https://github.com/tbenavi1/genomescope2.0) for polishing and assembly evaluation

## Running Merfin

A more detailed WIKI with examples can be found [here](https://github.com/arangrhie/merfin/wiki/Best-practices-for-Merfin). Please read it through before running Merfin.


Merfin can be used to:
* filter any variant calls for accurate genotyping (`-filter`: reference is from another individual, i.e. GRCh38)
* assess collapsed or duplicated region of the assembly (`-hist` or `-dump`)
* QV* for all scaffolds and the assembly (`-hist`)
* K* completeness (`-completeness`)
* filter variant calls for polishing (`-polish`: reference is from the same individual)


### Determine kmer copy numbers

Except for the `-filter` mode, a haploid/diploid peak estimate must be provided (`-peak`).


This can be obtained from the kmer histogram or computed using tools such as [GenomeScope 2.0](https://github.com/tbenavi1/genomescope2.0) (kcov value).


As a rule of thumb, the `-peak` should be chosen from:
- haploid peak: partially-phased (assembly has both primary and alternate haplotigs) or fully-phased assemblies (i.e. trio-binning, trio-hifiasm, ...).
- diploid (i.e. twice the haploid peak): pseudo-haploid assembly (haploid representation of a diploid genome)


In addition to the `-peak`, we recommend to provide a lookup table of kmer multiplicities with fitted copy numbers and probabilities (`-prob`).

The lookup table can be generated with `--fitted_hist` option we added to [GenomeScope 2.0](https://github.com/tbenavi1/genomescope2.0). The probability is estimated for 0 to 4-copy kmer multiplicity and is prioritized over the estimates from `-peak`. Providing the fitted probability significantly improves the accuracy of all analyses.


Once downloaded and installed, run GenomeScope using the following options:

```
Rscript genomescope.R -i <kmer_histogram>  -k <k_size> -o <output_folder> -p <ploidy> --fitted_hist [ploidy] [verbose]
kmer_histogram  tab-delimited, 2-column file with (same as for Genomescope2, usually generated by meryl hist, or jellyfish)
k_size          kmer length used for the histogram
ploidy          haploid/diploid (default = 2)
--fitted_hist     generates fitted hist plot and lookup_table.txt
```

### Assess collapses and duplications

The output of `-dump` can be further converted to `.wig/.bw` tracks for visualization on IGV or UCSC Genome Browser with:

```shell
awk 'BEGIN{print "track autoScale=on"}{if($1!=chr){chr=$1; print "variableStep chrom="chr" span=1"};if($3!=0){printf $2+1"\t"$5"\n"}}' $dump_output > $dump_output.wig
# Generate chromosome sizes:
awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' genome.fasta > chr.sizes
# Convert to bigWig:
wigToBigWig $dump_output.Wig chr.sizes $dump_output.bw
```


#### Assess correlation between data types
For multiple data types (e.g. HiFi and Illumina), the genome-wide K* agreement can be assessed by generating a cartesian plot.
Use the same output of `-dump` with the scripts under `scripts/cartesian_plot`.


### Assess per base QV
Merfin quickly produces [Merqury](https://github.com/marbl/merqury) QV estimates for each scaffold and genome-wide averages when `-hist` is used. Merqury QV estimate consider only kmers missing from the read sets. In addition, Merfin produces a QV* estimate, which accounts also for kmers that are seen in excess with respect to their expected multiplicity predicted from the reads.

These analyses can be further refined when the lookup table is provided (`-prob`),  Missing kmers then include plausible low frequency kmers. 0 to 4-copy kmer multiplicity estimates are weighted for the probability that the multiplicity estimate was correct. 

### Filter variant calls for polishing

The input `.vcf` can be supplied with the `-polish` option. The includes variants that passed the Merfin screening.

Once the filtered `.vcf` is generated, the assembly can be polished with:

```
bcftools view -Oz $merfin_output.polish.vcf > $merfin_output.polish.vcf.gz #bgzip merfin output
bcftools index $merfin_output.polish.vcf.gz
bcftools consensus $merfin_output.polish.vcf.gz -f assembly.fasta -H 1 > polished_assembly.fasta # -H 1 applies only first allele from GT at each position
```

Merfin is still under active development. Feel free to reach out to us if you have any question.


### Helper

```
cd ../build/bin/
./merfin

usage: merfin <report-type>        \
         -sequence <seq.fasta>     \
         -readmers <read.meryl>    \
         -peak     <haploid_peak>  \
         -prob     <lookup_table>  \
         -vcf      <input.vcf>     \
         -output   <output>        

  Predict the kmer consequences of variant calls <input.vcf> given the consensus sequence <seq.fasta>
  and lookup the k-mer multiplicity in the consensus sequence <seq.meryl> and in the reads <read.meryl>.

  Input -sequence and -vcf files can be FASTA or FASTQ; uncompressed, gz, bz2 or xz compressed

  Each readmers can be filtered by value.  More advanced filtering
  requires a new database to be constructed using meryl.
    -min     m     Ignore kmers with value below m
    -max     m     Ignore kmers with value above m
    -threads t     Multithreading for meryl lookup table construction, dump and hist.

  Memory usage can be limited, within reason, by sacrificing kmer lookup
  speed.  If the lookup table requires more memory than allowed, the program
  exits with an error.
    -memory  m     Don't use more than m GB memory for loading mers

  For k* based evaluation and polishing, -peak is required with optional -prob.
    -peak    m     Required input to hard set copy 1 and infer multiplicity to copy number (recommended).
    -prob    file  Optional input vector of probabilities. Adjust multiplicity to copy number
                   in case both -prob and -peak are provided, -prob takes higher priority
                   than -peak for multiplicity listed in the vector table.

  By default, <seq.fasta>.meryl will be generated unless -seqmers is provided.
    -seqmers seq.meryl  Optional input for pre-built sequence meryl db

  Exactly one report type must be specified.


  -filter
   Filter variants within distance k and their combinations by missing k-mers.
   Assumes the reference (-sequence) is from a different individual.
   Required: -sequence, -readmers, -vcf, and -output
   Optional: -comb <N>  set the max N of combinations of variants to be evaluated (default: 15)
             -nosplit   without this options combinations larger than N are split
             -debug     output a debug log, into <output>.THREAD_ID.debug.gz

   Output: <output>.filter.vcf : variants chosen.


  -polish
   Score each variant, or variants within distance k and their combinations by k*.
   Assumes the reference (-sequence) is from the same individual.

   Required: -sequence, -readmers, -peak, -vcf, and -output
   Optional: -comb <N>    set the max N of combinations of variants to be evaluated (default: 15)
             -nosplit     without this options combinations larger than N are split
             -prob <file> use probabilities to adjust multiplicity to copy number (recommended)
             -debug       output a debug log, into <output>.THREAD_ID.debug.gz

   Output: <output>.polish.vcf : variants chosen.
     use bcftools view -Oz <output>.polish.vcf and bcftools consensus -H 1 -f <seq.fata> to polish.
     first ALT in heterozygous alleles are usually better supported by avg. |k*|.


  -hist
   Generate a 0-centered k* histogram for sequences in <input.fasta>.
     Positive k* values are expected collapsed copies.
     Negative k* values are expected expanded  copies.
     Closer to 0 means the expected and found k-mers are well balenced, 1:1.

   Required: -sequence, -readmers, -peak, and -output.
   Optional: -prob <file>  use probabilities to adjust multiplicity to copy number (recommended)

   Output: k* <tab> frequency
           Reports QV at the end, in stderr.


  -dump
   Dump readK, asmK, and k* per bases (k-mers) in <input.fasta>.

   Required: -sequence, -readmers, -peak, and -output
   Optional: -skipMissing  skip the missing kmer sites to be printed
             -prob <file>  use probabilities to adjust multiplicity to copy number (recommended)

   Output: seqName <tab> seqPos <tab> readK <tab> asmK <tab> k*
      seqName    - name of the sequence this kmer is from
      seqPos     - start position (0-based) of the kmer in the sequence
      readK      - normalized read copies (read multiplicity / peak)
      asmK       - assembly copies as found in <seq.meryl>
      k*         - 0-centered k* value


  -completeness
   Compute kmer completeness using expected copy numbers for all kmers.

   Required: -seqmers (or -sequence), -readmers, -peak
   Optional: -prob <file>  use probabilities to adjust multiplicity to copy number (recommended)

   Output: total kmers in reads, number of kmers under the expected copy number, and completeness


  Optional output from -debug in -filter and -polish:
   <output>.THREAD_ID.debug.gz : some useful info for debugging.
      seqName <tab> varMerStart <tab> varMerEnd <tab> varMerSeq <tab> score <tab> path
      varMerID                - unique numbering, starting from 0
      varMerRange             - seqName:start-end. position (0-based) of the variant (s),
                                including sequences upstream and downstream of k-1 bp
      varMerSeq               - combination of variant sequence to evalute
      numMissings             - total number of missing kmers
      min k*                  - minimum of all |k*| for non-missing kmers. -1 when all kmers are missing.
      max k*                  - maximum of all |k*| for non-missing kmers. -1 when all kmers are missing.
      median k*               - median  of all |k*| for non-missing kmers. -1 when all kmers are missing.
      avg k*                  - average of all |k*| for non-missing kmers. -1 when all kmers are missing.
      avg ref-alt k*          - difference between reference and alternate average k*.
      delta kmer multiplicity - cumulative sum of kmer multiplicity variation.
                                positive values imply recovered kmers, while
                                negative values imply overrepresented kmers introduced.
      record                  - vcf record with <tab> replaced to <space>.
                                only non-reference alleles are printed with GT being 1/1.
```


## Acknowledgements
This code was developed as part of the T2T consortium chm13-polishing working group by the following individuals:
* Giulio Formenti
* Arang Rhie
* Brian Walenz
* Sergey Koren
* Adam Phillippy

With special thanks for their support to integrate `--fitted_hist` option in GenomeScope 2.0:
* Rhyker Ranallo-Benavidez
* Michael Schatz
