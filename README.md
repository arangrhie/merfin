# Merfin

k-mer-based assembly and variant calling evaluation for improved consensus accuracy.

## Installation

* Required: git v.2.12 or higher, OMP (only for parallelization)

```
git clone https://github.com/arangrhie/merfin.git
cd merfin/src
make -j 12
```

`meryl` and `meryl-utility` are installed as submodules.

## Running Merfin

Merfin can be used to assess collapsed or duplicated region of the assembly (`-hist`, `-dump`) or to evaluate variant calls (`-vmer`). QV estimates for all scaffolds will also be generated with `-hist` and `-dump`.

In all cases a haploid/diploid peak estimate must be provided (`-peak`), either from the kmer histogram, or computed using the GenomeScope 2.0 model available under `scripts/genomescope`.

As a rule of thumb, the `-peak` should be:
- haploid, if the reference used for read mapping and variant calling contains both the primary and the haplotigs, or both haplotypes of a trio
- diploid (i.e. twice the haploid peak), for haploid representations of diploid genomes

Optionally, a custom table of probabilities can be used as input (`-lookup`), also generated using the script under `scripts/genomescope`. This is still experimental.

### Assess collapses/duplications and per base QV ###

The output of `-dump` can be further converted to `.Wig/.bw` tracks for visualization on IGV with:

```
awk 'BEGIN{print "track autoScale=on"}{if($1!=chr){chr=$1; print "variableStep chrom="chr" span=1"};if($3!=0){printf $2+1"\t"$5"\n"}}' $dump_output > $dump_output.Wig
#Convert to bigWig:
wigToBigWig $dump_output.Wig $dump_output.bw
```

### Filter variant calls for polishing ###

The input `.vcf` can be supplied with the `-vmer` option. The ouput will include only the variants that passed the Merfin screening.

Once the filtered `.vcf` is generated, the assembly can be polished with:

```
bcftools -Oz merfin_output.polish.vcf > merfin_output.polish.vcf.gz #bgzip merfin output
bcftools index merfin_output.polish.vcf.gz
bcftools consensus merfin_output.polish.vcf.gz -f assembly.fasta -H 2 > polished_assembly.fasta # -H 2 applies only alt alleles at each position
```

Two set of similar scripts for further parallelization on HPC (slurm) are available under `scripts/parallel1` and `scripts/parallel2`.

Merfin is still under active development. Feel free to reach out to us if you have any question.

### Helper ###

```
cd ../build/bin/
./merfin

usage: ./merfin <report-type>      \
         -sequence <seq.fasta>     \
         -seqmers  <seq.meryl>     \
         -readmers <read.meryl>    \
         -peak     <haploid_peak>  \
         -lookup   <lookup_table>  \
         -vcf      <input.vcf>     \
         -output   <output>           

  Predict the kmer consequences of variant calls <input.vcf> given the consensus sequence <seq.fasta>
  and lookup the k-mer multiplicity in the consensus sequence <seq.meryl> and in the reads <read.meryl>.

  Input -sequence and -vcf files can be FASTA or FASTQ; uncompressed, gz, bz2 or xz compressed

  Each input database can be filtered by value.  More advanced filtering
  requires a new database to be constructed using meryl.
    -min   m    Ignore kmers with value below m
    -max   m    Ignore kmers with value above m
    -threads t  Multithreading for meryl lookup table construction, dump and hist.

  Memory usage can be limited, within reason, by sacrificing kmer lookup
  speed.  If the lookup table requires more memory than allowed, the program
  exits with an error.
  -memory1 m   Don't use more than m GB memory for loading seqmers
  -memory2 m   Don't use more than m GB memory for loading readmers
    
  -lookup optional input vector of probabilities.

  Exactly one report type must be specified.

  -hist
   Generate a 0-centered k* histogram for sequences in <input.fasta>.
   Positive k* values are expected collapsed copies.
   Negative k* values are expected expanded  copies.
   Closer to 0 means the expected and found k-mers are well balenced, 1:1.
   Reports QV at the end, in stderr.
   Required: -sequence, -seqmers, -readmers, -peak, and -output.
   Optional: -lookup <probabilities> use probabilities to adjust multiplicity to copy number

   Output: k* <tab> frequency


  -dump
   Dump readK, asmK, and k* per bases (k-mers) in <input.fasta>.
   Required: -sequence, -seqmers, -readmers, -peak, and -output
   Optional: -skipMissing will skip the missing kmer sites to be printed
             -lookup <probabilities> use probabilities to adjust multiplicity to copy number

   Output: seqName <tab> seqPos <tab> readK <tab> asmK <tab> k*
      seqName    - name of the sequence this kmer is from
      seqPos     - start position (0-based) of the kmer in the sequence
      readK      - normalized read copies (read multiplicity / peak)
      asmK       - assembly copies as found in <seq.meryl>
      k*         - 0-centered k* value


  -vmer
   Score each variant, or variants within distance k and their combinations by k*.
   Required: -sequence, -seqmers, -readmers, -peak, -vcf, and -output
   Optional: -comb <N>  set the max N of combinations of variants to be evaluated (default: 15)
             -nosplit   without this options combinations larger than N are split
             -by-kstar  output variants by kstar. *experimental*
                        if chosen, use bcftools to compress and index, and consensus -H 1 -f <seq.fata> to polish.
                        first ALT in heterozygous alleles are better supported by avg. |k*|.
             -lookup <probabilities> use probabilities to adjust multiplicity to copy number

   Output files: <output>.debug and <output>.polish.vcf
    <output>.debug : some useful info for debugging.
                     seqName <tab> varMerStart <tab> varMerEnd <tab> varMerSeq <tab> score <tab> path
      varMerID    - unique numbering, starting from 0
      varMerRange - seqName:start-end. position (0-based) of the variant (s), including sequences upstream and downstream of k-1 bp
      varMerSeq   - combination of variant sequence to evalute
      numMissings - total number of missing kmers
      min k*      - minimum of all |k*| for non-missing kmers. -1 when all kmers are missing.
      max k*      - maximum of all |k*| for non-missing kmers. -1 when all kmers are missing.
      avg k*      - average of all |k*| for non-missing kmers. -1 when all kmers are missing.
      median k*   - median  of all |k*| for non-missing kmers. -1 when all kmers are missing.
      record      - vcf record with <tab> replaced to <space>. only non-reference alleles are printed with GT being 1/1.

    <output>.polish.vcf : variants chosen.
     use bcftools view -Oz <output>.polish.vcf and bcftools consensus -H 2 -f <seq.fata> to polish.
     ALT alleles are favored with more support compared to the REF allele.
```

Example run:
```
// Get histogram and QV specific to chrX - QV is correct. Histogram will be biased for collapses
merfin -hist -memory1 2 -memory2 16 -sequence chrX.fasta -seqmers chrX.meryl -readmers chrX.read.meryl/ -peak 104 -output out.chrX.hist

// Get histogram and QV specific to chrX - Histogram is correct. QV will be biased
merfin -hist -memory1 2 -memory2 16 -sequence chrX.fasta -seqmers asm.meryl -readmers chrX.read.meryl/ -peak 104 -output out.chrX.hist

// Get histogram and QV for the full assembly - Histogram and QV correct
merfin -hist -memory1 12 -memory2 120 -sequence asm.fasta -seqmers asm.meryl -readmers read.meryl/ -peak 104 -output out.chrX.hist

// Dump readK, asmK, k* for each position of chrX
merfin -dump -memory1 2 -memory2 16 -sequence chrX.fasta -seqmers asm.meryl -readmers chrX.read.meryl/ -peak 104 -output out.dump.gz

// Score each variant call and sort
merfin -vmer -memory1 2 -memory2 16 -sequence chrX.fasta -seqmers asm.meryl -readmers chrX.read.meryl/ -peak 104 -vcf chrX.tiny.vcf -output out.dump.gz
```

## Acknowledgements
This code was developed as part of the T2T consortium chm13-polishing working group by the following individuals:
* Giulio Formenti
* Arang Rhie
* Brian Walenz
* Sergey Koren
* Adam Phillippy
