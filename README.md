# Merfin

Evaluate genome assembly and variant calls with k*

## Installation

* Required: git v.2.12 or higher

```
git clone https://github.com/arangrhie/merfin.git
cd merfin/src
make -j 12
```

`meryl` and `meryl-utility` are installed as submodules.


## Running Merfin

Merfin is still under active development. Currently, `-hist` and `-dump` is feature complete, but will be updated in the near future.

```
cd ../build/bin/
./merfin

usage: ./merfin <report-type> \
         -sequence <seq.fasta>   \
         -seqmers  <seq.meryl>   \
         -readmers <read.meryl>    \
         -peak     <haploid_peak>  \
         -vcf      <input.vcf>     \
         -output   <output>        

  Predict the kmer from <input.vcf> given sequence <seq.fasta>
  and lookup the k-mer multiplicity from sequence and reads.

  Input -sequence and -vcf files can be FASTA or FASTQ; uncompressed, gz, bz2 or xz compressed

  Each input database can be filtered by value.  More advanced filtering
  requires a new database to be constructed using meryl.
    -min   m    Ignore kmers with value below m
    -max   m    Ignore kmers with value above m
    -threads t  Number of threads to use when constructing lookup table.

  Memory usage can be limited, within reason, by sacrificing kmer lookup
  speed.  If the lookup table requires more memory than allowed, the program
  exits with an error.
    -memory m   Don't use more than m GB memory

  Exactly one report type must be specified.


  -hist
   Generate a 0-centered k* histogram for sequences in <input.fasta>.
   Positive k* values are expected collapsed copies.
   Negative k* values are expected expanded  copies.
   Closer to 0 means the expected and found k-mers are well balenced, 1:1.
   Reports QV at the end, in stderr.
   Required: -sequence, -seqmers, -readmers, -peak, and -output.

   Output: k* <tab> frequency


  -dump
   Dump readK, asmK, and k* per bases (k-mers) in <input.fasta>.
   Required: -sequence, -seqmers, -readmers, -peak, and -output

   Output: seqName <tab> seqPos <tab> readK <tab> asmK <tab> k*
      seqName    - name of the sequence this kmer is from
      seqPos     - start position (0-based) of the kmer in the sequence
      readK      - normalized read copies (read multiplicity / peak)
      asmK       - assembly copies as found in <seq.meryl>
      k*         - 0-centered k* value

  -vmer (experimental)
   Score each variant, or variants within distance k and its combination by k*.
   Required: -sequence, -seqmers, -readmers, -peak, -vcf, and -output

   Output:
    <output>.debug : some useful info for debugging.
                     seqName <tab> varMerStart <tab> varMerEnd <tab> varMerSeq <tab> score <tab> path
      seqName         - name of the sequence this kmer is from
      varMerStart     - start position (0-based) of the variant (s), including sequences upstream of k-1 bp
      varMerEnd       - end position (1-based) of the variant (s), including sequences downstream of k-1 bp
      varMerSeq       - combination of variant sequence to evalute
      score           - score, min k*? median k*? to be decided...
      path            - position of the variant (type of variant used: 0 = hap1, 1 = hap2, 2 = ref)
    <output>.pri.vcf  - variants chosen. use bcftools to polish <seq.fata> (Not implemented yet)
    <output>.alt.vcf  - variants not chosen but have high support (Not implemented yet)
```

Example run:
```
// Get histogram and QV specific to chrX - QV is correct. Histogram will be biased for collapses
merfin -hist -memory 16 -sequence chrX.fasta -seqmers chrX.meryl -readmers chrX.read.meryl/ -peak 104 -output out.chrX.hist

// Get histogram and QV specific to chrX - Histogram is correct. QV will be biased
merfin -hist -memory 16 -sequence chrX.fasta -seqmers asm.meryl -readmers chrX.read.meryl/ -peak 104 -output out.chrX.hist

// Get histogram and QV for the full assembly - Histogram and QV correct
merfin -hist -memory 120 -sequence asm.fasta -seqmers asm.meryl -readmers read.meryl/ -peak 104 -output out.chrX.hist

// Dump readK, asmK, k* for each position of chrX
merfin -dump -memory 16 -sequence chrX.fasta -seqmers asm.meryl -readmers chrX.read.meryl/ -peak 104 -output out.dump.gz

// Score each variant call and sort
merfin -vmer -memory 16 -sequence chrX.fasta -seqmers asm.meryl -readmers chrX.read.meryl/ -peak 104 -vcf chrX.tiny.vcf -output out.dump.gz
```

## Acknowledgements
This code was developed as part of the T2T consortium chm13-polishing working group by the following individuals:
* Giulio Formenti
* Arang Rhie
* Brian Walenz
* Sergey Koren
* Adam Phillippy
