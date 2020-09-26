# Merfin

k-mer-based assembly and variant calling evaluation for improved consensus accuracy.

## Installation

* Required: git v.2.12 or higher, OMP (for parallelization)

```
git clone https://github.com/arangrhie/merfin.git
cd merfin/src
make -j 12
```

`meryl` and `meryl-utility` are installed as submodules.


## Running Merfin

Merfin can be used to assess collapsed or duplicated region of the assembly (-hist, -dump) or to evaluate variant calls (-vmer). QV estimates for all scaffolds will also be generated with -hist and -dump.

In all cases a haploid peak estimate must be provided (-peak), either from the kmer histogram, or computed using the GenomeScope 2.0 model available under `scripts/genomescope`. 

Optionally, a custom table of probabilities can be used as input (-lookup), also generated using the script under `scripts/genomescope`.

Two set of similar scripts are available to run Merfin on a slurm cluster under `scripts/parallel1` and `scripts/parallel2`.

Merfin is still under active development. Feel free to reach out to us if you have any question.

```
cd ../build/bin/
./merfin

usage: ./merfin <report-type> \
         -sequence <seq.fasta>   \
         -seqmers  <seq.meryl>   \
         -readmers <read.meryl>    \
         -peak     <haploid_peak>  \
		 -lookup     <lookup_table> \
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
   Optional: -lookup

   Output: k* <tab> frequency

  -dump
   Dump readK, asmK, and k* per bases (k-mers) in <input.fasta>.
   Required: -sequence, -seqmers, -readmers, -peak, and -output
   Optional: -lookup

   Output: seqName <tab> seqPos <tab> readK <tab> asmK <tab> k*
      seqName    - name of the sequence this kmer is from
      seqPos     - start position (0-based) of the kmer in the sequence
      readK      - normalized read copies (read multiplicity / peak)
      asmK       - assembly copies as found in <seq.meryl>
      k*         - 0-centered k* value

  -vmer (experimental)
   Score each variant, or variants within distance k and its combination by k*.
   Required: -sequence, -seqmers, -readmers, -peak, -vcf, and -output
   Optional: -lookup

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
