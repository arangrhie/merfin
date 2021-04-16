
/******************************************************************************
 *
 *  This is a k-mer based variant evaluation tool for polishing assemblies.
 *
 *  This software is based on:
 *    'Meryl'                  (https://github.com/marbl/meryl)
 *
 *  This is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "merfin-globals.H"

#include "sweatShop.H"


//  sweatShop sequence loading method.
//   - make a new merfinInput object.
//   - try to load data into it.  if there is no data left
//     in the iput, return nullptr to indicate all input
//     has been loaded.
//   - otherwise, do any initilization needed.
//   - return the new object; the sweatShop will then pass
//     the object to a compute function
//
void *
loadSequence(void *G) {
  merfinGlobal  *g = (merfinGlobal *)G;
  merfinInput   *s = new merfinInput;

  //  Try to load a new input sequence.  If it fails, destroy the merfinInput
  //  object we're trying to load and signal that we're done.

  if (g->seqFile->loadSequence(s->seq) == false) {
    delete s;
    return(nullptr);
  }

  //  Loaded something, initialize for processing.

  s->kiter.addSequence(s->seq.bases(), s->seq.length());

  //  We could, but do not, allocate space for undr and over here.  There's
  //  no gain to allocating here; we'd just be reserving lots of unused
  //  memory.

  //fprintf(stderr, "Loaded sequence %s\n", s->seq.ident());
  return(s);
}


//  sweatShop compute and output functions.  All are in other files, just to
//  organize things a bit.
//
void processHistogram(void *G, void *T, void *S);
void processDump     (void *G, void *T, void *S);
void processVariants (void *G, void *T, void *S);

void outputHistogram (void *G, void *S);
void outputDump      (void *G, void *S);
void outputVariants  (void *G, void *S);

void computeCompleteness(merfinGlobal *G);




int
main(int32 argc, char **argv) {
  merfinGlobal  *G = new merfinGlobal;

  std::vector<const char *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if        (strcmp(argv[arg], "-sequence") == 0) {
      G->seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-seqmers") == 0) {
      G->seqDBname = argv[++arg];

    } else if (strcmp(argv[arg], "-readmers") == 0) {
      G->readDBname = argv[++arg];

    } else if (strcmp(argv[arg], "-peak") == 0) {
      G->peak = strtodouble(argv[++arg]);

    } else if (strcmp(argv[arg], "-lookup") == 0) {
      G->pLookupTable = argv[++arg];

    } else if (strcmp(argv[arg], "-vcf") == 0) {
      G->vcfName = argv[++arg];

    } else if (strcmp(argv[arg], "-output") == 0) {
      G->outName = argv[++arg];

    } else if (strcmp(argv[arg], "-min") == 0) {
      G->minV = strtouint64(argv[++arg]);

    } else if (strcmp(argv[arg], "-max") == 0) {
      G->maxV = strtouint64(argv[++arg]);

    } else if (strcmp(argv[arg], "-threads") == 0) {
      G->threads = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-memory") == 0) {
      G->maxMemory = strtodouble(argv[++arg]);

#if 0
    } else if (strcmp(argv[arg], "-memory1") == 0) {
      G->maxMemory += strtodouble(argv[++arg]);

    } else if (strcmp(argv[arg], "-memory2") == 0) {
      G->maxMemory += strtodouble(argv[++arg]);
#endif

    } else if (strcmp(argv[arg], "-nosplit") == 0) {
      G->nosplit = true;

    } else if (strcmp(argv[arg], "-disable-kstar") == 0) {
      G->bykstar = false;

    } else if (strcmp(argv[arg], "-hist") == 0) {
      G->reportType = OP_HIST;

    } else if (strcmp(argv[arg], "-dump") == 0) {
      G->reportType = OP_DUMP;

    } else if (strcmp(argv[arg], "-skipMissing") == 0) {
      G->skipMissing = true;

    } else if (strcmp(argv[arg], "-vmer") == 0) {
      G->reportType = OP_VAR_MER;

    } else if (strcmp(argv[arg], "-completeness") == 0) {
      G->reportType = OP_COMPL;

    } else if (strcmp(argv[arg], "-comb") == 0) {
      G->comb = strtouint32(argv[++arg]);

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }
  }

  //  Check inputs are present for the various modes.

  if ((G->reportType == OP_HIST) ||
      (G->reportType == OP_DUMP) ||
      (G->reportType == OP_VAR_MER)) {
    if (G->seqName == nullptr)   err.push_back("No input sequences (-sequence) supplied.\n");
    if (G->outName == nullptr)   err.push_back("No output (-output) supplied.\n");
  }

  if  (G->reportType == OP_VAR_MER) {
    if (G->vcfName == nullptr)   err.push_back("No variant call input (-vcf) supplied; mandatory for -vmer reports.\n");
  }

  if  (G->reportType == OP_NONE) {
    err.push_back("No report type (-hist, -dump, -vmer, -completeness) supplied.\n");
  }

  if (G->seqDBname == nullptr)   err.push_back("No sequence meryl database (-seqmers) supplied.\n");
  if (G->readDBname == nullptr)  err.push_back("No read meryl database (-readmers) supplied.\n");
  if (G->peak == 0)              err.push_back("Peak=0 or no haploid peak (-peak) supplied.\n");

  //  check reportType OP_VAR_MER has vcfName

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s <report-type>            \\\n", argv[0]);
    fprintf(stderr, "         -sequence <seq.fasta>     \\\n");
    fprintf(stderr, "         -seqmers  <seq.meryl>     \\\n");
    fprintf(stderr, "         -readmers <read.meryl>    \\\n");
    fprintf(stderr, "         -peak     <haploid_peak>  \\\n");
    fprintf(stderr, "         -lookup   <lookup_table>  \\\n");
    fprintf(stderr, "         -vcf      <input.vcf>     \\\n");
    fprintf(stderr, "         -output   <output>        \n\n");
    fprintf(stderr, "  Predict the kmer consequences of variant calls <input.vcf> given the consensus sequence <seq.fasta>\n");
    fprintf(stderr, "  and lookup the k-mer multiplicity in the consensus sequence <seq.meryl> and in the reads <read.meryl>.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Input -sequence and -vcf files can be FASTA or FASTQ; uncompressed, gz, bz2 or xz compressed\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Each input database can be filtered by value.  More advanced filtering\n");
    fprintf(stderr, "  requires a new database to be constructed using meryl.\n");
    fprintf(stderr, "    -min   m    Ignore kmers with value below m\n");
    fprintf(stderr, "    -max   m    Ignore kmers with value above m\n");
    fprintf(stderr, "    -threads t  Multithreading for meryl lookup table construction, dump and hist.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Memory usage can be limited, within reason, by sacrificing kmer lookup\n");
    fprintf(stderr, "  speed.  If the lookup table requires more memory than allowed, the program\n");
    fprintf(stderr, "  exits with an error.\n");
    fprintf(stderr, "    -memory1 m   Don't use more than m GB memory for loading seqmers\n");
    fprintf(stderr, "    -memory2 m   Don't use more than m GB memory for loading readmers\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -lookup file   Optional input vector of probabilities.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Exactly one report type must be specified.\n");
    fprintf(stderr, "\n\n");
    fprintf(stderr, "  -hist\n");
    fprintf(stderr, "   Generate a 0-centered k* histogram for sequences in <input.fasta>.\n");
    fprintf(stderr, "   Positive k* values are expected collapsed copies.\n");
    fprintf(stderr, "   Negative k* values are expected expanded  copies.\n");
    fprintf(stderr, "   Closer to 0 means the expected and found k-mers are well balenced, 1:1.\n");
    fprintf(stderr, "   Reports QV at the end, in stderr.\n");
    fprintf(stderr, "   Required: -sequence, -seqmers, -readmers, -peak, and -output.\n");
    fprintf(stderr, "   Optional: -lookup <probabilities> use probabilities to adjust multiplicity to copy number\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   Output: k* <tab> frequency\n");
    fprintf(stderr, "\n\n");
    fprintf(stderr, "  -dump\n");
    fprintf(stderr, "   Dump readK, asmK, and k* per bases (k-mers) in <input.fasta>.\n");
    fprintf(stderr, "   Required: -sequence, -seqmers, -readmers, -peak, and -output\n");
    fprintf(stderr, "   Optional: -skipMissing will skip the missing kmer sites to be printed\n");
    fprintf(stderr, "             -lookup <probabilities> use probabilities to adjust multiplicity to copy number\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   Output: seqName <tab> seqPos <tab> readK <tab> asmK <tab> k*\n");
    fprintf(stderr, "      seqName    - name of the sequence this kmer is from\n");
    fprintf(stderr, "      seqPos     - start position (0-based) of the kmer in the sequence\n");
    fprintf(stderr, "      readK      - normalized read copies (read multiplicity / peak)\n");
    fprintf(stderr, "      asmK       - assembly copies as found in <seq.meryl>\n");
    fprintf(stderr, "      k*         - 0-centered k* value\n");
    fprintf(stderr, "\n\n");
    fprintf(stderr, "  -vmer\n");
    fprintf(stderr, "   Score each variant, or variants within distance k and their combinations by k*.\n");
    fprintf(stderr, "   Required: -sequence, -seqmers, -readmers, -peak, -vcf, and -output\n");
    fprintf(stderr, "   Optional: -comb <N>  set the max N of combinations of variants to be evaluated (default: 15)\n");
    fprintf(stderr, "             -nosplit   without this options combinations larger than N are split\n");
    fprintf(stderr, "             -disable-kstar  only missing kmers are considered for filtering.\n");
    //fprintf(stderr, "                        if chosen, use bcftools to compress and index, and consensus -H 1 -f <seq.fata> to polish.\n");
    //fprintf(stderr, "                        first ALT in heterozygous alleles are better supported by avg. |k*|.\n");
    fprintf(stderr, "             -lookup <probabilities> use probabilities to adjust multiplicity to copy number\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -completeness\n");
    fprintf(stderr, "   Compute kmer completeness using expected copy numbers for all kmers.\n");
    fprintf(stderr, "\n");    
    fprintf(stderr, "   Output files: <output>.debug and <output>.polish.vcf\n");
    fprintf(stderr, "    <output>.debug : some useful info for debugging.\n");
    fprintf(stderr, "                     seqName <tab> varMerStart <tab> varMerEnd <tab> varMerSeq <tab> score <tab> path\n");
    fprintf(stderr, "      varMerID                - unique numbering, starting from 0\n");
    fprintf(stderr, "      varMerRange             - seqName:start-end. position (0-based) of the variant (s), including sequences upstream and downstream of k-1 bp\n");
    fprintf(stderr, "      varMerSeq               - combination of variant sequence to evalute\n");
    fprintf(stderr, "      numMissings             - total number of missing kmers\n");
    fprintf(stderr, "      min k*                  - minimum of all |k*| for non-missing kmers. -1 when all kmers are missing.\n");
    fprintf(stderr, "      max k*                  - maximum of all |k*| for non-missing kmers. -1 when all kmers are missing.\n");
    fprintf(stderr, "      median k*               - median  of all |k*| for non-missing kmers. -1 when all kmers are missing.\n");
    fprintf(stderr, "      avg k*                  - average of all |k*| for non-missing kmers. -1 when all kmers are missing.\n");
    fprintf(stderr, "      avg ref-alt k*          - difference between reference and alternate average k*.\n");
    fprintf(stderr, "      delta kmer multiplicity - cumulative sum of kmer multiplicity variation. Positive values imply recovered kmers. Negative values overrepresented kmers introduced.\n");
    fprintf(stderr, "      record          - vcf record with <tab> replaced to <space>. only non-reference alleles are printed with GT being 1/1.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    <output>.polish.vcf : variants chosen.\n");
    //fprintf(stderr, "     use bcftools view -Oz <output>.polish.vcf and bcftools consensus -H 2 -f <seq.fata> to polish.\n");
    //fprintf(stderr, "     ALT alleles are favored with more support compared to the REF allele.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }


  omp_set_num_threads(G->threads);


  //  Open read kmers, build a lookup table.

  G->load_Kmetric();
  G->load_Kmers();
  G->open_Inputs();

  //  Configure the sweatShop.

  sweatShop       *ss = nullptr;
  merfinThrData   *td = new merfinThrData [G->threads];

  //  Check report type

  if (G->reportType == OP_HIST) {
    fprintf(stderr, "-- Generate histogram of the k* metric to '%s'.\n", G->outName);
    ss = new sweatShop(loadSequence, processHistogram, outputHistogram);
    ss->setInOrderOutput(false);
  }

  if (G->reportType == OP_DUMP) {
    fprintf(stderr, "-- Dump per-base k* metric to '%s'.\n", G->outName);
    ss = new sweatShop(loadSequence, processDump, outputDump);
    ss->setInOrderOutput(true);
  }

  if (G->reportType == OP_VAR_MER) {
    fprintf(stderr, "-- Generate variant mers and score them.\n");
    ss = new sweatShop(loadSequence, processVariants, outputVariants);
    ss->setInOrderOutput(false);
  }

  if (G->reportType == OP_COMPL) {
    fprintf(stderr, "-- Compute completeness.\n");
    computeCompleteness(G);
  }

  //  Run the compute.

  if (ss) {
    ss->setNumberOfWorkers(G->threads);

    for (uint32 tt=0; tt<G->threads; tt++) {
      td[tt].threadID = tt;
      ss->setThreadData(tt, td + tt);
    }

    ss->setLoaderBatchSize(1);
    ss->setLoaderQueueSize(G->threads * 2);
    ss->setWorkerBatchSize(1);
    ss->setWriterQueueSize(16384);

    ss->run(G, false);
  }

  //  Call all the report methods.  If there are no results for some report,
  //  the report isn't emitted.

  G->reportHistogram();
  G->reportCompleteness();

  //  Cleanup and quit.
  delete [] td;
  delete    ss;

  delete    G;

  fprintf(stderr, "Bye!\n");
  return(0);
}
