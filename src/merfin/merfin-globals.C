
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
#include "strings.H"
#include <libgen.h>

//  Read probabilities lookup table for 1-4 copy kmers.
void
merfinGlobal::load_Kmetric(void) {

  if (pLookupTable == nullptr)    //  No input supplied, no lookup table to make.
    return;

  fprintf(stderr, "-- Loading probability table '%s'.\n\n", pLookupTable);

  if (fileExists(pLookupTable) == false) {
    fprintf(stderr, "ERROR: Probability table (-prob) file '%s' doesn't exist!\n", pLookupTable);
    exit(1);
  }

  compressedFileReader  F(pLookupTable);

  uint32   lineMax = 0;
  uint32   lineLen = 0;
  uint32   lineNum = 0;
  char    *line    = nullptr;

  while (AS_UTL_readLine(line, lineLen, lineMax, F.file())) {
    splitToWords  S(line, ',');

    if (S.numWords() == 2) {
      uint32  k = S.touint32(0);
      double  p = S.todouble(1);

      copyKmerK.push_back(k);
      copyKmerP.push_back(p);

      lineNum++;

      fprintf(stderr, "Copy-number: %u\t\tReadK: %u\tProbability: %f\n", lineNum, k, p);
    }

    else {
      fprintf(stderr, "Copy-number: invalid line %u:  '%s'\n", lineNum, line);
    }
  }

  delete [] line;
}



void
merfinGlobal::getK(kmvalu   seqValue,   //  Value of the kmer in the read database
                   kmvalu   asmValue,   //  Value of the kmer in the assembly database
                   double  &readK,
                   double  &asmK,
                   double  &prob) {

  //  Compute the default readK and prob based on the lookup values.
  //    kmer_value / peak == 0  ->  readK = 0
  //    kmer_value / peak  < 1  ->  readK = 1
  //    else                        readK = round(kmer_value / peak)
  //
  //  asmK is always just the asmValue.

  readK = 0.0;
  asmK  = asmValue;
  prob  = 1.0;

  if      (seqValue == 0)
    readK = 0;
  else if (seqValue < peak)
    readK = 1;
  else
    readK = round(seqValue / peak);

  //  If there are pre-loaded probabilities, use those.

  if ((seqValue > 0) &&
      (seqValue <= copyKmerK.size())) {
    readK = copyKmerK[seqValue-1];
    prob  = copyKmerP[seqValue-1];
  }
}


void
merfinGlobal::getK(kmer     fmer,
                   kmer     rmer,
                   double  &readK,
                   double  &asmK,
                   double  &prob) {
  getK(readLookup->value(fmer) + readLookup->value(rmer),
       asmLookup->value(fmer)  + asmLookup->value(rmer),
       readK, asmK, prob);
}



void
merfinGlobal::load_Kmers(void) {
  double            reqMemory = 0.0;

  //  Make readDB first so we know the k size
  merylFileReader*  readDB     = new merylFileReader(readDBname);
  
  //  Open sequence and build seqDBname if not provided
  load_Sequence();

  //  Since estimateMemoryUsage() is now including space for temporary
  //  buffers that are used only when loading, this estimate is significantly
  //  too large for small datasets.  If table1 and table2 need only 5 GB
  //  memory (each), the estimate for each will also include several GB for
  //  buffers (based on the number of threads); 16 threads = 8 GB buffers.
  //  So while the data needs 10 GB memory, meryl claims it needs 2x 13 GB =
  //  26 GB memory.  Since the tables are loaded sequentially, it really only
  //  needs 13 - 8 + 13 - 8 = 18 GB peak, 10 GB final.
#warning estimate is too high

  fprintf(stderr, "-- Estimating required space for loading '%s'\n", readDBname);
  readLookup = new merylExactLookup();
  reqMemory += readLookup->estimateMemoryUsage(readDB, maxMemory, 0, minV, maxV);

  fprintf(stderr, "-- Estimating required space for loading '%s'\n", seqDBname);
  merylFileReader*  asmDB      = new merylFileReader(seqDBname);
  asmLookup  = new merylExactLookup();
  reqMemory += asmLookup->estimateMemoryUsage(asmDB, maxMemory, 0);

  fprintf(stderr, "--\n");
  fprintf(stderr, "-- Memory needed: %.3f GB\n", reqMemory);
  fprintf(stderr, "-- Memory limit:  %.3f GB\n", maxMemory);
  fprintf(stderr, "--\n");

  if (reqMemory > maxMemory) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Not enough memory to load databases.  Increase -memory.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  fprintf(stderr, "-- Loading kmers from '%s' into lookup table.\n", readDBname);
  readLookup->load(readDB, maxMemory, 0, minV, maxV);

  fprintf(stderr, "-- Loading kmers from '%s' into lookup table.\n", seqDBname);
  asmLookup-> load(asmDB,  maxMemory, 0);

  delete readDB;    //  Not needed anymore.
  delete asmDB;
}

void
merfinGlobal::load_Sequence(void) {

  if (seqDBname == nullptr) {
    seqDBname = new char[FILENAME_MAX+1];
    snprintf(seqDBname, FILENAME_MAX, "%s.meryl", basename(seqName));
    fprintf(stderr, "-- No -seqmer given. Build sequence db as '%s'.\n", seqDBname);

    //  TODO: Replace this hacky counting when meryl has proper APIs
    char merylcount[FILENAME_MAX+1];
    char merylpath[FILENAME_MAX+1];

    if (strcmp(execName, "merfin") == 0)
      snprintf(merylpath, FILENAME_MAX, "meryl");
    else
      snprintf(merylpath, FILENAME_MAX, "%s/meryl", dirname(execName));

    snprintf(merylcount, FILENAME_MAX, "%s threads=%d count k=%d memory=%.3f %s output %s",
             merylpath, threads, kmer::merSize(), maxMemory, seqName, seqDBname);

    fprintf(stderr, "%s\n\n", merylcount);
    system(merylcount);
    fprintf(stderr, "\n");

  }

  //  Open input sequence.
  if (seqName != nullptr) {
    fprintf(stderr, "-- Opening sequences in '%s'.\n", seqName);
    seqFile = new dnaSeqFile(seqName);
  } 

}



void
merfinGlobal::open_Inputs(void) {

  if (reportType == OP_COMPL) {
    return;
  }

  //  Open VCF input.  This is only needed for reportType OP_VAR_MER.
  //  main() checks that we have a vcfName.

  if (reportType == OP_FILTER || reportType == OP_POLISH) {
    fprintf(stderr, "-- Opening vcf file '%s'.\n", vcfName);
    inVcf = new vcfFile(vcfName);

    fprintf(stderr, "Merge variants within %u-mer bases, splitting combinations greater than %u.\n",
            kmer::merSize(), comb);
    inVcf->mergeChrPosGT(kmer::merSize(), comb, nosplit);
  }
}
