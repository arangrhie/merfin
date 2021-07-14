
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

  double            minMem, minMemTotal = 0;
  double            optMem, optMemTotal = 0;
  bool              useOpt = false;
  bool              useMin = false;

  //  Make readDB first so we know the k size
  merylFileReader*  readDB     = new merylFileReader(readDBname);
  
  //  Open sequence and build seqDBname if not provided
  load_Sequence();

  fprintf(stderr, "-- Estimating required space for loading '%s'\n", readDBname);
  readLookup = new merylExactLookup();
  readLookup->estimateMemoryUsage(readDB, maxMemory, minMem, optMem, minV, maxV);
  minMemTotal += minMem;
  optMemTotal += optMem;

  fprintf(stderr, "-- Estimating required space for loading '%s'\n", seqDBname);
  merylFileReader*  asmDB      = new merylFileReader(seqDBname);
  asmLookup  = new merylExactLookup();
  asmLookup->estimateMemoryUsage(asmDB, maxMemory, minMem, optMem, minV, maxV);
  minMemTotal += minMem;
  optMemTotal += optMem;

  if      (optMemTotal <= maxMemory)
    useOpt = true;
  else if (minMemTotal <= maxMemory)
    useMin = true;

  fprintf(stderr, "--\n");
  fprintf(stderr, "-- Minimal memory needed: %.3f GB%s\n", minMemTotal, (useMin) ? "  enabled" : "");
  fprintf(stderr, "-- Optimal memory needed: %.3f GB%s\n", optMemTotal, (useOpt) ? "  enabled" : "");
  fprintf(stderr, "-- Memory limit           %.3f GB\n",   maxMemory);
  fprintf(stderr, "--\n");

  if ((useMin == false) &&
      (useOpt == false)) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Not enough memory to load databases.  Increase -memory.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  fprintf(stderr, "-- Loading kmers from '%s' into lookup table.\n", readDBname);
  readLookup->load(readDB, maxMemory, useMin, useOpt, minV, maxV);

  fprintf(stderr, "-- Loading kmers from '%s' into lookup table.\n", seqDBname);
  asmLookup-> load(asmDB,  maxMemory, useMin, useOpt);

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

    snprintf(merylcount, FILENAME_MAX, "%s count k=%d memory=%.3f %s output %s", merylpath, kmer::merSize(), maxMemory, seqName, seqDBname);

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
