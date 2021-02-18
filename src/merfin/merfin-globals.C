/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "merfin-globals.H"
#include "strings.H"


//  Read probabilities lookup table for 1-4 copy kmers.
void
merfinGlobal::load_Kmetric(void) {

  if (pLookupTable == nullptr)    //  No input supplied, no lookup table to make.
    return;

  fprintf(stderr, "-- Loading copy-number lookup table '%s'.\n\n", pLookupTable);

  if (fileExists(pLookupTable) == false) {
    fprintf(stderr, "ERROR: Lookup table (-lookup) file '%s' doesn't exist!\n", pLookupTable);
    exit(1);
  }

  compressedFileReader  F(pLookupTable);

  uint32   lineMax = 0;
  uint32   lineLen = 0;
  uint32   lineNum = 0;
  char    *line    = nullptr;

  while (AS_UTL_readLine(line, lineLen, lineMax, F.file())) {
    splitToWords  S(line, splitLetter, ',');

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
merfinGlobal::getK(kmer     fmer,
                   kmer     rmer,
                   double  &readK,
                   double  &asmK,
                   double  &prob) {

  //  A first lookup to get the read values.

  uint64 value = readLookup->value(fmer) + readLookup->value(rmer);

  //  Compute the default readK and prob based on the lookup values.
  //
  //    kmer_value / peak == 0  ->  readK = 0
  //    kmer_value / peak  < 1  ->  readK = 1
  //    else                        readK = round(kmer_value / peak)

  readK = 0.0;
  prob  = 1.0;

  if      (value == 0)
    readK = 0;
  else if (value < peak)
    readK = 1;
  else
    readK = round(value / peak);

  //  But if there are pre-loaded probabilities, use those.

  if ((value > 0) &&
      (value <= copyKmerK.size())) {
    readK = copyKmerK[value-1];
    prob  = copyKmerP[value-1];
  }

  //  Another lookup to get the assembly values.

  asmK = asmLookup->value(fmer) + asmLookup->value(rmer);
}



void
merfinGlobal::load_Kmers(void) {

  fprintf(stderr, "-- Loading kmers from '%s' into lookup table.\n", readDBname);

  merylFileReader*  readDB     = new merylFileReader(readDBname);
  merylFileReader*  asmDB      = new merylFileReader(seqDBname);

  double            minMem, minMemTotal = 0;
  double            optMem, optMemTotal = 0;
  bool              useOpt = false;
  bool              useMin = false;

  readLookup = new merylExactLookup();
  readLookup->estimateMemoryUsage(readDB, maxMemory, minMem, optMem, minV, maxV);
  minMemTotal += minMem;
  optMemTotal += optMem;

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

  fprintf(stderr, "-- Loading kmers from '%s' into lookup table.\n", seqDBname);

  readLookup->load(readDB, maxMemory, useMin, useOpt, minV, maxV);
  asmLookup-> load(asmDB,  maxMemory, useMin, useOpt);

  delete readDB;    //  Not needed anymore.
  delete asmDB;
}




void
merfinGlobal::open_Inputs(void) {

  //  Open input sequences.

  if (seqName != nullptr) {
    fprintf(stderr, "-- Opening sequences in '%s'.\n", seqName);
    seqFile = new dnaSeqFile(seqName);
  }

  //  Open VCF input.  This is only needed for reportType OP_VAR_MER,
  //  but we don't check anything here.

  if (vcfName != nullptr) {
    fprintf(stderr, "-- Opening vcf file '%s'.\n", vcfName);
    inVcf = new vcfFile(vcfName);
  }

  //  Process the vcf.  Only needded for -vmer mode.

  if (inVcf) {
    fprintf(stderr, "Merge variants within %u-mer bases, splitting combinations greater than %u.\n",
            kmer::merSize(), comb);
    inVcf->mergeChrPosGT(kmer::merSize(), comb, nosplit);
  }
}
