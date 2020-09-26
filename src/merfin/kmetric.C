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

#include "kmetric.H"

uint64 
peak = 0;
           
double
getKmetricDef(
  		   merylExactLookup   *rlookup,
           merylExactLookup   *alookup,
           kmer                fmer,
           kmer                rmer,
           vector<string>     &copyKmerDict,
           double             &readK,
           double             &asmK,
           double			  &prob) {

  uint64 fValue = 0;
  uint64 rValue = 0;
  double kMetric;

  rlookup->exists(fmer, fValue);
  rlookup->exists(rmer, rValue);

  readK = (double) (fValue + rValue) / peak;

  if (0 < readK && readK < 1) {
     readK = 1;
  } else {
     readK = round(readK);
  }

  fValue = 0;
  rValue = 0;
  alookup->exists(fmer, fValue);
  alookup->exists(rmer, rValue);

  asmK  = (double) (fValue + rValue);

  if ( readK == 0 ) {
    kMetric = 0;
  } else if ( asmK > readK ) {
    kMetric = asmK / readK - 1;
    kMetric = kMetric * -1;
  } else { // readK > asmK
    kMetric = readK / asmK - 1;
  }
  return kMetric;
}

double
getKmetricProb(
		   merylExactLookup   *rlookup,
           merylExactLookup   *alookup,
           kmer                fmer,
           kmer                rmer,
           vector<string>	  &copyKmerDict,
           double             &readK,
           double             &asmK,
           double			  &prob) {

  uint64 fValue = 0;
  uint64 rValue = 0;
  uint64 tValue = 0;
  double kMetric;

  rlookup->exists(fmer, fValue);
  rlookup->exists(rmer, rValue);
  
  tValue = fValue + rValue;

  if (0 < tValue && tValue < copyKmerDict.size()) {
  
    std::string s = copyKmerDict[tValue];
    std::string delimiter = ",";
    readK = (int) stod(s.substr(0, s.find(delimiter)));
    prob = (double) stod(s.erase(0, s.find(delimiter) + delimiter.length()));

  } else {
  
  	prob = (double) 1;
    readK = round((double) prob / peak);

  }

  fValue = 0;
  rValue = 0;
  
  alookup->exists(fmer, fValue);
  alookup->exists(rmer, rValue);
  
  tValue = fValue + rValue;

  asmK  = (double) tValue;
  
  if ( asmK >= readK ) {
    kMetric = ((asmK - readK) / asmK) * prob;
  } else { // readK > asmK
    kMetric = 1 - ((readK - asmK) / readK) * prob;
  }

  return kMetric;
           
}