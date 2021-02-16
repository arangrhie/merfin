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

double 
peak = 0;
           
void
getK(
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

  rlookup->exists(fmer, fValue);
  rlookup->exists(rmer, rValue);

  getreadK(fValue, rValue, copyKmerDict, readK, prob);

  fValue = 0;
  rValue = 0;
  alookup->exists(fmer, fValue);
  alookup->exists(rmer, rValue);

  asmK  = (double) (fValue + rValue);

}

double
getKmetric(double readK, double asmK) {
  double kMetric;
  if ( readK == 0 ) {
    kMetric = 0;
  } else if ( asmK > readK ) {
    kMetric = (asmK / readK - 1) * -1;
  } else { // readK > asmK
    kMetric = readK / asmK - 1;
  }
  return kMetric;
}

double
getreadKdef(
  		   uint64   		   fValue,
           uint64   		   rValue,
           vector<string>     &copyKmerDict,
           double             &readK,
           double			  &prob
		   ) {

  readK = (double) (fValue + rValue) / peak;

  if (0 < readK && readK < 1) {
     readK = 1;
  } else {
     readK = round(readK);
  }
  
  prob = (double) 1;
  return readK;
}

double
getreadKprob(
  		   uint64   		   fValue,
           uint64   		   rValue,
           vector<string>     &copyKmerDict,
           double             &readK,
           double			  &prob
		   ) {

  double tValue = fValue + rValue;

  if (0 < tValue && tValue <= copyKmerDict.size()) {
  
    std::string s = copyKmerDict[tValue-1];
    std::string delimiter = ",";
    readK = (int) stod(s.substr(0, s.find(delimiter)));
    prob = (double) stod(s.erase(0, s.find(delimiter) + delimiter.length()));

  } else {
  
  	prob = (double) 1;
    readK = round((double) tValue / peak);

  }

  return readK;
  
}
