
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

#ifndef AS_UTL_VARMER_H
#define AS_UTL_VARMER_H

#include "runtime.H"
#include "types.H"
#include "kmers.H"
#include "varMer.H"
#include <string>
#include <algorithm>
#include <map>
#include <vector>
#include <cmath>
using namespace std;


class varMer {

public:

  varMer(posGT* posGt){
    this->posGt = posGt;
  };
  ~varMer();

  static double getKmetric(merylExactLookup *rlookup, merylExactLookup *alookup,
                           kmer fmer, kmer rmer, double &readK, double  &asmK);

  void addSeqPath(string seq, vector<int> idxPath, vector<uint32> refIdxList, vector<uint32> refLenPath);

  void score(merylExactLookup *rlookup, merylExactLookup *alookup);

  string  bestVariant();
  string  getHomRecord(int idx);
  string  getHetRecord(int idx1, int idx2);
  double  getMinAbsK(int idx);
  double  getMaxAbsK(int idx);
  double  getAvgAbsK(int idx);  // avg. k*
  double  getMedAbsK(int idx);  // median k*

  static void   setPeak(uint64 p = 0) { peak   = p;  };
  static uint64 getPeak()         { return peak; };

public:
  vector< vector<int> >     gtPaths;  //  contains idx
  vector< vector<uint32> >  idxPaths;  //  contains idx
  vector< vector<uint32> >  lenPaths;  //  contains idx
  vector< string >          seqs;   //  sequence combinations
  vector< uint32 >          numMs;  //  num missing kmers. follows the order of seqs
  vector< vector<double> >  kstrs;  //  k* values, kstrs.at(0) is the k* value of the first seq
  multimap<double, int >    avgKs;  //  avg |k*|,  sorted from lowest to greatest, with idx on the second
  posGT*               posGt;
  vector< uint32 >     order;

private:
  static uint64        peak;
};

#endif  //  AS_UTL_VARMER_H
