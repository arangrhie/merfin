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

#include "runtime.H"

#include "kmers.H"
#include "system.H"
#include "sequence.H"
#include "bits.H"
#include "strings.H"
#include "files.H"
#include "vcf.H"
#include "varMer.H"
#include "types.H"

#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

uint64 varMer::peak = 0;

double
varMer::getKmetric(merylExactLookup   *rlookup,
           merylExactLookup   *alookup,
           kmer                fmer,
           kmer                rmer,
           double             &readK,
           double             &asmK ) {

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

void
varMer::addSeqPath(string seq, vector<int> idxPath, vector<uint32> varIdxPath, vector<uint32> varLenPath) {

  vector<string>::iterator it = find (seqs.begin(), seqs.end(), seq);
  if ( it != seqs.end() ) { return; }

  // only insert elements if seq is a new sequence
  seqs.push_back(seq);
  gtPaths.push_back(idxPath);      // 0 = ref, 1 = alt1, 2 = alt2, ...
  idxPaths.push_back(varIdxPath);  // 0-base index where the var start is in the seq
  lenPaths.push_back(varLenPath);  // 0-base index where the var start is in the seq
  return;
}

void
varMer::score(merylExactLookup *rlookup, merylExactLookup *alookup) {

  //  iterate through each base and get kmer
  uint32 numM;  // num. missing kmers
  double minK;
  double maxK;
  string seq;
  double readK;
  double asmK;
  double kMetric;
  vector<double> m_ks;
  bool  evaluate = true;

  uint32 idx = 0;       // var index in the seqe

  //  get scores at each kmer pos and minimum read multiplicity
  for ( int ii = 0; ii < seqs.size(); ii++ ) {
    minK  = DBL_MAX;
    maxK  = -2;
    numM  = 0;

    evaluate = true;

    seq    = seqs.at(ii);
    m_ks.clear();
    // fprintf(stderr, "%s:%u-%u\t%s", posGt->_chr, posGt->_rStart, posGt->_rEnd, seq.c_str());

    idx = 0;

    kmerIterator kiter((char*) seq.c_str(), seq.size());
    while (kiter.nextBase()) {
      readK = 0;
      asmK  = 0;

      if (kiter.isValid()) {
        //  we only need readK and asmK, no need to get the kMetric here yet
        getKmetric(rlookup, alookup, kiter.fmer(), kiter.rmer(), readK, asmK);
      }

      //  is the idx anywhere close to idxPath?
      for ( int jj = 0; jj < idxPaths.at(ii).size(); jj++) {
        uint32 idxPath = idxPaths.at(ii).at(jj);
        uint32 lenPath = lenPaths.at(ii).at(jj);
        int    gtPath  = gtPaths.at(ii).at(jj);
        // fprintf(stderr, "\tidxPath:%u len:%u gt:%d", idxPath, lenPath, gtPath);
        if ( gtPath > 0 && idxPath + 1 - kmer::merSize() <= idx && idx < idxPath + lenPath + kmer::merSize()) {
          asmK++; // +1 as we are introducing a new kmer
          break;  // add only once
        }
      }

      //  re-define k* given rounded readK and asmK, in absolute values
      if (readK == 0) {
        kMetric = -1;  // use 0 if we are using non-abs k*
        numM++;

        // no need to evaluate min or max K
        evaluate = false;

      } else if (readK > asmK) {
        kMetric = readK / asmK - 1;

      } else {
        kMetric = asmK  / readK - 1;
      }

      if (evaluate && minK > kMetric) { minK = kMetric; };
      if (evaluate && maxK < kMetric) { maxK = kMetric; };

      m_ks.push_back(kMetric);
      // fprintf(stderr, "\tidx:%u Kr:%.3f Ka:%.0f K*:%.3f", idx, readK, asmK, kMetric);

      idx++;
    }

    numMs.push_back(numM);
    kstrs.push_back(m_ks);

    //  only evaluate when no missing kmers are generated
    if (evaluate) {
      minKs.insert(pair<double, int>(minK, ii));      // Automatically sorted by min value
      maxKs.insert(pair<double, int>(maxK, ii));      // Automatically sorted by min value
    }

    // fprintf(stderr, "\n");
  }
  return;
}


double
varMer::getMinAbsK(int idx) {

  // fprintf(stderr, "getMinAbsK(%d) called.\n", idx);
  double minAbsK = DBL_MAX;
  double absK;

  vector<double> kstr = kstrs.at(idx);
  for (int i = 0; i < kstr.size(); i++) {
    absK = kstr.at(i);
    if ( absK < 0 )  continue;  //  ignore the missings
    if ( absK < minAbsK ) { minAbsK = absK; }
  }
  //  if all kmers are 'missings'
  if (minAbsK == DBL_MAX) return -1;
  return minAbsK;
}


double
varMer::getMaxAbsK(int idx) {
  // fprintf(stderr, "getMaxAbsK(%d) called.\n", idx);

  double maxAbsK = -2;
  double absK;

  vector<double> kstr = kstrs.at(idx);
  for (int i = 0; i < kstr.size(); i++) {
    absK = kstr.at(i);
    if ( absK > maxAbsK ) { maxAbsK = absK; }
  }
  return maxAbsK;
}

double
varMer::getAvgAbsK(int idx) {
  // fprintf(stderr, "getAvgAbsK(%d) called.\n", idx);

  double sum = 0;
  double absK;

  vector<double> kstr = kstrs.at(idx);
  for (int i = 0; i < kstr.size(); i++) {
    absK = kstr.at(i);
    if (absK >= 0)
      sum += absK;
  }

  // no k* >= 0
  if ( kstr.size() == numMs.at(idx) )
    return -1;
  else
    return sum / ( kstr.size() - numMs.at(idx) );
}

double
varMer::getMedAbsK(int idx) {
  // fprintf(stderr, "getMedAbsK(%d) called.\n", idx);

  vector<double> kstr = kstrs.at(idx);
  sort (kstr.begin(), kstr.end());
  int i = 0;
  for (; i < kstr.size(); i++) {
    if ( kstr.at(i) >= 0 ) { break; }
  }

  //  no k* >= 0
  if (i == kstr.size())
     return -1;
  else
    return kstr.at(i + ((kstr.size() - i)/2));
}

