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
#include <list>

//  Add a new path if the variant 'seq' is not known already.
//
void
varMer::addSeqPath(string seq, vector<int> idxPath, vector<uint32> varIdxPath, vector<uint32> varLenPath) {
  if (find(seqs.begin(), seqs.end(), seq) == seqs.end()) {
    seqs.push_back(seq);
    gtPaths.push_back(idxPath);      // 0 = ref, 1 = alt1, 2 = alt2, ...
    idxPaths.push_back(varIdxPath);  // 0-base index where the var start is in the seq
    lenPaths.push_back(varLenPath);  // 0-base index where the var start is in the seq
  }
}


void
varMer::score(merfinGlobal *g) {

  //  iterate through each base and get kmer
  uint32 numM;  // num. missing kmers
  string seq;
  double prob;
  double readK;
  double asmK;
  double oDeltak;
  double nDeltak;
  double kMetric;
  vector<double> m_ks;
  vector<double> m_dks;

  uint32 idx = 0;       // var index in the seqe

  //  get scores at each kmer pos and minimum read multiplicity
  for ( int ii = 0; ii < seqs.size(); ii++ ) {
    numM  = 0;

    seq    = seqs.at(ii);
    m_ks.clear();
    m_dks.clear();
    // fprintf(stderr, "\n[ DEBUG ]:: score %d th combination :: %s:%u-%u\t%s\n", ii, posGt->_chr, posGt->_rStart, posGt->_rEnd, seq.c_str());

    idx = 0;

    kmerIterator kiter((char*) seq.c_str(), seq.size());
    while (kiter.nextBase()) {
      readK = 0;
      asmK  = 0;

      if (kiter.isValid()) {
        //  we only need readK and asmK, no need to get the kMetric here yet
        //fprintf(stderr, "[ DEBUG ] :: idx %u -- has a valid kmer. getKmetric()..\n", idx);
        g->getK(kiter.fmer(), kiter.rmer(), readK, asmK, prob);
        //fprintf(stderr, "[ DEBUG ] :: idx %u -- readK=%.0f , asmK=%.0f\n", idx, readK, asmK);
      }

      if (readK == 0) {
        numM++;
      }
    
      //  We only use the num missings in -filter mode
      if (g->reportType == OP_FILTER) {
        idx++;
        continue;
      }

      // store difference in kmer count accounting for uncertainty in the estimate of readK
      oDeltak = abs(readK - asmK) * prob;

      //  fprintf(stderr, "[ DEBUG ] :: is the idx a newly introduced kmer? Check the idx (%u) falls in any of the %lu idxPaths.at(%d)\n", idx, idxPaths.at(ii).size(), ii);
      //  fprintf(stderr, "[ DEBUG ] :: size of idxPaths.at(%d)=%lu | lenPaths.at(%d)=%lu | gtPath.at(%d)=%lu\n", ii, idxPaths.at(ii).size(), ii, lenPaths.at(ii).size(), ii, gtPaths.at(ii).size());
      for ( int jj = 0; jj < idxPaths.at(ii).size(); jj++) {
        uint32 idxPath = idxPaths.at(ii).at(jj);
        uint32 lenPath = lenPaths.at(ii).at(jj);
        int    gtPath  = gtPaths.at(ii).at(jj);
        //  fprintf(stderr, "\tidxPath:%u len:%u gt:%d", idxPath, lenPath, gtPath);
        if ( gtPath > 0 && idxPath + 1 - kmer::merSize() <= idx && idx < idxPath + lenPath + kmer::merSize()) {
          asmK++; // +1 as we are introducing a new kmer
          break;  // add only once
        }
      }
      //  fprintf(stderr, "\n");

      //  re-define k* given rounded readK and asmK, in absolute values
      if (readK == 0) {
        kMetric = -1;  // use 0 if we are using non-abs k*

      } else if (readK > asmK) {
        kMetric = readK / asmK - 1;

      } else {
        kMetric = asmK / readK - 1;
      }
      // recompute difference in kmer count if the variant is introduced
      nDeltak = abs(readK - asmK) * prob;

      // store new k*
      m_ks.push_back(kMetric);
      // fprintf(stderr, "[ DEBUG ] :: push_back kMetric for idx: %u. Kr: %.3f Ka: %.0f K*: %.3f\n\n", idx, readK, asmK, kMetric);
      // store the delta k* when the variant is applied
      m_dks.push_back(oDeltak-nDeltak);
      
      idx++;     
    }

    numMs.push_back(numM);
    kstrs.push_back(m_ks);
    dkstrs.push_back(m_dks);

    // avgKs.insert(pair<double, int>(getAvgAbsK(ii), ii));      // Automatically sorted by min value

    //  fprintf(stderr, "\n");
  }
}

/***
 * Get the best variant combination in -filter mode, but print the original vcf 
 ***/
vector<vcfRecord*>
varMer::bestFilter() {
  uint32 numMissing = UINT32_MAX;  //  actual minimum number of missing kmers in the combination with minimum missings
  vector<int> idxs;
  vector<vcfRecord*> records;      // empty vector with no records

  for ( int ii = 0; ii < numMs.size(); ii++ ) {
    //  ignore when all kmers are 'missings'
    if ( numMs.at(ii)  == seqs.at(ii).size() - kmer::merSize() + 1)  continue;

    // 1) path with no missings? add to report
    if ( numMs.at(ii) == 0 ) {
      idxs.push_back(ii);
      numMissing = 0;
    }

    // 2) lowest num. missings? add to report
    //    only if no paths satisfying 1) are found
    if ( numMs.at(ii) < numMissing ) {
      numMissing = numMs.at(ii);
      idxs.clear();
      idxs.push_back(ii);

    } else if ( numMs.at(ii) == numMissing ) {
      //  has the same numMissing
      idxs.push_back(ii);
    } // else : ignore

  }

  //  ignore if all kmers only increase missing kmers
  if ( idxs.size() == 0 ) return records;

  //  report all the rest
  list<int> gtIdxs;
  for ( int ii = 0; ii < idxs.size(); ii++ ) {
    int idx = idxs.at(ii);
    for ( int i = 0; i < gtPaths.at(idx).size(); i++) {
      if (gtPaths.at(idx).at(i) > 0) gtIdxs.push_back(i);
    }
  }

  gtIdxs.sort();
  gtIdxs.unique();

  for (list<int>::iterator it=gtIdxs.begin(); it!=gtIdxs.end(); ++it)
    records.push_back(posGt->_gts[*it]->_record);

  return records;
}

/***
 * Get the better variant, report HOM only if the variant reduces num. missing kmers
 ***/
string
varMer::betterVariant() {
  uint32 numMissing = UINT32_MAX;  //  actual minimum number of missing kmers in the combination with minimum missings
  vector<int> idxs;
  //  fprintf(stderr, "[ DEBUG ] :: Better mode called\n");

  //  add reference path as default
  if ( numMs.size() == 0 ) return "";
  uint32 refMissing = numMs.at(0);
  numMissing = refMissing;
  
  //  only add the variant if the numMs are less or equal to the ref path
  for ( int ii = 0; ii < numMs.size(); ii++ ) {

    // lowest num. missings? add to report
    if ( numMs.at(ii) < numMissing ) {
      numMissing = numMs.at(ii);
      idxs.clear();
      idxs.push_back(ii);

    } else if ( numMs.at(ii) == numMissing && numMs.at(ii) < refMissing ) {
      //  has the same numMissing, less than refMissing
      idxs.push_back(ii);

    } // else : ignore

  }

  //  ignore if all kmers only increase missing kmers
  if ( idxs.size() == 0 ) return "";

  int idx = idxs.at(0); // path index
  //  report only one best
  //  select the longest alt if there are multiple bests
  //
  //  1) only one combination has the minimum num. of missings
  if ( idxs.size() == 1 ) {
    // return records as HOM
    return getHomRecord(idx);
  }
  // 2) multiple combinations with equal num. of missings
  else {
    //  get the path idx with longest seq length
    uint32 seqLenMax = seqs.at(idx).size();
    for ( int ii = 1; ii < idxs.size(); ii++ ) {
      uint32 seqLen  = seqs.at(idxs.at(ii)).length();
      if ( seqLen > seqLenMax ) {
        seqLenMax = seqLen;
        idx = idxs.at(ii);
      }
    }
    // return records as HOM
    return getHomRecord(idx);
  }
}

string
varMer::strictPolish() {
  uint32 numMissing = UINT32_MAX;  //  actual minimum number of missing kmers in the combination with minimum missings
  vector<int> idxs;

  //  add reference path as default
  if ( numMs.size() == 0 ) return "";
  uint32 refMissing = numMs.at(0);
  numMissing = refMissing;
  
  //  only add the variant if the numMs are less or equal to the ref path
  for ( int ii = 0; ii < numMs.size(); ii++ ) {

    // lowest num. missings? add to report
    if ( numMs.at(ii) < numMissing ) {
      numMissing = numMs.at(ii);
      idxs.clear();
      idxs.push_back(ii);

    } else if ( numMs.at(ii) == numMissing && numMs.at(ii) < refMissing ) {
      //  has the same numMissing, less than refMissing
      idxs.push_back(ii);

    } // else : ignore

  }

  //  ignore if all kmers only increase missing kmers
  if ( idxs.size() == 0 ) return "";

  list<int> gtIdxs;
  int idx = idxs.at(0); // path index
  //  report only one best
  //  select the longest alt if there are multiple bests
  //
  //  1) only one combination has the minimum num. of missings
  if ( idxs.size() == 1 ) {
    // return records as HOM
    return getHomRecord(idx);
  }
  // 2) multiple combinations with equal num. of missings
  else {
    //  TODO: Adjust this part
    //  get the path idx with longest seq length
    uint32 seqLenMax = seqs.at(idx).size();
    for ( int ii = 1; ii < idxs.size(); ii++ ) {
      uint32 seqLen  = seqs.at(idxs.at(ii)).length();
      if ( seqLen > seqLenMax ) {
        seqLenMax = seqLen;
        idx = idxs.at(ii);
      }
    }
    // return records as HOM
    return getHomRecord(idx);
  }
}

string
varMer::loosePolish() {

  //  fprintf(stderr, "[ DEBUG ] :: Loose mode called\n");
  
  uint32 numMissing = UINT32_MAX;  //  actual minimum number of missing kmers in the combination with minimum missings
  vector<int> idxs;

  //  add reference path as default
  if ( numMs.size() == 0 ) return "";

  uint32 refMissing = numMs.at(0);
  numMissing = refMissing;
  
  //  only add the variant if the numMs are less or equal to the ref path
  for ( int ii = 0; ii < numMs.size(); ii++ ) {

    // lowest num. missings? add to report
    if ( numMs.at(ii) < numMissing ) {
      numMissing = numMs.at(ii);
      idxs.clear();
      idxs.push_back(ii);

    } else if ( numMs.at(ii) == numMissing && numMs.at(ii) <= refMissing ) {
      //  has the same numMissing, less than refMissing
      idxs.push_back(ii);

    } // else : ignore

  }
  
  //  fprintf(stderr, "[ DEBUG ] :: num. combinations picked up: %lu at Mmin %u\n", idxs.size(), numMissing);

  //  ignore if all kmers only increase missing kmers
  if ( idxs.size() == 0 ) return "";

  int idx = idxs.at(0); // path index
  //  report all, with warnings for multiple pathes
  
  //  1) only one combination has the minimum num. of missings
  if ( idxs.size() == 1 ) {
    // return records as HOM
    return getHomRecord(idx);
  }
  // 2) multiple combinations with equal num. of missings
  else {

    //  only 2 pathes exist, first is the REF path
    if (idxs.at(0) == 0 && idxs.size() == 2)
      return getHomRecord(idxs.at(1));
    
    //  find the path with most ALT variants
    int maxVars = 0;
    int maxIdx  = idx;

    //  for each pathes
    for ( int ii = 1; ii < idxs.size(); ii++) {
      int count = 0;
      idx = idxs.at(ii);

      // check if each gt is ALT
      for ( uint64 i = 0; i < gtPaths.at(idx).size(); i++) {
        int altIdx = gtPaths[idx][i];
        if (altIdx > 0) {
          count++;
        }
      }

      if (count > maxVars) {
        maxVars = count;
        maxIdx  = idx;
      }
    }
    
    fprintf(stderr, "[ WARNING ] :: Multiple (%lu) alternate pathes detected in a path beginning with variant : %s", idxs.size(), posGt->_gts[0]->_record->save().c_str());
    fprintf(stderr, "[ WARNING ] :: Max. %d ALT variants selected\n", maxVars);
    return getHomRecord(maxIdx);
  }
}

/***
 * -polish mode. Use k*
 */
string
varMer::bestVariant() {

  uint32 numMissing = UINT32_MAX;  //  actual minimum number of missing kmers in the combination with minimum missings
  vector<int> idxs;

  for ( int ii = 0; ii < numMs.size(); ii++ ) {
    //  ignore when all kmers are 'missings'
    if ( numMs.at(ii)  == seqs.at(ii).size() - kmer::merSize() + 1)  continue;

    //  found a smaller numMissing
    if ( numMs.at(ii) < numMissing ) {
      numMissing = numMs.at(ii);
      idxs.clear();
      idxs.push_back(ii);

    } else if ( numMs.at(ii) == numMissing ) {
      //  has the same numMissing
      idxs.push_back(ii);
      
    } // else : ignore
  }

  //  ignore if all kmers creates only missings
  if ( numMissing == UINT32_MAX ) return "";

  //  only one combination has the minimum num. of missings
  if ( idxs.size() == 1 ) {
    // get idx of the combination
    int idx = idxs.at(0);

    // make a vcf record and return
    return getHomRecord(idx);

  } else if ( idxs.size() > 1) {
    // found multiple combination equally having the minimum missing kmers.
    // sort by the Avg. Absolute K* closest to 0
    // and the second best as het
    for ( int i = 0; i < idxs.size(); i++ ) {
      int idx = idxs.at(i);
      //avgKs.insert(pair<double, int>(getAvgAbsK(idx), idx));
      avgKs.insert(make_pair(getTotdK(idx), idx));
    }

    multimap<double, int>::iterator it = avgKs.begin();
    double  avgK1 = (*it).first;
    int      idx1 = (*it).second;
    it++;
    double avgK2 = (*it).first;
    int     idx2 = (*it).second;

    double kDiff = avgK1 - avgK2;

    if ( avgK1 == avgK2 ) {
      // if equaly scored, chose the longer allele as hap1
      if ( seqs.at(idx1).length() >= seqs.at(idx2).length() ) {
        return getHetRecord(idx1, idx2);
      } else {
        return getHetRecord(idx2, idx1);
      }
    } else {
      return getHomRecord(idx1);
    }
  }

  // no idxs?
  return "";
}

/***
 * Experimental
 */
string
varMer::getHetRecord(int idx1, int idx2) {
  string records;

  for ( int i = 0; i < gtPaths.at(idx1).size(); i++) {
    int altIdx1 = gtPaths.at(idx1).at(i);
    int altIdx2 = gtPaths.at(idx2).at(i);

    // alt1 == ref && alt2 == ref: ignore
    if ( altIdx1 + altIdx2 > 0 ) {
      //  Two problems here.  First, this isn't really the correct qual
      //  value.  Second, qual can be a floating point value and so we input
      //  it as such, but printing it as one results in 6 digits of precision
      //  with no easy way to change it.  So we cast it to an integer.
      string qualStr = to_string((int)posGt->_gts[i]->_qual);

      records = records +
                posGt->_chr + "\t" + 
                to_string(posGt->_gts[i]->_pos+1) + "\t" + 
                "." + "\t" +
                posGt->_gts[i]->_alleles[0] + "\t";

      // alt1 == alt2: 1/1
      // Sometimes, there are cases where path is different but the allele chosen overlaps
      if ( altIdx1 == altIdx2 ) {
        records = records +
                  posGt->_gts[i]->_alleles[altIdx1] + "\t" +
                  qualStr + "\t" +
                  "PASS\t.\tGT\t1/1\n";


      }
      // alt1 == ref && alt2 == alt: 0/1
      else if ( altIdx1 == 0 && altIdx2 > 0 ) {
        records = records +
                  posGt->_gts[i]->_alleles[altIdx2] + "\t" +
                  qualStr + "\t" +
                  "PASS\t.\tGT\t0/1\n";

      }
      // alt1 == alt && alt2 == alt: 1/2
      else if ( altIdx1 > 0 && altIdx2 > 0 ) {
        records = records +
                  posGt->_gts[i]->_alleles[altIdx1] + "," + posGt->_gts[i]->_alleles[altIdx2] + "\t" +
                  qualStr + "\t" +
                  "PASS\t.\tGT\t1/2\n";
      }
      // alt1 == alt && alt2 == ref: 1/0
      else if ( altIdx1 > 0 && altIdx2 == 0 ) {
        records = records +
                  posGt->_gts[i]->_alleles[altIdx1] + "\t" +
                  qualStr + "\t" +
                  "PASS\t.\tGT\t1/0\n";
      }
    }
  }
  return records;
}

string
varMer::getHomRecord(int idx) {
  string records;
  for ( int i = 0; i < gtPaths.at(idx).size(); i++) {
    int altIdx = gtPaths.at(idx).at(i);
    if ( altIdx > 0 ) {
      string qualStr = to_string((int)posGt->_gts[i]->_qual);  //  See comment above.

      // altIdx is the reference allele: ignore
      records = records +
                posGt->_chr + "\t" +
                to_string(posGt->_gts[i]->_pos+1) + "\t.\t" +
                posGt->_gts[i]->_alleles[0] + "\t" +
                posGt->_gts[i]->_alleles[altIdx] + "\t" +
                qualStr + "\t" +
                "PASS\t.\tGT\t1/1\n";
    }
  }
  return records;
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

double
varMer::getAvgAbsdK(int idx, double kstr_ref) {
  // fprintf(stderr, "getAvgAbsdK(%d) called.\n", idx);

  double sum = 0;
  double absK;

  vector<double> kstr = kstrs.at(idx);
  for (int i = 0; i < kstr.size(); i++) {
    absK = kstr.at(i);
    if (absK >= 0)
      sum += absK;
  }

  // no absk* >= 0
  if ( kstr.size() == numMs.at(idx) )
    return -1;
  else
    return (sum / ( kstr.size() - numMs.at(idx)) - kstr_ref);
}

double
varMer::getTotdK(int idx) {
  // fprintf(stderr, "getTotdK(%d) called.\n", idx);

  double sum = 0;

  vector<double> dkstr = dkstrs.at(idx);
  for (int i = 0; i < dkstr.size(); i++) {
    sum += dkstr.at(i);
  }

  return sum;
}
