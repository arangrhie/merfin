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
#include <sstream>

uint64 varMer::peak = 0;

double
varMer::getKmetric(
		   merylExactLookup    *rlookup,
           merylExactLookup    *alookup,
           kmer                fmer,
           kmer                rmer,
           map<int, string>    pValuesDict,
           double              &readK,
           double              &asmK,
           double 			   &pValue
    		) {

  uint64 fValue = 0;
  uint64 rValue = 0;
  uint64 tValue = 0;
  double kMetric;

  rlookup->exists(fmer, fValue);
  rlookup->exists(rmer, rValue);
  
  tValue = fValue + rValue;
  
  fprintf(stderr, "tValue: '%lu'\n", tValue);

  if (0 < tValue && tValue < pValuesDict.size()) {

    fprintf(stderr, "pValuesDict.at: '%s'\n", pValuesDict.at(tValue).c_str());
  
    std::string s = pValuesDict.at(tValue).c_str();
    std::string delimiter = ",";
    readK = stod(s.substr(0, s.find(delimiter)));
    pValue = (double) stod(s.erase(0, s.find(delimiter) + delimiter.length()));

  } else {
  
  	pValue = (double) 1;
    readK = round((double) tValue / peak);

  }

	fprintf(stderr, "pValue: '%f'\n", pValue);
	fprintf(stderr, "ReadK: '%f'\n", readK);

  fValue = 0;
  rValue = 0;
  
  alookup->exists(fmer, fValue);
  alookup->exists(rmer, rValue);
  
  tValue = fValue + rValue;

  asmK  = (double) tValue;
  
  fprintf(stderr, "AsmK: '%f'\n", asmK);
  
  if ( asmK >= readK ) {
    kMetric = ((asmK - readK) / asmK) * pValue;
  } else { // readK > asmK
    kMetric = 1 - ((readK - asmK) / readK) * pValue;
  }
  fprintf(stderr, "kMetric: '%f'\n", kMetric);

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
varMer::score(
			  merylExactLookup *rlookup,
			  merylExactLookup *alookup,
			  map<int, string> pValuesDict
			  ) {

  //  iterate through each base and get kmer
  string seq;
  double readK;
  double asmK;
  double pValue;
  double kMetric;
  vector<double> m_ks;

  uint32 idx = 0;       // var index in the seqe

  //  get scores at each kmer pos and minimum read multiplicity
  for ( int ii = 0; ii < seqs.size(); ii++ ) {

    seq = seqs.at(ii);
    m_ks.clear();
    // fprintf(stderr, "%s:%u-%u\t%s", posGt->_chr, posGt->_rStart, posGt->_rEnd, seq.c_str());

    idx = 0;

    kmerIterator kiter((char*) seq.c_str(), seq.size());
    while (kiter.nextBase()) {
      readK = 0;
      asmK  = 0;
      pValue = 0;

      if (kiter.isValid()) {
        //  we only need readK and asmK, no need to get the kMetric here yet
        getKmetric(rlookup, alookup, kiter.fmer(), kiter.rmer(), pValuesDict, readK, asmK, pValue);
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

	  if ( asmK > readK ) {
		kMetric = ((asmK - readK) / asmK) * pValue;
	  } else { // readK >= asmK
		kMetric = ((readK - asmK) / readK) * pValue;
	  }

      m_ks.push_back(kMetric);
      // fprintf(stderr, "\tidx:%u Kr:%.3f Ka:%.0f K*:%.3f", idx, readK, asmK, kMetric);

      idx++;
    }

    kstrs.push_back(m_ks);

    // avgKs.insert(pair<double, int>(getAvgAbsK(ii), ii));      // Automatically sorted by min value

    // fprintf(stderr, "\n");
  }
  return;
}


string
varMer::bestVariant() {

    for ( int i = 0; i < kstrs.size(); i++ ) {

      avgKs.insert(pair<double, int>(getAvgAbsK(i), i));
      
    }
    
    multimap<double, int >::iterator it = avgKs.begin();
    double  avgK1 = (*it).first;
    int      idx1 = (*it).second;
    it++;
    double avgK2 = (*it).first;
    int     idx2 = (*it).second; 

    if ( avgK1 == avgK2 ) {
       
      // if equaly scored, choose the longer allele as hap1
      if ( seqs.at(idx1).length() >= seqs.at(idx2).length() ) {
        return getHetRecord(idx1, idx2);
      } else {
        return getHetRecord(idx2, idx1);
      }
    } else {

      return getHomRecord(idx1);
    }

  // no idxs?
  return "";
}

string
varMer::getHetRecord(int idx1, int idx2) {

  string records;
  for ( int i = 0; i < gtPaths.at(idx1).size(); i++) {
    int altIdx1 = gtPaths.at(idx1).at(i);
    int altIdx2 = gtPaths.at(idx2).at(i);

    // alt1 == ref && alt2 == ref: ignore
    if ( altIdx1 + altIdx2 > 0 ) {
      records = records + posGt->_chr + "\t" + 
                to_string(posGt->_gts->at(i)->_pos+1) + "\t.\t" +
                posGt->_gts->at(i)->alleles->at(0) + "\t";

      // alt1 == alt2: 1/1
      // Sometimes, there are cases where path is different but the allele chosen overlaps
      if ( altIdx1 == altIdx2 ) {
        records = records + posGt->_gts->at(i)->alleles->at(altIdx1) + "\t.\t" +
                  "PASS\t.\tGT\t1/1\n";


      }
      // alt1 == ref && alt2 == alt: 0/1
      else if ( altIdx1 == 0 && altIdx2 > 0 ) {
        records = records + posGt->_gts->at(i)->alleles->at(altIdx2) + "\t.\t" +
                  "PASS\t.\tGT\t0/1\n";

      } 
      // alt1 == alt && alt2 == alt: 1/2
      else if ( altIdx1 > 0 && altIdx2 > 0 ) {
        records = records + posGt->_gts->at(i)->alleles->at(altIdx1) +
                  "," +
                  posGt->_gts->at(i)->alleles->at(altIdx2) +
                  "\t.\t" +
                  "PASS\t.\tGT\t1/2\n";
      }
      // alt1 == alt && alt2 == ref: 1/0
      else if ( altIdx1 > 0 && altIdx2 == 0 ) {
        records = records + posGt->_gts->at(i)->alleles->at(altIdx1) + "\t.\t" +
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
      // altIdx is the reference allele: ignore
      records = records + posGt->_chr + "\t" +
                to_string(posGt->_gts->at(i)->_pos+1) + "\t.\t" +
                posGt->_gts->at(i)->alleles->at(0) + "\t" +
                posGt->_gts->at(i)->alleles->at(altIdx) + "\t.\t" +
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
    sum += absK;
  }

  return sum / kstr.size();
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

  return kstr.at(i + ((kstr.size() - i)/2));
}

