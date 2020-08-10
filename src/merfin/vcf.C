
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

#include "runtime.H"
#include "files.H"
#include "vcf.H"
#include "arrays.H"
#include <algorithm>
#include <vector>



gtAllele::gtAllele(vcfRecord *record) {
    
    _pos  = record->_pos - 1;
    
    splitFields GT((*record->_arr_samples)[0], '/');
    alleles = new vector<char*>;
    alleles->push_back(record->_ref);   //  make the alleles.at(0) always be the ref allele
    _refLen    = strlen(record->_ref);
    parseGT(GT[0], record->_ref, record->_arr_alts, alleles);
    parseGT(GT[1], record->_ref, record->_arr_alts, alleles);
    
    // _maxRefLen = (strlen(_hap1) > strlen(_hap2)) ? strlen(_hap1) : strlen(_hap2);
}


void
gtAllele::parseGT(char* gt, char *ref, splitFields *alts, vector <char*> *alleles) {

    // is "." or ref allele
    if ( strcmp(gt, ".") == 0 || strcmp(gt, "0") == 0 ) {
      //  Do nothing
    } else {
      uint32  altIdx  = strtouint32(gt) - 1;    // gt starts from 1, index from 0
      char* hap = (*alts)[altIdx];

      vector<char*>::iterator it = find (alleles->begin(), alleles->end(), hap);
      if ( it == alleles->end() )  //  only add if hap wasn't already added
        alleles->push_back(hap);
    }
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
varMer::score(merylExactLookup *rlookup, merylExactLookup *alookup, uint64 peak) {

  //  iterate through each base and get kmer
  uint64 fValue = 0;
  uint64 rValue = 0;
  uint64 freq;
  uint32 numM;  // num. missing kmers
  uint64 minF;
  double minK;
  double maxK;
  string seq;
  bool   fExists = false;
  bool   rExists = false;
  double readK;
  double asmK;
  double kMetric;
  vector<double> m_ks;
  bool  evaluate = true;

  uint32 idx = 0;	// var index in the seqe

  //  get scores at each kmer pos and minimum read multiplicity
  for ( int ii = 0; ii < seqs.size(); ii++ ) {
    minF  = UINT64_MAX;
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
      freq = 0;
      fValue = 0;
      rValue = 0;
      if (kiter.isValid()) {
        fExists = rlookup->exists(kiter.fmer(), fValue);
        rExists = rlookup->exists(kiter.rmer(), rValue);
      }

      if ( fExists || rExists ) {
          freq = fValue + rValue;
      }

      if (minF > freq) { minF = freq; };

      readK = (double) freq / peak;

      fValue = 0;
      rValue = 0;
      if (kiter.isValid()) {
        alookup->exists(kiter.fmer(), fValue);
        alookup->exists(kiter.rmer(), rValue);
      }

      asmK  = (double) (fValue + rValue);
      
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

      if (freq == 0) {
        kMetric = -1;  // use 0 if we are using non-abs k*
        numM++;

        // no need to evaluate or give scores
        evaluate = false;
        // continue;
        
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
      minFs.insert(pair<uint64, int>(minF, ii));      // Automatically sorted by min value
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

vcfRecord::vcfRecord() {
  _chr         = NULL;
  _pos         = UINT32_MAX;
  _id          = NULL;
  _ref         = NULL;
  _alts        = NULL;
  _qual        = 0;
  _filter      = NULL;
  _info        = NULL;
  _formats     = NULL;
  _samples     = NULL;
  _size_alts   = 0;
  _size_format = 0;
}


vcfRecord::vcfRecord(char *inLine) {
  load(inLine);
}


vcfRecord::~vcfRecord() {
  delete [] _chr;
  delete [] _ref;
  delete [] _alts;
  delete [] _filter;
  delete [] _info;
  delete [] _formats;
  delete [] _samples;
}


void
vcfRecord::load(char *inLine) {

  // Skip header lines
  if (inLine[0] == '#')  {
    
    return;
  }
  
  splitFields W(inLine, '\t');

  _chr      = new char [strlen(W[0]) + 1];
  _pos      = W.toint32(1);
  _id       = new char [strlen(W[2]) + 1];
  _ref      = new char [strlen(W[3]) + 1];
  _alts     = new char [strlen(W[4]) + 1];
  _qual     = W.todouble(5);
  _filter   = new char [strlen(W[6]) + 1];
  _info     = new char [strlen(W[7]) + 1];
  _formats  = new char [strlen(W[8]) + 1];
  _samples  = new char [strlen(W[9]) + 1];

  strcpy(_chr,    W[0]);
  strcpy(_id,     W[2]);
  strcpy(_ref,    W[3]);
  strcpy(_alts,   W[4]);
  strcpy(_filter, W[6]);
  strcpy(_info,   W[7]);
  strcpy(_formats,W[8]);
  strcpy(_samples,W[9]);

  _arr_alts    = new splitFields(_alts, ',');
  _arr_formats = new splitFields(_formats, ':');
  _arr_samples = new splitFields(_samples, ':');

}


void
vcfRecord::save(FILE *outFile) {
  fprintf(outFile, "%s\t%d\t%s\t%s\t%s\t%f\t%s\t%s\t%s\t%s\n",
          _chr, _pos, _id, _ref, _alts, _qual, _filter, _info, _formats, _samples);
}

vcfFile::vcfFile(char *inName) {
  _numChr     = 0;
  _fName      = inName;
  _mapChrPosGT = new map<string, vector<posGT*>*>();
  loadFile(inName);
}


vcfFile::~vcfFile() {
  for (uint32 ii=0; ii<_records.size(); ii++)
    delete _records[ii];
}


bool
vcfFile::loadFile(char *inName) {
  char  *L    = NULL;
  uint32 Llen = 0;
  uint32 Lmax = 0;

  FILE *F = AS_UTL_openInputFile(inName);

  vcfRecord *record = NULL;
  string chr;
  string prevChr;

  while (AS_UTL_readLine(L, Llen, Lmax, F)) {

    // Header line?
    if (L[0] == '#') {
      // keep header lines
      _headers.push_back(string(L));

      // Count unique CHR ids
      if(strncmp(L, "##contig=<ID", strlen("##contig=<ID")) == 0) {
        _numChr++;
        
      }
      continue;	// No need to handle the header lines
    }
    record = new vcfRecord(L);
    _records.push_back(record);
    chr = string(record->_chr);

    if ( prevChr.empty() || chr.compare(prevChr) != 0 ) {
      _mapChrPosGT->insert(pair<string, vector<posGT*> *> (chr, new vector<posGT*>()));
      // fprintf(stderr, "[ DEBUG ] :: Made a new posGT list for %s\n", chr.c_str());
    }
    prevChr = chr;
  }

  AS_UTL_closeFile(F, inName);

  delete [] L;

  fprintf(stderr, "   Collected " F_SIZE_T " header lines.\n", _headers.size());
  fprintf(stderr, "   Loaded " F_SIZE_T " records with %lu unique contig(s) from  %u contig IDs.\n\n", _records.size(), _mapChrPosGT->size(), _numChr);

  // Iterate through the records and get per chr posGTs
  for (uint32 ii = 0; ii < _records.size(); ii++) {
    if (_records[ii]) {
      //  fprintf(stderr, " DEBUG :: _records[ii]->_chr : %s", _records[ii]->_chr);
      chr = string(_records[ii]->_chr);
      _mapChrPosGT->at(chr)->push_back(new posGT(_records[ii]));
    }
  }

  return(true);
}


bool
vcfFile::saveFile(char *outName) {

  FILE *F = AS_UTL_openOutputFile(outName);

  for (uint32 ii=0; ii<_records.size(); ii++)
    if (_records[ii])
      _records[ii]->save(F);

  AS_UTL_closeFile(F, outName);

  return(true);
}




bool
vcfFile::mergeChrPosGT(uint32 ksize) {

  uint32 K_OFFSET  =  ksize--;	// ksize - 1

  map<string , vector<posGT*> *>::iterator it;

  string  chr;

  vector<posGT*>  *posGTlist;
  vector<gtAllele*>     *gts;

  uint32  start  = 0;
  uint32  end    = 0;
  uint32  iStart = 0;	// start at iterating pos
  uint32  iEnd   = 0;   // end   at iterating pos

  int    posGtSizeB = 0;
  int    posGtSizeA = 0;
  int    removed   = 0;

  for ( it = _mapChrPosGT->begin(); it != _mapChrPosGT->end(); it++ ) {
    //  for each chromosome - posGTlist

    //  Initialize variables
    removed    = 0;
    chr        = it->first;
    posGTlist  = it->second;
    posGtSizeB = posGTlist->size();	// Before
    posGtSizeA = posGTlist->size();	// After

    if ( posGtSizeB == 1 ) {
      fprintf(stderr, "%s : Nothing to merge. Only 1 variant found.\n", chr.c_str());
      continue;
    }

    //  Get first start and end
    start     = posGTlist->at(0)->_rStart;
    end       = posGTlist->at(0)->_rEnd;

    //  for each posGT,
    //  keep order when iterating,
    //  begin from [1] not [0]
    for (int ii=1; ii<posGtSizeA;) {
      iStart = posGTlist->at(ii)->_rStart;
      iEnd   = posGTlist->at(ii)->_rEnd;


      //  fprintf(stderr, "[ DEBUG ] :: start=%u, end=%u, iStart=%u, iEnd=%u\n", start, end, iStart, iEnd);

      //  Is ii overlapping with ii-1 th posGT ?
      //  Move - K_OFFSET from left to right (or vice versa)
      //  to prevent overflow
      if (( iStart < end + (2 * K_OFFSET) && start < iStart)
         || // In case not sorted
          ( iEnd + (2 * K_OFFSET) > start && iEnd < end )) {

        //  Add gts to previous posGT
        gts = posGTlist->at(ii)->_gts;
        for (uint32 jj=0; jj < gts->size(); jj++) {
            posGTlist->at(ii-1)->addGtAllele(gts->at(jj));
        //  fprintf(stderr, "[ DEBUG ] :: Adding allele at %u to %u\n", gts->at(jj)->_pos, posGTlist->at(ii-1)->_gts->at(0)->_pos);
        }
     
        //  fprintf(stderr, "[ DEBUG ] :: Total gts at pos %u: %u\n", posGTlist->at(ii-1)->_gts->at(0)->_pos, posGTlist->at(ii-1)->_gts->size());
        //  Remove posGTlist[ii]
        //  fprintf(stderr, "[ DEBUG ] :: Remove allele at %u\n", ii);
        posGTlist->erase(posGTlist->begin() + (ii));
        //  fprintf(stderr, "[ DEBUG ] :: Now allele at %u is %u\n", ii, posGTlist->at(ii)->_gts->at(0)->_pos);
        posGtSizeA--;
        removed++;

      } else {
        start = posGTlist->at(ii)->_rStart;
        end   = posGTlist->at(ii)->_rEnd;
        ii++;
      }
    }
    fprintf(stderr, "%s : Merged %d variants from %d to %d\n", chr.c_str(), removed, posGtSizeB, posGtSizeA);
  }
  return(true);
}



splitFields::splitFields(const char *string, char delim) {
  _wordsLen  = 0;
  _wordsMax  = 0;
  _words     = NULL;

  _charsLen = 0;
  _charsMax = 0;
  _chars    = NULL;

  if (string)
    split(string, delim);

}


splitFields::~splitFields() {
  delete [] _chars;
  delete [] _words;
}


void
splitFields::split(const char *line, char delim) {

  _wordsLen = 0;        //  Initialize to no words
  _charsLen = 0;        //  and no characters.

  if (line == NULL)     //  Bail if there isn't a line to process.
    return;

  while (line[_charsLen] != 0)
    if (line[_charsLen++] == delim)
      _wordsLen++;

  resizeArray(_words, 0, _wordsMax, _wordsLen + 1, resizeArray_doNothing);
  resizeArray(_chars, 0, _charsMax, _charsLen + 1, resizeArray_doNothing);

  memset(_words, 0,    sizeof(char *) * (_wordsLen + 1));
  memcpy(_chars, line, sizeof(char)   * (_charsLen + 1));

  _wordsLen = 0;

  for (uint32 st=1, ii=0; ii < _charsLen; ii++) {
    if (line[ii] == delim) {                //  If the character is a word
      _chars[ii] = 0;                       //  separator, convert to NUL,
      st         = true;                    //  and flag the next character
    }                                       //  as the start of a new word.

    else if (st) {                          //  Otherwise, if this is the
      _words[_wordsLen++] = _chars + ii;    //  start of a word, make
      st                  = false;          //  a new word.
    }
  }
}


