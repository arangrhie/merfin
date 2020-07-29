
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

#ifndef AS_UTL_VCF_H
#define AS_UTL_VCF_H

#include "runtime.H"
#include "types.H"
#include "kmers.H"

#include <string>
#include <algorithm>
#include <map>
using namespace std;



class splitFields {
public:
  splitFields(const char *string=NULL, char delim='\t');
  ~splitFields();

public:
  void    split(const char *line, char delim='\t');  //  split by delim

  uint32  numWords(void)        { return(_wordsLen); };

  char   *operator[](uint32 i)  { return(first(i)); };

  char   *first(uint32 i=0)     { return((_wordsLen <= i) ? NULL : _words[i]);  };
  char   *shift(void)           {
    if (_wordsLen == 0)                     //  If no words, nothing to return.
      return(NULL);

    for (uint32 ii=1; ii<_wordsLen; ii++)   //  Shift all words down one place, moving
      swap(_words[ii-1], _words[ii]);       //  the word to shift off to the end.

    return(_words[--_wordsLen]);            //  Return the word we shifted out.
  };

  char   *last(uint32 i=0)      { return((_wordsLen == 0) ? NULL : _words[_wordsLen - i - 1]); };
  char   *pop(void)             { return((_wordsLen == 0) ? NULL : _words[--_wordsLen]);       };

  int32   toint32(uint32 i)     { return(strtoint32 (_words[i])); };
  uint32  touint32(uint32 i)    { return(strtouint32(_words[i])); };
  int64   toint64(uint32 i)     { return(strtoint64 (_words[i])); };
  uint64  touint64(uint32 i)    { return(strtouint64(_words[i])); };
  double  todouble(uint32 i)    { return(strtodouble(_words[i])); };

private:
  uint32    _wordsLen;
  uint32    _wordsMax;
  char    **_words;

  uint32    _charsLen;
  uint32    _charsMax;
  char     *_chars;
};


class vcfRecord {
public:
  vcfRecord();
  vcfRecord(char *inLine);
  ~vcfRecord();

  void    load(char *inLine);
  void    save(FILE *outFile);

public:
  char*         _chr;          //  CHR
  int32         _pos;           //  POS
  char*         _id;           //  ID
  char*         _ref;          //  REF
  char*         _alts;         //  ALT
  float         _qual;          //  QUAL
  char*         _filter;       //  FILTER
  char*         _info;         //  INFO
  char*         _formats;      //  FORMAT
  char*         _samples;      //  SAMPLE
  int           _size_alts;     //  Num. alt variants
  int           _size_format;   //  Num. format values
  splitFields*  _arr_alts;
  splitFields*  _arr_formats;
  splitFields*  _arr_samples;
};


/**********************************************
 * Parse gt and infer haplotype 1 and 2.
 * haplotype 1 is the first allele (REF),
 * haplotype 2 the second allele (ALT).
 *********************************************/
class gtAllele {
public:
  gtAllele(vcfRecord *record) {

    _pos  = record->_pos - 1;

    //  BPW - all three of these are pointers to objects, but the functions
    //  are wanting just objects, so we need to dereference the pointer to
    //  get to the object.
    //
    //  (*record->_arr_samples)[1] is saying dereference _arr_samples to get
    //  the object, then apply operator[] to return the string.
    //
    //  no simple fixes (except what's here) - you're in this ugly part of
    //  c++ where pointers-to-objects are 'needed' to create new objects, but
    //  then reference-to-object (aka just the object itself) is needed for
    //  access (with operator[]).
    //
    //  unfortunately, I'm used to the pointer-to-object style and just avoid
    //  using objects directly.  but you're java background is wanting to use
    //  objects directly.

    splitFields GT((*record->_arr_samples)[0], '/');
    _hasRef = false;	// will be set true if found in parseGT()
    _ref    = record->_ref;
    _hap1   = parseGT(GT[0], record->_ref, record->_arr_alts);
    _hap2   = parseGT(GT[1], record->_ref, record->_arr_alts);

    //  fprintf(stderr, "hap1=%s\n", _hap1);
    //  fprintf(stderr, "hap2=%s\n", _hap2);
    _refLen    = strlen(record->_ref);
    _maxRefLen = (strlen(_hap1) > strlen(_hap2)) ? strlen(_hap1) : strlen(_hap2);
    // fprintf(stderr, "_maxRefLen=%u\n", _maxRefLen);
  };
  ~gtAllele() {
    delete [] _hap1;
    delete [] _hap2;
  };

public:
  uint32      _pos;
  uint32      _refLen;
  uint32      _maxRefLen;
  char*       _hap1;	// inferred haplotype 1
  char*       _hap2;    // inferred haplotype 2
  char*       _ref;
  bool        _hasRef;	// has reference allele

private:
  char  *parseGT(char* gt, char *ref, splitFields *alts) {
    
    // is it "."?
    if ( strcmp(gt, ".") == 0 ) {
      return(".");
    } else if ( strcmp(gt, "0") == 0 ) {
      _hasRef = true;
      return(ref);
    } else {
      uint32  altIdx  = strtouint32(gt) - 1;	// gt starts from 1, index from 0
      return((*alts)[altIdx]);
    }
  };
};



/****************************************************
 *  _rStart : 0-base, ref coord. POS-1from vcf
 *  _rEnd   : 1-base, ref coord.
 *  _rStart ~ _rEnd : ref bases affected by this var
 *  
 ****************************************************/
class posGT {
public:
  posGT(vcfRecord *record) {
    _chr         = record->_chr;
    _rStart      = record->_pos - 1;
    _rEnd        = _rStart + 1;
    _gts         = new vector<gtAllele*>();
    gtAllele *gt = new gtAllele(record);

    addGtAllele(gt);

  };
  ~posGT();

public:
  uint32             _rStart;
  uint32             _rEnd;
  vector<gtAllele*> *_gts;

public:
  void    addGtAllele(gtAllele *gt) {
    _gts->push_back(gt);

    uint32 rStart = gt->_pos; // already 0-based
    uint32 rEnd   = rStart + gt->_refLen;	// use refLen to get actual ref coords

    _rStart = ( rStart < _rStart ) ? rStart : _rStart;
    _rEnd   = ( rEnd   > _rEnd   ) ? rEnd   : _rEnd;

    return;
  };

public:
  char* _chr;
};


class vcfFile {
public:
  vcfFile(char *inName);
  ~vcfFile();

  bool    loadFile(char *inName);
  bool    saveFile(char *outName);
  bool    mergeChrPosGT(uint32 ksize);	//  Merge ChrPosGT when POS are within ksize

public:
  char                         *_fName;
  int32                         _numChr;  //  num. of CHR entries
  vector<vcfRecord *>           _records;
  map<string, vector<posGT *>* > *_mapChrPosGT; 

};


class varMer {

public:

  varMer(posGT* posGt){
    this->posGt = posGt;
  };
  ~varMer();

  void addSeqPath(string seq, vector<int> idxPath) {
    vector<string>::iterator it;

    it = find (seqs.begin(), seqs.end(), seq);
    if ( it != seqs.end() ) { return; }

    // only insert elements if seq is a new sequence
    seqs.push_back(seq);
    paths.push_back(idxPath);
    return;
  };

  void score(merylExactLookup *rlookup, merylExactLookup *alookup, uint64 peak) {

    //  iterate through each base nad get kmer
    uint64 fValue = 0;
    uint64 rValue = 0;
    uint64 freq;
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

    //  get scores at each kmer pos and minimum read multiplicity
    for ( int jj = 0; jj < seqs.size(); jj++ ) {
      minF  = UINT64_MAX;
      minK  = DBL_MAX;
      maxK  = -2;
      evaluate = true;

      seq    = seqs.at(jj);
      m_ks.clear();
      fprintf(stderr, "%s\t%u\t%u\t%s\t", posGt->_chr, posGt->_rStart, posGt->_rEnd, seq.c_str());

      kmerIterator kiter((char*) seq.c_str(), seq.size());
      while (kiter.nextMer()) {
        freq = 0;
        fValue = 0;
        rValue = 0;

        fExists = rlookup->exists(kiter.fmer(), fValue);
        rExists = rlookup->exists(kiter.rmer(), rValue);
        if ( fExists || rExists ) {
            freq = fValue + rValue;
        }

        if (minF > freq) { minF = freq; };
        //  fprintf(stderr, "%lu ", freq);

        readK = (double) freq / peak;

        fValue = 0;
        rValue = 0;
        alookup->exists(kiter.fmer(), fValue);
        alookup->exists(kiter.rmer(), rValue);

        asmK  = (double) (fValue + rValue) + 1;  // + 1 as we are newly introducing a kmer

        if (freq == 0) {
          kMetric = -1;            // set to negative so we can filter out < 0 values
          // no need to evaluate or give scores
          evaluate = false;
          continue;
        } else if (readK > asmK) {
          kMetric = readK / asmK;  // we don't need to shift this to 0 here
        } else {
          kMetric = asmK  / readK;
        }

        m_ks.push_back(kMetric);
        fprintf(stderr, "%.2f ", kMetric);

        if (minK > kMetric) { minK = kMetric; };
        if (maxK < kMetric) { maxK = kMetric; };

      }

      //  only evaluate when no missing kmers are generated
      if (evaluate) {
        kstrs.push_back(m_ks);
        minFs.insert(pair<uint64, int>(minF, jj));	// Automatically sorted by min value
        minKs.insert(pair<double, int>(minK, jj));	// Automatically sorted by min value
        maxKs.insert(pair<double, int>(maxK, jj));	// Automatically sorted by min value
      }

      fprintf(stderr, "\n");
    }
    return;
  };

  multimap<uint64, int> getMinFreqs() { return minFs; };
  multimap<double, int> getMinKs()    { return minKs; };
  multimap<double, int> getMaxKs()    { return maxKs; };

public:
  vector< vector<int> >     paths;  //  contains idx
  vector< string >          seqs;
  vector< vector<double> >  kstrs;  //  k* values, kstrs.at(0) is the k* value of the first seq
  multimap<uint64, int >    minFs;  //  min freq, sorted from lowest to greatest, with idx on the second
  multimap<double, int >    minKs;  //  min k*,   sorted from lowest to greatest, with idx on the second
  multimap<double, int >    maxKs;  //  max k*,   sorted from lowest to greatest, with idx on the second
  posGT*               posGt;
  vector< uint32 >     order;
};


#endif  //  AS_UTL_VCF_H
