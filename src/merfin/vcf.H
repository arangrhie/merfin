
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
  bool    isInvalid() { return !isValid;  };
  void    save(FILE *outFile);

public:
  char*         _chr;          //  CHR
  uint32        _pos;          //  POS
  char*         _id;           //  ID
  char*         _ref;          //  REF
  char*         _alts;         //  ALT
  float         _qual;         //  QUAL
  char*         _filter;       //  FILTER
  char*         _info;         //  INFO
  char*         _formats;      //  FORMAT
  char*         _samples;      //  SAMPLE
  int           _size_alts;    //  Num. alt variants
  int           _size_format;  //  Num. format values
  splitFields*  _arr_alts;
  splitFields*  _arr_formats;
  splitFields*  _arr_samples;
  bool           isValid = true;  //  has non-sense GT?
};


/**********************************************
 * Parse gt and infer haplotype 1 and 2 to alleles.
 * alleles[0] is the first allele (REF),
 * alleles[1] is the first non-reference allele (ALT).
 *********************************************/
class gtAllele {
public:
  gtAllele(vcfRecord *record);
  ~gtAllele();

public:
  uint32        _pos;
  uint32        _refLen;
  vector<char*> *alleles;

private:
  void parseGT(char* gt, char *ref, splitFields *alts, vector <char*> *alleles, bool &isValid);
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
  int    size() {
  
  return _gts->size();
  
  };

public:
  void    addGtAllele(gtAllele *gt) {
    _gts->push_back(gt);

    uint32 rStart = gt->_pos;                // already 0-based
    uint32 rEnd   = rStart + gt->_refLen;    // use refLen to get actual ref coords

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
  bool    mergeChrPosGT(uint32 ksize, uint32 comb, bool nosplit);	//  Merge ChrPosGT when POS are within ksize
  vector<string> getHeaders()  { return _headers; };

public:
  char                           *_fName;
  int32                           _numChr;  //  num. of CHR entries
  vector<string>                  _headers;
  vector<vcfRecord *>             _records;
  map<string, vector<posGT *>* > *_mapChrPosGT; 

};


#endif  //  AS_UTL_VCF_H
