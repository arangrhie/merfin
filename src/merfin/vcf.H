
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

#include "strings.H"

#include <string>
#include <algorithm>
#include <map>
using namespace std;




class vcfRecord {
public:
  vcfRecord(char *inLine) {
    load(inLine);
  }

  vcfRecord() {};
  ~vcfRecord() {}

  bool    isValid(void) {
    return(_isValid);
  }

  void    invalidate(void) {
    _isValid = false;
  }

  void    load(char *inLine) {
    _words.split(inLine, splitLetter, '\t');

    _isValid = (_words.numWords() >= 10);

    //  If not valid, release any memory in splitToWords.
    if (_isValid == false) {
      _words.erase();
    }

    //  If valid, populate our pointers.
    else {
      _chr      = _words[0];
      _pos      = _words.touint32(1);
      _id       = _words[2];
      _ref      = _words[3];
      _alts     = _words[4];
      _qual     = _words.todouble(5);
      _filter   = _words[6];
      _info     = _words[7];
      _formats  = _words[8];
      _samples  = _words[9];

      _arr_alts.   split(_alts,    splitLetter, ',');
      _arr_formats.split(_formats, splitLetter, ':');
      _arr_samples.split(_samples, splitLetter, ':');
    }
  }

  void    save(compressedFileWriter* outFile) {
    fprintf(outFile->file(), "%s\t%d\t%s\t%s\t%s\t%.1f\t%s\t%s\t%s\t%s\n",
            _chr, _pos, _id, _ref, _alts, _qual, _filter, _info, _formats, _samples);
  }

public:
  char const   *get_chr(void)          { return(_chr);         }
  uint32        get_pos(void)          { return(_pos);         }
  //char const   *get_id(void)           { return(_id);          }
  char const   *get_ref(void)          { return(_ref);         }
  //char         *get_alts(void)         { return(_alts);        }
  //float         get_qual(void)         { return(_qual);        }
  //char         *get_filter(void)       { return(_filter);      }
  //char         *get_info(void)         { return(_info);        }
  //char         *get_formats(void)      { return(_formats);     }
  //char         *get_samples(void)      { return(_samples);     }
  //int           get_size_alts(void)    { return(_size_alts);   }
  //int           get_size_format(void)  { return(_size_format); }

private:
  splitToWords  _words;                    //  Input line, the rest are pointers into this.

  char         *_chr         = nullptr;    //  CHR
  uint32        _pos         = uint32max;  //  POS
  char         *_id          = nullptr;    //  ID
  char         *_ref         = nullptr;    //  REF
  char         *_alts        = nullptr;    //  ALT
  float         _qual        = 0.0;        //  QUAL
  char         *_filter      = nullptr;    //  FILTER
  char         *_info        = nullptr;    //  INFO
  char         *_formats     = nullptr;    //  FORMAT
  char         *_samples     = nullptr;    //  SAMPLE
  int           _size_alts   = 0;          //  Num. alt variants
  int           _size_format = 0;          //  Num. format values

public:
  splitToWords  _arr_alts;
  splitToWords  _arr_formats;
  splitToWords  _arr_samples;

public:
  bool          _isValid = false;          //  has non-sense GT?
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
  vcfRecord            *_record;
  uint32                _pos;
  uint32                _refLen;
  vector<char const *> *_alleles;

private:
  void parseGT(char const *gt);
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
    _chr         = record->get_chr();
    _rStart      = record->get_pos() - 1;
    _rEnd        = record->get_pos();

    _gts         = new vector<gtAllele*>;
    gtAllele *gt = new gtAllele(record);

    addGtAllele(gt);
  };

  ~posGT();

public:
  char const        *_chr;
  uint32             _rStart;
  uint32             _rEnd;
  vector<gtAllele*> *_gts;

public:
  int    size() {
    return(_gts->size());
  };

public:
  void    addGtAllele(gtAllele *gt) {
    _gts->push_back(gt);

    _rStart = min(_rStart, gt->_pos);
    _rEnd   = max(_rEnd,   gt->_pos + gt->_refLen);
  };
};


class vcfFile {
public:
  vcfFile(char *inName) {
    loadFile(inName);
  }

  ~vcfFile() {
    for (uint32 ii=0; ii<_records.size(); ii++)
      delete _records[ii];
  }

  bool    loadFile(char *inName);

  bool    saveFile(compressedFileWriter *outFile) {
    for (uint32 ii=0; ii<_records.size(); ii++)
      if (_records[ii])
        _records[ii]->save(outFile);

    return(true);
  }

  //  Merge ChrPosGT when POS are within ksize
  bool            mergeChrPosGT(uint32 ksize, uint32 comb, bool nosplit);

  vector<string> getHeaders(void)  {
    return _headers;
  };

public:
  int32                           _numChr      = 0;  //  num. of CHR entries
  vector<string>                  _headers;
  vector<vcfRecord *>             _records;
  map<string, vector<posGT *>* >  _mapChrPosGT; 
};


#endif  //  AS_UTL_VCF_H
