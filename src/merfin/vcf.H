
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

#include "vcfRecord.H"


//  A single variant position, with alleles listed.
//    _allele[0] is the reference allele.
//    _allele[1] is the first alternate allele.
//
class gtAllele {
public:
  gtAllele(vcfRecord *record);
  ~gtAllele();

public:
  vcfRecord            *_record;    //  Pointer to the record this allele came from
  uint32                _pos;       //  Position of this allele in the reference
  uint32                _refLen;    //  Length of the reference allele
  double                _qual;      //  QUAL of the input record
  vector<char const *>  _alleles;   //  List of alleles; pointers into _record.
};


//  A collection of variant positions that are 'close enough' to
//  each other to be considered as one.
//
//  The start/end coordinates are just min/max of all the gtAllele
//  in this set.
//
class posGT {
public:
  posGT(vcfRecord *record) {
    gtAllele *gt = new gtAllele(record);

    _chr         = record->get_chr();

    _gts.push_back(gt);

    _rStart      = gt->_pos;
    _rEnd        = gt->_pos + gt->_refLen;
  };
  ~posGT();

public:
  void    addGtAllele(gtAllele *gt) {
    _gts.push_back(gt);

    _rStart = min(_rStart, gt->_pos);
    _rEnd   = max(_rEnd,   gt->_pos + gt->_refLen);
  };

public:
  char const        *_chr;      //  Pointer to the CHROM string in the input record.
  uint32             _rStart;   //  Reference start of this variant
  uint32             _rEnd;     //  Reference end of this variant
  vector<gtAllele*>  _gts;      //  List of alleles for this variant
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

  vector<string>  getHeaders(void)  {
    return _headers;
  };

public:
  int32                           _numChr      = 0;  //  num. of CHR entries
  vector<string>                  _headers;
  vector<vcfRecord *>             _records;

  //  Maps a chromosome name (string) to a vector.
  //  When loading the file, 
  map<string, vector<posGT *>* >  _mapChrPosGT; 
};


#endif  //  AS_UTL_VCF_H
