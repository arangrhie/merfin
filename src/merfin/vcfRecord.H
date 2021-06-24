
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

#ifndef VCFRECORD_H
#define VCFRECORD_H

#include "types.H"
#include "strings.H"


class vcfRecord {
public:
  vcfRecord()  {}
  ~vcfRecord() {}

  bool    isValid(void) {
    return(_isValid);
  }

  void    invalidate(void) {
    _isValid = false;
  }

  //  Parse a VCF input record line and save pointers to the values in the line.
  //  Returns true/false if the line is a valid record.
  //  The _words member holds a copy of the input line.
  //
  //  A record is:
  //    [---------------------------- mandatory ---------------------------]    [genotype data, any number of IND columns]
  //    #CHROM          POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  IND
  //    scaffold_607    16044   .       T       TA      3       q40     .       GT:DP   1/1:6
  //
  //  GT - genotype with alleles separated by / (unphased) or | (phased)
  //       1 for first allele listed in ALT,
  //       2 for second allele, etc
  //  DP - read depth
  //
  bool    load(char *inLine) {
    _words.split(inLine, '\t');

    if (_words.numWords() < 10) {
      _words.erase();
      return(false);
    }

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

    _arr_alts.   split(_alts,    ',');
    _arr_formats.split(_formats, ':');
    _arr_samples.split(_samples, ':');

    _isValid  = true;

    return(true);
  }

  void    save(compressedFileWriter* outFile) {
    fprintf(outFile->file(), "%s\t%d\t%s\t%s\t%s\t%.1f\t%s\t%s\t%s\t%s\n",
            _chr, _pos, _id, _ref, _alts, _qual, _filter, _info, _formats, _samples);
  }

  string  save(void) {
    uint32 len = (strlen(_chr) +
                  12 +
                  strlen(_id) +
                  strlen(_ref) +
                  strlen(_alts) +
                  12 +
                  strlen(_filter) +
                  strlen(_info) +
                  strlen(_formats) +
                  strlen(_samples));
    char  line[len + 128];

    sprintf(line, "%s\t%d\t%s\t%s\t%s\t%.1f\t%s\t%s\t%s\t%s\n",
            _chr, _pos, _id, _ref, _alts, _qual, _filter, _info, _formats, _samples);

    return(string(line));
  }

  //  To enforce that this is a read-only object, there are only accessors
  //  that don't allow modifications.
  //
  //  Not all data is even accessed.
  //  
public:
  char const   *get_chr(void)          { return(_chr);         }
  uint32        get_pos(void)          { return(_pos);         }
  //char const   *get_id(void)           { return(_id);          }
  char const   *get_ref(void)          { return(_ref);         }
  //char         *get_alts(void)         { return(_alts);        }
  double        get_qual(void)         { return(_qual);        }
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
  double        _qual        = 0.0;        //  QUAL
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

private:
  bool          _isValid = false;          //  has non-sense GT?
};

#endif  //  VCFRECORD_H
