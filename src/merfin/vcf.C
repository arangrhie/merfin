
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
  if (inLine[0] == '#')  return;
  
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

  fprintf(stderr, "  Loaded " F_SIZE_T " records with %lu unique contig(s) from  %u contig IDs.\n\n", _records.size(), _mapChrPosGT->size(), _numChr);

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


