
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
#include <vector>


gtAllele::gtAllele(vcfRecord *record) {
    
    _pos  = record->_pos - 1;
    alleles = new vector<char*>;
    
    //  fprintf(stderr, "[ DEBUG ] :: (*record->_arr_samples)[0] = %s\n", (*record->_arr_samples)[0]);

    if ( strncmp((*record->_arr_samples)[0], "./.", 3) == 0 
      || strncmp((*record->_arr_samples)[0], "0/0", 3) == 0) {
      record->isValid = false;
      //  fprintf(stderr, "[ DEBUG ] :: %s is invalid.\n", (*record->_arr_samples)[0]);
      return;
    }

    splitFields GT((*record->_arr_samples)[0], '/');
    alleles->push_back(record->_ref);   //  make the alleles.at(0) always be the ref allele
    _refLen    = strlen(record->_ref);
    parseGT(GT[0], record->_ref, record->_arr_alts, alleles, record->isValid);
    parseGT(GT[1], record->_ref, record->_arr_alts, alleles, record->isValid);

    /*  test
    //  fprintf(stderr, "[ DEBUG ] :: gtAlleles at pos=%u are %lu :", _pos, alleles->size());
    for ( int i = 0; i < alleles->size(); i++) {
      fprintf(stderr, " %d = %s", i, alleles->at(i));
    }
    fprintf(stderr, "\n");
    */
}


void
gtAllele::parseGT(char* gt, char *ref, splitFields *alts, vector <char*> *alleles, bool &isValid) {

    // is "." or ref allele
    if ( strcmp(gt, ".") == 0 ) {
      isValid = false;
      //  Do nothing
    } else if (strcmp(gt, "0") == 0 ) {
      //  Do nothing
    } else {
      uint32  altIdx  = strtouint32(gt) - 1;    // gt starts from 1, index from 0
      char* hap = (*alts)[altIdx];

      vector<char*>::iterator it = find (alleles->begin(), alleles->end(), hap);
      if ( it == alleles->end() )  //  only add if hap wasn't already added
        alleles->push_back(hap);
    }
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
  delete [] _id;
  delete [] _ref;
  delete [] _alts;
  delete [] _filter;
  delete [] _info;
  delete [] _formats;
  delete [] _samples;
}


void
vcfRecord::load(char *inLine) {

  splitFields W(inLine, '\t');
  
  if ( W.numWords() < 10 ) {
    isValid = false;
    return;
  }

  _chr      = new char [strlen(W[0]) + 1];
  _pos      = W.touint32(1);
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

  uint64 excluded = 0;

  while (AS_UTL_readLine(L, Llen, Lmax, F)) {

    // Header line?
    if (L[0] == '#') {
      // keep header lines
      _headers.push_back(L);

      // Count unique CHR ids
      if(strncmp(L, "##contig=<ID", strlen("##contig=<ID")) == 0) {
        _numChr++;
        
      }

      continue;	// No need to handle the header lines
    }

    record = new vcfRecord(L);

    //  check the record is valid
    //  exclude: 0/0, ./.
    if ( record->isInvalid() ) {
      excluded++;
      continue;
    }

    _records.push_back(record);

    chr = record->_chr;
    // fprintf(stderr, "[ DEBUG ] :: L = %s\n", L);

    //  is a new chromosome?
    if ( prevChr.empty() || chr.compare(prevChr) != 0 ) {
      _mapChrPosGT->insert(pair<string, vector<posGT*> *> (chr, new vector<posGT*>()));
      // fprintf(stderr, "[ DEBUG ] :: Made a new posGT list for %s\n", chr.c_str());
    }
    prevChr = chr;
  }

  AS_UTL_closeFile(F, inName);

  delete [] L;

  fprintf(stderr, "   Collected " F_SIZE_T " header lines.\n", _headers.size());
  fprintf(stderr, "   Loaded " F_SIZE_T " records with %lu unique contig(s) from  %u contig IDs.\n", _records.size(), _mapChrPosGT->size(), _numChr);
  fprintf(stderr, "   while excluding %lu invalid records\n\n", excluded);

  // Iterate through the records and get per chr posGTs
  for (uint32 ii = 0; ii < _records.size(); ii++) {
    //  fprintf(stderr, "[ DEBUG ] :: _records.at(ii)->_pos : %u\n", _records.at(ii)->_pos);
    chr = string(_records[ii]->_chr);
    _mapChrPosGT->at(chr)->push_back(new posGT(_records[ii]));
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
vcfFile::mergeChrPosGT(uint32 ksize, uint32 comb, bool nosplit) {

  uint32 K_OFFSET  =  ksize--;	// ksize - 1

  map<string , vector<posGT*> *>::iterator it;

  #pragma omp parallel private(it)
  {
  for ( it = _mapChrPosGT->begin(); it != _mapChrPosGT->end(); it++ ) {
    //  for each chromosome - posGTlist
    #pragma omp single nowait
    {

    //  Initialize variables
    int removed   = 0;
    string chr    = it->first;
    vector<posGT*> *posGTlist  = it->second;
    int posGtSizeB = posGTlist->size();	// Before
    int posGtSizeA = posGTlist->size();	// After

    // fprintf(stderr, "[ DEBUG ] :: Merge variants in %s ... \n", chr.c_str());
    if ( posGtSizeB == 1 ) {
      fprintf(stderr, "%s : Nothing to merge. Only 1 variant found.\n", chr.c_str());
    }else{

    //  Get first start and end
    uint32 start     = posGTlist->at(0)->_rStart;
    uint32 end       = posGTlist->at(0)->_rEnd;

    //  for each posGT,
    //  keep order when iterating,
    //  begin from [1] not [0]
    for (int ii=1; ii<posGtSizeA;) {
      uint32 iStart = posGTlist->at(ii)->_rStart;
      uint32 iEnd   = posGTlist->at(ii)->_rEnd;

      //  ignore gts with 0 alleles (0/0 or ./.)
      vector<gtAllele*> *gts = posGTlist->at(ii)->_gts;
      if (gts->at(0)->alleles->size() == 0 ) {
        // fprintf(stderr, "[ DEBUG ] :: erase posGTlist->(%d) \n", ii);
        posGTlist->erase(posGTlist->begin() + ii);
        posGtSizeA--;
        removed++;
        continue;
      }

      //  fprintf(stderr, "[ DEBUG ] :: start=%u, end=%u, iStart=%u, iEnd=%u\n", start, end, iStart, iEnd);

      //  Is ii overlapping with ii-1 th posGT ?
      //  Move - K_OFFSET from left to right (or vice versa)
      //  to prevent overflow
      if (
          (
           ( iStart < end + (2 * K_OFFSET) && start < iStart)
           || // In case not sorted
           ( iEnd + (2 * K_OFFSET) > start && iEnd < end )
          )
          && 
          (
           (posGTlist->at(ii-1)->size() <= comb)
           ||
           nosplit
          )
         ) {

        //  Add gts to previous posGT
        //  fprintf(stderr, "[ DEBUG ] :: Num. alleles = %u. Allele at %u = %s\n", gts->at(0)->alleles->size(), 0, gts->at(0)->alleles->at(0));
        posGTlist->at(ii-1)->addGtAllele(gts->at(0));
        // fprintf(stderr, "[ DEBUG ] :: Adding allele at %u to %u\n", gts->at(0)->_pos, posGTlist->at(ii-1)->_gts->at(0)->_pos);
     
        //  fprintf(stderr, "[ DEBUG ] :: Total gts at pos %u: %u. Total gts = %u\n", posGTlist->at(ii-1)->_gts->at(0)->_pos, posGTlist->at(ii-1)->_gts->size(), posGTlist->size());
        //  Remove posGTlist[ii]
        // fprintf(stderr, "[ DEBUG ] :: Remove allele at %u\n", ii);
        posGTlist->erase(posGTlist->begin() + ii);
        //  fprintf(stderr, "[ DEBUG ] :: Now allele at %u is %u. Total gts = %u\n", ii, posGTlist->at(ii)->_gts->at(0)->_pos, posGTlist->size());
        posGtSizeA--;
        removed++;

				//  Extend end
				if ( iEnd > end )
					end = iEnd;

      } else if (posGTlist->at(ii-1)->size() > comb && !nosplit) {
      
        fprintf(stderr, "---%s : More than %u variants in the combination when variant at position %u is included. Splitting. Consider filtering the vcf upfront.\n", chr.c_str(), posGTlist->at(ii-1)->size(), posGTlist->at(ii)->_rStart);
        start = posGTlist->at(ii)->_rStart;
        end   = posGTlist->at(ii)->_rEnd;
        ii++;
      
      } else {
        start = posGTlist->at(ii)->_rStart;
        end   = posGTlist->at(ii)->_rEnd;
        ii++;
      }
    }
    fprintf(stderr, "%s : Reduced %d variants down to %d combinations for evaluation (merged %d)\n", chr.c_str(), posGtSizeB, posGtSizeA, removed);
  }
  }
  }
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


