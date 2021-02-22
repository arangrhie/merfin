
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
#include "files.H"
#include "vcf.H"
#include "arrays.H"
#include <vector>


void
gtAllele::parseGT(char const *gt) {

  //  Is "." or ref allele?  Do nothing.
  if (strcmp(gt, ".") == 0) {
    _record->invalidate();
    return;
  }

  //  Add if hap isn't already there (and if it's a valid index).
  //  gt is 1-based while alts[] is 0-based.

  int32   altIdx = strtoint32(gt);

  if (altIdx > 0) {
    char   *hap = _record->_arr_alts[altIdx - 1];

    if (find(_alleles->begin(), _alleles->end(), hap) == _alleles->end())
      _alleles->push_back(hap);
  }
}



gtAllele::gtAllele(vcfRecord *r) {

  _record = r;
  _pos    = _record->get_pos() - 1;

  //  fprintf(stderr, "[ DEBUG ] :: (*record->_arr_samples)[0] = %s\n", (*record->_arr_samples)[0]);

  if ((strncmp(_record->_arr_samples[0], "./.", 3) == 0) ||
      (strncmp(_record->_arr_samples[0], "0/0", 3) == 0)) {
    //  fprintf(stderr, "[ DEBUG ] :: %s is invalid.\n", (*record->_arr_samples)[0]);
    _record->invalidate();
    return;
  }

  splitToWords GT(_record->_arr_samples[0], splitLetter, '/');

  _alleles = new vector<char const *>;
  _alleles->push_back(_record->get_ref());   //  make the alleles.at(0) always be the ref allele

  _refLen    = strlen(_record->get_ref());

  parseGT(GT[0]);  //, record->get_ref(), record->_arr_alts, alleles, record->isValid);
  parseGT(GT[1]);  //, record->get_ref(), record->_arr_alts, alleles, record->isValid);

  //fprintf(stderr, "[ DEBUG ] :: gtAlleles at pos=%u are %lu :", _pos, alleles->size());
  //for ( int i = 0; i < alleles->size(); i++)
  //  fprintf(stderr, " %d = %s", i, alleles->at(i));
  //fprintf(stderr, "\n");
}





bool
vcfFile::loadFile(char *inName) {
  char  *L    = NULL;
  uint32 Llen = 0;
  uint32 Lmax = 0;

  compressedFileReader F(inName);

  uint64 excluded = 0;

  while (AS_UTL_readLine(L, Llen, Lmax, F.file())) {

    //  Save all the header lines as is - push_back is copying L into a
    //  string.  Count the number of unique CHR ids.  Then get another line.

    if (L[0] == '#') {
      _headers.push_back(L);

      if (strncmp(L, "##contig=<ID", 12) == 0)
        _numChr++;

      continue;
    }

    //  Convert this line into a vcfRecord.  If it's not valid, destroy the
    //  record and get another line.

    vcfRecord *record = new vcfRecord(L);

    if (record->isValid() == false) {
      excluded++;
      delete record;
      continue;
    }

    //  A good record!  Save it and make a new vector of posGT's if this is
    //  the first time we've seen the CHR.

    _records.push_back(record);

    string  chr = record->get_chr();

    if (_mapChrPosGT.count(chr) == 0)
      _mapChrPosGT[chr] = new vector<posGT *>;

    //  Add a new posGT for this record.

    _mapChrPosGT[chr]->push_back(new posGT(record));
  }

  delete [] L;

  fprintf(stderr, "   Collected " F_SIZE_T " header lines.\n", _headers.size());
  fprintf(stderr, "   Loaded " F_SIZE_T " records:\n", _records.size());
  fprintf(stderr, "      %-8lu unique contig%s\n", _mapChrPosGT.size(), (_mapChrPosGT.size() == 1) ? "" : "s");
  fprintf(stderr, "      %-8u contig IDs\n", _numChr);
  fprintf(stderr, "   Excluded %lu invalid records\n", excluded);
  fprintf(stderr, "\n");

  return(true);
}




bool
vcfFile::mergeChrPosGT(uint32 ksize, uint32 comb, bool nosplit) {

  uint32 K_OFFSET  =  ksize--;	// ksize - 1

    //  for each chromosome - posGTlist
  for (auto it = _mapChrPosGT.begin(); it != _mapChrPosGT.end(); it++ ) {
    string          chr       =  it->first;
    vector<posGT*> &posGTlist = *it->second;   //  A reference to the vector

    //  Nothing to do if there is only one thing on the list!
    if (posGTlist.size() == 1)
      continue;

    //  Initialize variables
    int removed   = 0;
    int posGtSizeB = posGTlist.size();	// Before
    int posGtSizeA = posGTlist.size();	// After

    // fprintf(stderr, "[ DEBUG ] :: Merge variants in %s ... \n", chr.c_str());

    //  Get first start and end
    uint32 start     = posGTlist.at(0)->_rStart;
    uint32 end       = posGTlist.at(0)->_rEnd;

    //  for each posGT,
    //  keep order when iterating,
    //  begin from [1] not [0]
    for (int ii=1; ii<posGtSizeA;) {
      uint32 iStart = posGTlist[ii]->_rStart;
      uint32 iEnd   = posGTlist[ii]->_rEnd;

      //  ignore gts with 0 alleles (0/0 or ./.)
      vector<gtAllele*> *gts = posGTlist[ii]->_gts;
      if (gts->at(0)->_alleles->size() == 0 ) {
        // fprintf(stderr, "[ DEBUG ] :: erase posGTlist.(%d) \n", ii);
        posGTlist.erase(posGTlist.begin() + ii);
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
           (posGTlist.at(ii-1)->size() < comb)
           ||
           nosplit
           )
          ) {

        //  Add gts to previous posGT
        //  fprintf(stderr, "[ DEBUG ] :: Num. alleles = %u. Allele at %u = %s\n", gts->at(0)->alleles->size(), 0, gts->at(0)->alleles->at(0));
        posGTlist.at(ii-1)->addGtAllele(gts->at(0));
        // fprintf(stderr, "[ DEBUG ] :: Adding allele at %u to %u\n", gts->at(0)->_pos, posGTlist.at(ii-1)->_gts->at(0)->_pos);
     
        //  fprintf(stderr, "[ DEBUG ] :: Total gts at pos %u: %u. Total gts = %u\n", posGTlist.at(ii-1)->_gts->at(0)->_pos, posGTlist.at(ii-1)->_gts->size(), posGTlist.size());
        //  Remove posGTlist[ii]
        // fprintf(stderr, "[ DEBUG ] :: Remove allele at %u\n", ii);
        posGTlist.erase(posGTlist.begin() + ii);
        //  fprintf(stderr, "[ DEBUG ] :: Now allele at %u is %u. Total gts = %u\n", ii, posGTlist[ii]->_gts->at(0)->_pos, posGTlist.size());
        posGtSizeA--;
        removed++;

        //  Extend end
        if ( iEnd > end )
          end = iEnd;

      } else if (posGTlist.at(ii-1)->size() >= comb && !nosplit) {
      
        fprintf(stderr, "---%s : More than %u variants in the combination when variant at position %u is included. Splitting. Consider filtering the vcf upfront.\n", chr.c_str(), posGTlist.at(ii-1)->size(), posGTlist[ii]->_rStart);
        start = posGTlist[ii]->_rStart;
        end   = posGTlist[ii]->_rEnd;
        ii++;
      
      } else {
        start = posGTlist[ii]->_rStart;
        end   = posGTlist[ii]->_rEnd;
        ii++;
      }
    }
    fprintf(stderr, "%s : Reduced %d variants down to %d combinations for evaluation (merged %d)\n", chr.c_str(), posGtSizeB, posGtSizeA, removed);
  }

  return(true);
}
