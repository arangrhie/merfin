
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


gtAllele::gtAllele(vcfRecord *r) {

  //  Initialize the easy stuff.

  _record = r;
  _pos    = _record->get_pos() - 1;
  _refLen = strlen(_record->get_ref());
  _qual   = _record->get_qual();

  //  

  if ((strncmp(_record->_arr_samples[0], "./.", 3) == 0) ||
      (strncmp(_record->_arr_samples[0], "0/0", 3) == 0)) {
    //  fprintf(stderr, "[ DEBUG ] :: %s is invalid.\n", (*record->_arr_samples)[0]);
    _record->invalidate();
    return;
  }

  splitToWords GT(_record->_arr_samples[0], splitLetter, '/');

  _alleles.push_back(_record->get_ref());   //  _alleles[0] is ALWAYS the reference allele.

  //  Add alternate alleles to the list, as long as they aren't already there.

  for (uint32 ii=0; ii<GT.numWords(); ii++) {
    int32   altIdx = strtoint32(GT[ii]);

    if (altIdx <= 0) {        //  This handles the case of gt being the empty string,
      _record->invalidate();  //  gt being "0" (or "00", etc), or gt being non-numeric,
      continue;               //  (and even invalid stuff like negative numbers).
    }

    //  Add the alternate allele for this variant to the list of alleles if it
    //  isn't already there.  Since the alleles come from the same source, we
    //  can just compare pointers.
    //
    //  Well, no, we cannot just compare pointers, since the reference allele will
    //  have a guaranteed different address to any alternate allele.

    char const *hap = _record->_arr_alts[altIdx - 1];

    //  Search for pointer-to-pointer matches between any alternate alleles.
    if (hap != nullptr) {
      for (uint32 ii=0; ii<_alleles.size(); ii++)
        if (_alleles[ii] == hap)
          hap = nullptr;
    }

    //  If not found, search for string matches between the reference and alternate.

    if (hap != nullptr) {
      if (strcmp(_alleles[0], hap) == 0)
        hap = nullptr;
    }

    if (hap != nullptr)           //  If it's null, we found it
      _alleles.push_back(hap);   //  on the list already.
  }
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

    //  Attempt to convert this line into a vcfRecord.  If it fails, delete the
    //  incomplete record; otherwise save it onto the master list of records
    //  and add a posGT to the per-chromosome list.

    vcfRecord *record = new vcfRecord;

    if (record->load(L) == false) {
      excluded++;
      delete record;
    }
    else {
      _records.push_back(record);

      string  chr = record->get_chr();

      if (_mapChrPosGT.count(chr) == 0)
        _mapChrPosGT[chr] = new vector<posGT *>;

      _mapChrPosGT[chr]->push_back(new posGT(record));
    }
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




//  The goal is to merge posGT elements in the list that are within 2(k-1) of each other.
//
bool
vcfFile::mergeChrPosGT(uint32 ksize, uint32 comb, bool nosplit) {

  uint32 K_OFFSET =  2 * ksize;

    //  for each chromosome - posGTlist
  for (auto it = _mapChrPosGT.begin(); it != _mapChrPosGT.end(); it++ ) {
    string          chr       =  it->first;
    vector<posGT*> &inlist    = *it->second;         //  A reference to the vector
    vector<posGT*> *otlist    = new vector<posGT*>;  //  A new list with the merged posGT's

    uint32          removed   = 0;
    uint32          split     = 0;
    uint32          merged    = 0;

    //  Nothing to do if there is only one thing on the list!
    //if (inlist.size() == 1)
    //  continue;

    //  Sort the original list by begin position.
    auto byBeginCoord = [](posGT * const &A, posGT * const &B) { return(A->_rStart < B->_rStart); };

    sort(inlist.begin(), inlist.end(), byBeginCoord);

    //  Push the first posGT onto the output list.

    otlist->push_back(inlist[0]);

    //  Iterate over all the other input posGT's.

    for (uint32 ii=1; ii < inlist.size(); ii++) {

      //  Silently skip input any posGT that has no alleles.
      if (inlist[ii]->_gts.size() == 0) {
        removed++;
        continue;
      }

      //  We must be sorted.
      assert(otlist->back()->_rStart <= inlist[ii]->_rStart);

      //  If this input record does not overlap with the most recent output record,
      //  or the most recent output record has more than 'comb' items and splitting is allowed,
      //  make a new output.

      bool  overlapping = (inlist[ii]->_rStart < otlist->back()->_rEnd + K_OFFSET);
      bool  toomany     = (otlist->back()->_gts.size() >= comb);

      if (overlapping == false) {     //  Boring, just no overlap between variants.
        //fprintf(stderr, "%s : No overlap to previous variant at position %u-%u; make new cluster starting at position %u-%u.\n",
        //        chr.c_str(),
        //        otlist->back()->_rStart, otlist->back()->_rEnd,
        //        inlist[ii]->_rStart, inlist[ii]->_rEnd);
        otlist->push_back(inlist[ii]);
        continue;
      }

      if ((toomany == true) &&        //  Exciting!  A big pile of variants at this position,
          (nosplit == false)) {       //  and we're allowed to split large clusters.
        //fprintf(stderr, "%s : More than %u variants at position %u-%u; split variants starting at position %u-%u into a new set.\n",
        //        chr.c_str(),
        //        comb,
        //        otlist->back()->_rStart, otlist->back()->_rEnd,
        //        inlist[ii]->_rStart, inlist[ii]->_rEnd);
        otlist->push_back(inlist[ii]);
        split++;
        continue;
      }

      //  Otherwise, we're overlapping AND allowed to merge, so do that.

      //fprintf(stderr, "%s : Merge variant at %u into cluster at %u-%u.\n",
      //        chr.c_str(),
      //        inlist[ii]->_rStart,
      //        otlist->back()->_rStart, otlist->back()->_rEnd);
      otlist->back()->addGtAllele(inlist[ii]->_gts[0]);
      merged++;
    }

    fprintf(stderr, "%s : Reduced %lu variants down to %lu combinations for evaluation:\n", chr.c_str(), inlist.size(), otlist->size());

    if (removed > 0)   fprintf(stderr, "%s :   Removed %u empty alleles.\n", chr.c_str(), removed);
    if (split   > 0)   fprintf(stderr, "%s :   Split   %u complicated combinations.\n", chr.c_str(), split);
    if (merged  > 0)   fprintf(stderr, "%s :   Merged  %u variants into combinations.\n", chr.c_str(), merged);

    delete _mapChrPosGT[chr];
    _mapChrPosGT[chr] = otlist;
  }

  return(true);
}
