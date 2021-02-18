
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

#include "merfin-globals.H"

#include "varMer.H"




uint64 getIndex(vector<string> v, string K)
{
  auto it = find(v.begin(), v.end(), K);
  int index = distance(v.begin(), it);
  return index;
}

char*
concat(const char *s1, const char *s2) {
  char *result = (char*) malloc(strlen(s1) + strlen(s2) + 1); // +1 for the null-terminator
  // check for errors in malloc
  strcpy(result, s1);
  strcat(result, s2);
  return result;
}


string
replace(string a, int i, int len, char* rep) {
  a.replace(i, len, rep);
  return a;
}



string
traverse(uint32          idx,
         vector<uint32>  refIdxList,
         vector<uint32>  refLenList,
         map<int, vector<char*> >  posHaps,
         string          candidate,
         vector<int>     path,
         varMer*         seqMer) {

  if (idx < refIdxList.size()) {
    // fprintf(stderr, "[ DEBUG ] :: idx = %u | candidate = %s\n", idx, candidate.c_str());
    vector<char*> haps = posHaps.at(idx);
    uint32 refLen = refLenList.at(idx);    // Keep refLen so we revert in case there are > 1 ALTs

    //  for each haplotype sequence, make a combination
    //  j = 0 is always the REF allele
    for (int j = 0; j < haps.size(); j++) {
      path.push_back(j);
      char* hap = haps.at(j);
      string replaced = candidate;
      int delta = 0;
      int skipped = 0;
      bool overlaps = false;

      //  Make the replaced ver. Skip if reference allele.
      if ( j > 0) {
        refLenList.at(idx) = refLen;  // initialize in case it was modified by a previous j
        replaced = replace(candidate, refIdxList.at(idx), refLenList.at(idx), hap);
        //  fprintf(stderr, "[ DEBUG ] :: replaced = %s | refIdxList.at(%u) = %u | refLenList.at(%u) = %u, hap = %s\n", replaced.c_str(), idx, refIdxList.at(idx), idx, refLenList.at(idx), hap);

        //  Apply to the rest of the positions, after skipping overlaps
        //  refIdx in overlaps should remain as they were as we are using ref allele at these sites anyway
        //  This has to be done *before* skipping as idx changes.
        delta = strlen(hap) - refLenList.at(idx);

        //  Affected base right most index due to this replacement
        uint32 refAffected = refIdxList.at(idx) + refLenList.at(idx);

        //  Change refLen
        refLenList.at(idx) = strlen(hap);

        //  If the next variant overlaps, skip it
        if (idx + 1 < refIdxList.size()) {
          for (uint32 i = idx + 1; i < refIdxList.size(); i++) {
            // fprintf(stderr, "[ DEBUG ] :: Does %u overlap %u? hap = %d\n", refIdxList.at(i), refAffected, j);
            if ( refIdxList.at(i) < refAffected ) {
              //  fprintf(stderr, "[ DEBUG ] :: Yes, %u overlaps %u. skip %u.\n", refIdxList.at(i), refAffected, refIdxList.at(i));
              overlaps = true;
              idx++;             //  force skip by jumping the idx
              path.push_back(0); //  take the ref allele for 'no change'
              skipped++;
            } else {
              //  fprintf(stderr, "[ DEBUG ] :: No, %u DOES NOT overlaps %u. keep %u\n", refIdxList.at(i), refAffected, refIdxList.at(i));
              break;
            }
          }
        }

        //  if this is the end of the list, we are done
        if ( overlaps && idx == refIdxList.size() - 1) {
          seqMer->addSeqPath(replaced, path, refIdxList, refLenList);
          /******* DEBUG *********
                   fprintf(stderr, "[ DEBUG ] :: We are done with overlapped %u. replaced = %s : hap = %d , skipped = %d\n", refIdxList.at(idx), replaced.c_str(), j, skipped);
                   fprintf(stderr, "[ DEBUG ] ::  path (size=%lu) =", path.size());
                   for (int hh = 0; hh < path.size(); hh++) {
                   fprintf(stderr, "\t%d,", path.at(hh));
                   }
                   fprintf(stderr, "\n");
                   fprintf(stderr, "[ DEBUG ] ::  refIdxList (size=%lu) =", refIdxList.size());
                   for (int hh = 0; hh < refIdxList.size(); hh++) {  fprintf(stderr, "\t%d,", refIdxList.at(hh));  }
                   fprintf(stderr, "\n");
                   fprintf(stderr, "[ DEBUG ] ::  refLenList (size=%lu) =", refLenList.size());
                   for (int hh = 0; hh < refLenList.size(); hh++) {  fprintf(stderr, "\t%d,", refLenList.at(hh));  }
                   fprintf(stderr, "\n\n");
          ************************/
          for (int k = 0; k < skipped; k++) {
            path.pop_back();
            idx--;
            //  fprintf(stderr, "[ DEBUG ] :: Pop! idx = %u\n", idx);
          }
          path.pop_back();  // pop the current j
          continue;
        }

        for (uint32 i = idx + 1; i < refIdxList.size(); i++) {
          refIdxList.at(i) += delta;
          //  fprintf(stderr, "[ DEBUG ] :: Shift refIdxList.at(%u) + %d = %d \n", i, delta, refIdxList.at(i));
        }

      } // Done with what needs to be done with ALT

      //  Traverse
      replaced = traverse(idx + 1, refIdxList, refLenList, posHaps, replaced, path, seqMer);

      //  Only do something when all nodes are visited
      if (idx == refIdxList.size() - 1) {
        seqMer->addSeqPath(replaced, path, refIdxList, refLenList);
        /********* DEBUG **********
                   fprintf(stderr, "[ DEBUG ] :: We are done with %u. replaced = %s : hap = %d \n", refIdxList.at(idx), replaced.c_str() , j) ;
                   fprintf(stderr, "[ DEBUG ] ::  path (size=%lu) =", path.size());
                   for (int hh = 0; hh < path.size(); hh++) {
                   fprintf(stderr, "\t%d,", path.at(hh));
                   }
                   fprintf(stderr, "\n");
                   fprintf(stderr, "[ DEBUG ] ::  refIdxList (size=%lu) =", refIdxList.size());
                   for (int hh = 0; hh < refIdxList.size(); hh++) {    fprintf(stderr, "\t%d,", refIdxList.at(hh));    }
                   fprintf(stderr, "\n");
                   fprintf(stderr, "[ DEBUG ] ::  refLenList (size=%lu) =", refLenList.size());
                   for (int hh = 0; hh < refLenList.size(); hh++) {    fprintf(stderr, "\t%d,", refLenList.at(hh));    }
                   fprintf(stderr, "\n");
        ***************************/
      }

      //  fprintf(stderr, "[ DEBUG ] :: idx = %u\n", idx);
      if ( idx + 1 < refIdxList.size()) {
        for (uint32 i = idx + 1; i < refIdxList.size(); i++) {
          refIdxList.at(i) -= delta;
          //  fprintf(stderr, "[ DEBUG ] :: Shift back refIdxList.at(%u) - %d = %d \n", i, delta, refIdxList.at(i));
        }
      }

      if ( skipped > 0 ) {  // put the idx and skipped back
        for (int k = 0; k < skipped; k++) {
          //  fprintf(stderr, "[ DEBUG ] :: Pop! %d\n", k);
          path.pop_back();
          idx--;
        }
      }

      path.pop_back();  // pop the current j
      //  fprintf(stderr, "[ DEBUG ] :: End of idx = %u | j = %d\n\n", idx, j);
    }
  }
  idx--;
  //fprintf(stderr, "traverse() returns candidate idx %d '%s'\n", idx, candidate.c_str());
  return candidate;
}



void
processVariants(void *G, void *T, void *S) {
  merfinGlobal  *g = (merfinGlobal *)G;
  merfinThrData *t = (merfinThrData *)T;
  merfinInput   *s = (merfinInput *)S;

  fprintf(stderr, "Processing sequence %s for variants\n", s->seq.ident());

  //  If no variants for this sequence, ignore it.

  map<string, vector<posGT*>*> *mapChrPosGT = g->inVcf->_mapChrPosGT;

  if (mapChrPosGT->find(string(s->seq.ident())) == mapChrPosGT->end()) {
    fprintf(stderr, "No variants in vcf for contig '%s'. Skipping.\n", s->seq.ident());
    return;
  }

  //  Initialize the per-thread data if needed.

  if (t->oDebug == nullptr) {
    char  name[FILENAME_MAX+1];

    snprintf(name, FILENAME_MAX, "%s.%02d.debug", g->outName, t->threadID);
    t->oDebug = new compressedFileWriter(name);

    snprintf(name, FILENAME_MAX, "%s.%02d.vcf", g->outName, t->threadID);
    t->oVcf = new compressedFileWriter(name);

    for (string header : g->inVcf->getHeaders())
      fprintf(t->oVcf->file(), "%s\n", header.c_str());
  }

  //  Get chromosome specific posGTs, and iterate over.

  vector<posGT *> *posGTlist = mapChrPosGT->at(s->seq.ident());

  vector<uint32>           refIdxList;
  vector<uint32>           refLenList;
  vector<int>              path;
  map<int, vector<char*> > mapPosHap;

  for (uint64 posGtIdx = 0; posGtIdx < posGTlist->size(); posGtIdx++) {
    posGT               *posGt  = posGTlist->at(posGtIdx);
    uint32               rStart = posGt->_rStart;   //  0-based!
    uint32               rEnd   = posGt->_rEnd;     //  1-based!
    vector<gtAllele *>  *gts    = posGt->_gts;

    uint32               K_PADD = kmer::merSize() - 1;

    if (rStart > K_PADD)                   rStart -= K_PADD;
    else                                   rStart  = 0;

    if (rEnd < s->seq.length() - K_PADD)   rEnd   += K_PADD;
    else                                   rEnd    = s->seq.length();

    refIdxList.clear();
    refLenList.clear();
    path.clear();
    mapPosHap.clear();

    //  Debug report the mapPosHap

    //fprintf(stderr, "\n");
    //fprintf(stderr, "[ DEBUG ] :: %s : %u - %u\n", s->seq.ident(), rStart, rEnd);
    //fprintf(stderr, "[ DEBUG ] :: gts->size = %lu | ", gts->size());
    //for (uint32 i = 0; i < gts->size(); i++)
    //  fprintf(stderr, "gt->_pos = %u ",  gts->at(i)->_pos);
    //fprintf(stderr, "\n");

    //  Load mapPosHap

    for (uint32 i = 0; i < gts->size(); i++) {
      gtAllele *gt = gts->at(i);

      refIdxList.push_back(gt->_pos - rStart);
      refLenList.push_back(gt->_refLen);

      //  add alleles. alleles.at(0) is always the ref allele
      mapPosHap.insert(pair<int, vector<char*> >(i, *(gt->alleles)));
    }

    //  Load original sequence from rStart to rEnd.

    char *refTemplate = new char[kmer::merSize() * 2 + rEnd - rStart];

    if (s->seq.copy(refTemplate, rStart, rEnd, true ) == false) {
      fprintf(stderr, "PANIC : Invalid region specified: %s : %u - %u\n", s->seq.ident(), rStart, rEnd);
      continue;
    }

    if (refIdxList.size() > g->comb) {
      fprintf(stderr, "PANIC : Combination %s:%u-%u has too many variants ( found %lu > %u ) to evaluate. Consider filtering the vcf upfront. Skipping...\n",
              s->seq.ident(), rStart, rEnd, gts->size(), g->comb);
      continue;
    }

    //fprintf(stderr, "%s\n", refTemplate);

    //

    varMer* seqMer = new varMer(posGt);

    //  Traverse through each gt combination
    //  fprintf(stderr, "[ DEBUG ] :: traverse begin\n");
    traverse(0, refIdxList, refLenList, mapPosHap, refTemplate, path, seqMer);
    //  fprintf(stderr, "[ DEBUG ] :: traverse done\n");

    //  Score each combination
    //  fprintf(stderr, "[ DEBUG ] :: score begin\n");
    seqMer->score(g);
    //  fprintf(stderr, "[ DEBUG ] :: score completed\n");

    //  Store the avgK of the reference to compute the delta
    //RefAvgK = seqMer->getAvgAbsK(0);

    //  Save debug info.

    for (uint64 idx = 0; idx < seqMer->seqs.size(); idx++) {
      fprintf(t->oDebug->file(), "%lu\t%s:%u-%u\t%s\t%u\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t",
              t->varMerId++,
              s->seq.ident(),
              rStart,
              rEnd,
              seqMer->seqs.at(idx).c_str(),  //  seq
              seqMer->numMs.at(idx),         //  missing
              seqMer->getMinAbsK(idx),
              seqMer->getMaxAbsK(idx),
              seqMer->getMedAbsK(idx),
              seqMer->getAvgAbsK(idx),
              //seqMer->getAvgAbsdK(idx, RefAvgK),
              seqMer->getTotdK(idx));

      //  new vcf records
      //  fprintf(stderr, "%s:%u-%u seqMer->gtPaths.at(idx).size() %d\n", s->seq.ident(), rStart, rEnd, seqMer->gtPaths.at(idx).size());

      if  ( seqMer->gtPaths.at(idx).size() > 0 ) {
        for (uint64 i = 0; i < seqMer->gtPaths.at(idx).size(); i++) {
          // Ignore the ref-allele (0/0) GTs
          // print only the non-ref allele variants for fixing
          int altIdx = seqMer->gtPaths.at(idx).at(i);
          if (altIdx > 0) {
            fprintf(t->oDebug->file(), "%s %u . %s %s . PASS . GT 1/1  ",
                    s->seq.ident(),
                    (gts->at(i)->_pos+1),
                    gts->at(i)->alleles->at(0),
                    gts->at(i)->alleles->at(altIdx)
                    );
          }
        }
      }

      fprintf(t->oDebug->file(), "\n");
    }

    //  Generate output VCFs.

    // Experimental: output vcf according to k*
    if (g->bykstar) {
      fprintf(t->oVcf->file(), "%s", seqMer->bestVariant().c_str());
    }

    // Filter vcf and print as it was in the original vcf, conservatively
    else {
      vector<vcfRecord*> records = seqMer->bestVariantOriginalVCF();

      for (uint64 i = 0; i < records.size(); i++)
        records.at(i)->save(t->oVcf);
    }

    //  Cleanup.

    delete   seqMer;
    delete[] refTemplate;
  }  //  Over posGTlist.
}






void
outputVariants(void *G, void *S) {
  merfinGlobal  *g = (merfinGlobal *)G;
  merfinInput   *s = (merfinInput *)S;

  fprintf(stderr, "Output sequence %s\n", s->seq.ident());

  //  Open the output file and write headers.

  if (g->oVCF == nullptr) {
    g->oVCF = new compressedFileWriter(concat(g->outName, ".polish.vcf"));

    for (string header : g->inVcf->getHeaders())
      fprintf(g->oVCF->file(), "%s\n", header.c_str());
  }

  //



  //  Bye.

  delete s;
}
