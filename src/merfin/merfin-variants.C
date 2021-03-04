
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



string
traverse(uint32                          idx,
         vector<uint32>                 &refIdxList,   //  This CAN be a reference.
         vector<uint32>                  refLenList,   //  This CANNOT be a reference.
         map<int, vector<char const*> >  posHaps,
         string                          candidate,    //  This cannot be a ref; won't compile
         vector<int>                    &path,         //  This can be a reference
         varMer                         *seqMer) {

  //fprintf(stderr, "enter   traverse(%2d) - IdxList ", idx);  for (uint32 ii=0; ii<refIdxList.size(); ii++)  fprintf(stderr, " %4u", refIdxList[ii]);  fprintf(stderr, "\n");

  assert(idx < refIdxList.size());

  vector<char const *> &haps   = posHaps[idx];
  uint32                refLen = refLenList[idx];    // Keep refLen so we revert in case there are > 1 ALTs

  //  for each haplotype sequence, make a combination
  //  j = 0 is always the REF allele
  for (int j = 0; j < haps.size(); j++) {
    path.push_back(j);

    char const * hap       = haps[j];
    string       replaced  = candidate;
    int          skipped   = 0;
    bool         overlaps  = false;
    int          delta     = 0;

    //  Make the replaced ver. Skip if reference allele.
    if ( j > 0) {
      refLenList[idx] = refLen;  // initialize in case it was modified by a previous j

      replaced = candidate;
      replaced.replace(refIdxList[idx], refLenList[idx], hap);

      //fprintf(stderr, "REPLACE candidate '%s' ->\n", candidate.c_str());
      //fprintf(stderr, "REPLACE replaced  '%s' by change %u-%u to '%s'\n", replaced.c_str(), refIdxList[idx], refLenList[idx], hap);

      //  Apply to the rest of the positions, after skipping overlaps
      //  refIdx in overlaps should remain as they were as we are using ref allele at these sites anyway
      //  This has to be done *before* skipping as idx changes.
      delta = strlen(hap) - refLenList[idx];

      //  Affected base right most index due to this replacement
      uint32 refAffected = refIdxList[idx] + refLenList[idx];

      //  Change refLen
      refLenList[idx] = strlen(hap);

      //  If the next variant overlaps, skip it
      for (uint32 i = idx + 1; i < refIdxList.size(); i++) {
        if (refIdxList[i] >= refAffected)
          break;

        overlaps = true;
        idx++;             //  force skip by jumping the idx
        path.push_back(0); //  take the ref allele for 'no change'
        skipped++;
      }

      //  if this is the end of the list, we are done
      if ( overlaps && idx == refIdxList.size() - 1) {
        seqMer->addSeqPath(replaced, path, refIdxList, refLenList);

        for (int k = 0; k < skipped; k++) {
          path.pop_back();
          idx--;
        }

        path.pop_back();  // pop the current j
        continue;
      }

      for (uint32 i = idx + 1; i < refIdxList.size(); i++)
        refIdxList[i] += delta;
    } // j>0 - Done with what needs to be done with ALT

    //fprintf(stderr, "goto    traverse(%2d) - IdxList ", idx+1);  for (uint32 ii=0; ii<refIdxList.size(); ii++)  fprintf(stderr, " %4u", refIdxList[ii]);  fprintf(stderr, "\n");

    if (idx+1 < refIdxList.size())
      replaced = traverse(idx + 1, refIdxList, refLenList, posHaps, replaced, path, seqMer);
    //else
    //  fprintf(stderr, "        traverse(%2d) - at the end\n", idx);

    //fprintf(stderr, "back in traverse(%2d) - IdxList ", idx+0);  for (uint32 ii=0; ii<refIdxList.size(); ii++)  fprintf(stderr, " %4u", refIdxList[ii]);  fprintf(stderr, "\n");

    //  Only do something when all nodes are visited
    if (idx == refIdxList.size() - 1)
      seqMer->addSeqPath(replaced, path, refIdxList, refLenList);


    for (uint32 i = idx + 1; i < refIdxList.size(); i++)
      refIdxList[i] -= delta;

    for (int k = 0; k < skipped; k++) {
      path.pop_back();
      idx--;
    }

    path.pop_back();  // pop the current j
  }

  //fprintf(stderr, "        traverse(%2d) - finished\n", idx);
  return(candidate);
}




void
processVariants(void *G, void *T, void *S) {
  merfinGlobal  *g = (merfinGlobal *)G;
  merfinThrData *t = (merfinThrData *)T;
  merfinInput   *s = (merfinInput *)S;

  fprintf(stderr, "Processing sequence %s for variants\n", s->seq.ident());

  //  If no variants for this sequence, ignore it.

  map<string, vector<posGT*>*> &mapChrPosGT = g->inVcf->_mapChrPosGT;

  if (mapChrPosGT.find(string(s->seq.ident())) == mapChrPosGT.end()) {
    fprintf(stderr, "No variants in vcf for contig '%s'. Skipping.\n", s->seq.ident());
    return;
  }

  //  Initialize the per-thread data if needed.

  if (t->oDebug == nullptr) {
    char  name[FILENAME_MAX+1];

    snprintf(name, FILENAME_MAX, "%s.%02d.debug.gz", g->outName, t->threadID);
    t->oDebug = new compressedFileWriter(name);
  }

  //  Get chromosome specific posGTs, and iterate over.

  vector<posGT *> *posGTlist = mapChrPosGT[s->seq.ident()];

  vector<uint32>           refIdxList;
  vector<uint32>           refLenList;
  vector<int>              path;
  map<int, vector<char const *> > mapPosHap;

  for (uint64 posGtIdx = 0; posGtIdx < posGTlist->size(); posGtIdx++) {
    posGT               *posGt  = posGTlist->at(posGtIdx);
    uint32               rStart = posGt->_rStart;   //  0-based!
    uint32               rEnd   = posGt->_rEnd;     //  1-based!
    vector<gtAllele *>  &gts    = posGt->_gts;

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
    //fprintf(stderr, "[ DEBUG ] :: gts.size = %lu | ", gts.size());
    //for (uint32 i = 0; i < gts.size(); i++)
    //  fprintf(stderr, "gt->_pos = %u ",  gts[i]->_pos);
    //fprintf(stderr, "\n");

    //  Load mapPosHap

    for (uint32 i = 0; i < gts.size(); i++) {
      gtAllele *gt = gts[i];

      refIdxList.push_back(gt->_pos - rStart);
      refLenList.push_back(gt->_refLen);

      //  add alleles. alleles[0] is always the ref allele
#warning is this copying the whole vector?
      mapPosHap[i] = gt->_alleles;
    }

    //  Load original sequence from rStart to rEnd.

    char *refTemplate = new char[kmer::merSize() * 2 + rEnd - rStart];

    if (s->seq.copy(refTemplate, rStart, rEnd, true ) == false) {
      fprintf(stderr, "PANIC : Invalid region specified: %s : %u - %u\n", s->seq.ident(), rStart, rEnd);
      continue;
    }

    if (refIdxList.size() > g->comb) {
      fprintf(stderr, "PANIC : Combination %s:%u-%u has too many variants ( found %lu > %u ) to evaluate. Consider filtering the vcf upfront. Skipping...\n",
              s->seq.ident(), rStart, rEnd, gts.size(), g->comb);
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
              seqMer->seqs[idx].c_str(),  //  seq
              seqMer->numMs[idx],         //  missing
              seqMer->getMinAbsK(idx),
              seqMer->getMaxAbsK(idx),
              seqMer->getMedAbsK(idx),
              seqMer->getAvgAbsK(idx),
              //seqMer->getAvgAbsdK(idx, RefAvgK),
              seqMer->getTotdK(idx));

      //  new vcf records
      //  fprintf(stderr, "%s:%u-%u seqMer->gtPaths[idx].size() %d\n", s->seq.ident(), rStart, rEnd, seqMer->gtPaths[idx].size());

      if  ( seqMer->gtPaths[idx].size() > 0 ) {
        for (uint64 i = 0; i < seqMer->gtPaths[idx].size(); i++) {
          // Ignore the ref-allele (0/0) GTs
          // print only the non-ref allele variants for fixing
          int altIdx = seqMer->gtPaths[idx][i];
          if (altIdx > 0) {
            fprintf(t->oDebug->file(), "%s %u . %s %s . PASS . GT 1/1  ",
                    s->seq.ident(),
                    gts[i]->_pos+1,
                    gts[i]->_alleles[0],
                    gts[i]->_alleles[altIdx]
                    );
          }
        }
      }

      fprintf(t->oDebug->file(), "\n");
    }

    //  Generate output VCFs.

    // Experimental: output vcf according to k*
    if (g->bykstar) {
      //fprintf(t->oVcf->file(), "%s", seqMer->bestVariant().c_str());
      s->result += seqMer->bestVariant();
    }

    // Filter vcf and print as it was in the original vcf, conservatively
    else {
      vector<vcfRecord*> records = seqMer->bestVariantOriginalVCF();

      for (uint64 i = 0; i < records.size(); i++)
        //records[i]->save(t->oVcf);
        s->result += records[i]->save();
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
    char  name[FILENAME_MAX+1];

    snprintf(name, FILENAME_MAX, "%s.polish.vcf", g->outName);

    g->oVCF = new compressedFileWriter(name);

    for (string header : g->inVcf->getHeaders())
      fprintf(g->oVCF->file(), "%s\n", header.c_str());
  }

  //  Dump the output string to the output file.  Flush it just
  //  so we can get a reliable record of where we are.

  fputs(s->result.c_str(), g->oVCF->file());
  fflush(g->oVCF->file());

  //  Bye.

  delete s;
}
