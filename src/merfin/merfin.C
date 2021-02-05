
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

#include "kmers.H"
#include "system.H"
#include "sequence.H"
#include "bits.H"
#include "strings.H"
#include "files.H"
#include "vcf.H"
#include "kmetric.H"
#include "varMer.H"
#include "types.H"

#include <vector>
#include <map>
#include <cmath>

// to read the lookup table
#include <iostream>
#include <fstream>
#include <sstream>

#define OP_NONE       0
#define OP_HIST       1
#define OP_DUMP       2
#define OP_VAR_MER    3

func_t getreadK;

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
    uint32 refLen = refLenList.at(idx);	// Keep refLen so we revert in case there are > 1 ALTs

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
          for (int hh = 0; hh < refIdxList.size(); hh++) {	fprintf(stderr, "\t%d,", refIdxList.at(hh));	}
          fprintf(stderr, "\n");
          fprintf(stderr, "[ DEBUG ] ::  refLenList (size=%lu) =", refLenList.size());
          for (int hh = 0; hh < refLenList.size(); hh++) {	fprintf(stderr, "\t%d,", refLenList.at(hh));	}
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
  return candidate;
}

void
dumpKmetric(char               *outName,
    			  char      			   *seqName,
            dnaSeqFile         *sfile,
            merylExactLookup   *rlookup,
            merylExactLookup   *alookup,
            bool                skipMissings,
            vector<string>      copyKmerDict,
            int                 threads) {

  double prob;
  uint64 tot_missing = 0;
  map<string, uint64> order;
  string tmp = tmpnam(nullptr);

  compressedFileWriter *k_dump   = new compressedFileWriter(outName);

  fprintf(stderr, "\nGenerating fasta index.\n");  
  sfile->generateIndex();

  int ctgn = sfile->numberOfSequences();
  
  //use at most threads equal to the number of sequences
  if (threads > ctgn) {
  	threads = ctgn;
  }
  
  fprintf(stderr, "\nNumber of contigs: %u\n", ctgn);
    
#pragma omp parallel for reduction(+:tot_missing) num_threads(threads) schedule(dynamic)
  for (uint64 seqId=0; seqId<ctgn;seqId++)
  {
    dnaSeq seq;
    double asmK;
    double readK;
    double kMetric;
    uint64 fValue = 0;
    uint64 rValue = 0;
#pragma omp critical
    {
      sfile->loadSequence(seq);      
    }

    kmerIterator kiter(seq.bases(), seq.length());
    uint64 missing = 0;
    uint64 tot = 0;

    char filename[64];

    FILE *out;
    sprintf(filename, "%s_%lu.dump", tmp.c_str(), seqId);

    auto p = order.insert(pair<string, uint64>(seq.ident(), seqId));

    if (!p.second)
    {
      fprintf(stderr, "\nSequence name used twice: %s\nPlease use only unique names.\n", seq.ident());
      exit (-1);	  
    }

    out = fopen(filename, "w");

    while (kiter.nextBase()) {
      if (kiter.isValid() == true) {
        tot++;
        getK(rlookup, alookup, kiter.fmer(), kiter.rmer(), copyKmerDict, readK, asmK, prob);
        kMetric = getKmetric(readK, asmK);
        if ( readK == 0 ){
          missing++;
        }

        if ( skipMissings )  continue;

        fprintf(out, "%s\t%lu\t%.2f\t%.2f\t%.2f\n",
            seq.ident(),
            kiter.position(),
            readK,
            asmK,
            kMetric
            );
      }
    }
    fclose(out);
#pragma omp critical
    {
      tot_missing+=missing;
      fprintf(stderr, "%s\t%lu\t%lu\t%lu\n",
          seq.ident(),
          missing,
          tot_missing,
          tot
          );
    }
  }

  ofstream ofile(outName, ios::out | ios::app); 
  char filename[64];

  sfile = new dnaSeqFile(seqName);
  dnaSeq seq;
  for (uint32 seqId=0; sfile->loadSequence(seq); seqId++) {
    sprintf(filename, "%s_%lu.dump", tmp.c_str(), order.at(seq.ident()));
    order.erase (seq.ident());
    ifstream ifile(filename, ios::in); 
    ofile << ifile.rdbuf(); 
    remove(filename);
  }
}

void
histKmetric(char               *outName,
            char		      	   *seqName,
            dnaSeqFile         *sfile,
            merylExactLookup   *rlookup,
            merylExactLookup   *alookup,
            vector<string>      copyKmerDict,
            int                 threads) {

  double prob = 1;
  uint64 tot_missing = 0;
  uint64 tot_kasm = 0;
  uint32 ksize = kmer::merSize();
  
  //  compressedFileWriter *k_values = new compressedFileWriter(concat(outName, ".gz"));
  compressedFileWriter *k_hist   = new compressedFileWriter(outName);


  //  variables for generating histogram
  uint64   histMax  = 32 * 1024 * 1024;
  uint64 * overHist = new uint64[histMax];	// positive k* values, overHist[0] = bin 0.0 ~ 0.2
  uint64 * undrHist = new uint64[histMax];	// negative k* values
  double   roundedReadK = 0;
  double   overcpy  = 0;

  for (uint64 ii = 0; ii < histMax; ii++) {
    overHist[ii] = 0;
    undrHist[ii] = 0;
  }

  fprintf(stderr, "\nGenerating fasta index.\n");  
  sfile->generateIndex();

  int ctgn = sfile->numberOfSequences();
  
  //use at most threads equal to the number of sequences
  if (threads > ctgn) {
  	threads = ctgn;
  }
  
  fprintf(stderr, "\nNumber of contigs: %u\n", ctgn);

#pragma omp parallel for reduction (+:overcpy) num_threads(threads) schedule(dynamic)
    for (uint32 seqId=0; seqId<ctgn;seqId++)
    {

    dnaSeq seq;
    double asmK;
    double readK;
    double kMetric;
    uint64 fValue = 0;
    uint64 rValue = 0;

#pragma omp critical
      {
        sfile->loadSequence(seq);
      }

      kmerIterator kiter(seq.bases(), seq.length());
      uint64 missing = 0;
      uint64 kasm = 0;
      double err;
      double qv;
      uint64 * undrHist_pvt = new uint64[histMax];
      uint64 * overHist_pvt = new uint64[histMax];

      while (kiter.nextBase()) {
        if (kiter.isValid() == true) {
          kasm++;
          getK(rlookup, alookup, kiter.fmer(), kiter.rmer(), copyKmerDict, readK, asmK, prob);
          kMetric = getKmetric(readK, asmK);

          if ( readK == 0 ) {
            missing++;
          } else if ( readK < asmK ) {
            undrHist_pvt[(uint64) (((-1 * kMetric) + 0.1) / 0.2)]++;
            //  TODO: Check if this kmer was already counted. Only if not,
            //  overcpy += (asmK - readK)
            overcpy += (double) (1 - readK / asmK) * prob;  //  (asmK - readK) / asmK
          } else { // readK > asmK
            overHist_pvt[(uint64) ((kMetric + 0.1 ) / 0.2)]++;
          }
        }
      }

#pragma omp critical
      {
        tot_missing+=missing;
        tot_kasm+=kasm;

        for(uint64 ii = histMax - 1; ii > 0; ii--) {
          undrHist[ii] += undrHist_pvt[ii];
        }
        undrHist[0] += undrHist_pvt[0];
        overHist[0] += overHist_pvt[0];
        for (uint64 ii = 1; ii < histMax; ii++) {
          overHist[ii] += overHist_pvt[ii];
        }

        err = 1 - pow((1-((double) missing) / kasm), (double) 1/ksize);
        qv = -10*log10(err);

        fprintf(stderr, "t%s\t%lu\t%lu\t%lu\t%.2f\n",
            seq.ident(),
            missing,
            tot_missing,
            kasm,
            qv
            );
      }
    }

  for (uint64 ii = histMax - 1; ii > 0; ii--) {
    if (undrHist[ii] > 0)  fprintf(k_hist->file(), "%.1f\t%lu\n", ((double) ii * -0.2), undrHist[ii]);
  }
  fprintf(k_hist->file(), "%.1f\t%lu\n", 0.0, (undrHist[0]+overHist[0]));
  for (uint64 ii = 1; ii < histMax; ii++) {
    if (overHist[ii] > 0)  fprintf(k_hist->file(), "%.1f\t%lu\n", ((double) ii * 0.2), overHist[ii]);
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "K-mers not found in reads (missing) : %lu\n", tot_missing);
  fprintf(stderr, "K-mers overly represented in assembly: %.2f\n", overcpy);
  fprintf(stderr, "K-mers found in the assembly: %lu\n", tot_kasm);
  double err = 1 - pow((1-((double) tot_missing) / tot_kasm), (double) 1/ksize);
  double qv = -10*log10(err);
  fprintf(stderr, "Merqury QV: %.2f\n", qv);
  tot_missing += (uint64) ceil(overcpy);
  err = 1 - pow((1-((double) tot_missing) / tot_kasm), (double) 1/ksize);
  qv = -10*log10(err);
  fprintf(stderr, "Adjusted QV: %.2f\n", qv);
  fprintf(stderr, "*** Note this QV is valid only if -seqmer was generated with -sequence ***\n\n");
  fprintf(stderr, "*** Merqury QV only considers missing kmers as errors. Merfin QV includes overrepresented kmers. ***\n\n");
}


void
varMers(char			       *dnaSeqName,
        dnaSeqFile       *sfile,
        vcfFile          *vfile,
        merylExactLookup *rlookup,
        merylExactLookup *alookup,
        char             *out,
        uint32			      comb,
        bool			        nosplit,
        vector<string>    copyKmerDict,
        bool              bykstar,
        int				        threads) {

  //  output file
  ofstream oVcf(concat(out, ".polish.vcf"));
  ofstream oDebug(concat(out, ".debug"));

  //  write vcf headers
  for ( string header : vfile->getHeaders()) {
    oVcf << header << "\n";
  }

  //  What is the kmer size?
  uint32 ksize = kmer::merSize();

  //  Merge posGTlist for each chr within ksize
  fprintf(stderr, "Merge variants within %u-mer bases, splitting combinations greater than %u.\n", ksize, comb);
  vfile->mergeChrPosGT(ksize, comb, nosplit);

  // print CHR rStart rEnd POS HAP1 HAP2 minHAP1 minHAP2 to out.debug
  
  map<string, vector<posGT*>*> *mapChrPosGT = vfile->_mapChrPosGT;

  uint32   pos;
  uint32   K_PADD = ksize - 1;
  
  fprintf(stderr, "\nGenerating fasta index.\n");  
  sfile->generateIndex();

  uint32 ctgn = sfile->numberOfSequences();

  //  threads = omp_get_max_threads();
  uint32 chunks = ctgn;

  //  Run chunks at maximum of MAX_FILES. Setting this too high will break network storag
  const int64 openMax = sysconf(_SC_OPEN_MAX) / 2; // prevent too many open files

  if (chunks > openMax) { chunks = openMax;  }
  uint32 chunkLeft = chunks;
  string tmpPrefix = string("tmp_") + sfile->filename() + string("_");
  splitToWords fp(out, splitPaths);
  tmpPrefix += fp.last() + string("_");

  fprintf(stderr, "\nTotal of %u sequences found. Will be processed over %d threads, with maximum %u intermediate %s# .debug and .vcf files\n",
      ctgn, threads, chunks, tmpPrefix.c_str());

  for (uint32 ii = 0; ii < ctgn; ii += chunks) {

    if ( ii + chunks > ctgn ) chunkLeft = ctgn - ii;
    fprintf(stderr, "Reading sequence %u - %u  ... \n", ii, ii+chunkLeft);
    
#pragma omp parallel for num_threads(threads) schedule(dynamic)
    for (uint32 seqId = ii; seqId < ii + chunkLeft; seqId++) {
    //  for (uint32 seqId=0; seqId<ctgn; seqId++)  {
      dnaSeq   seq;
      uint64   varMerId = 0;
#pragma omp critical
      {
        sfile->findSequence(seqId);
        sfile->loadSequence(seq);
        //  printf(stderr, "Load\tseqID: %u\tseq: %s\n", seqId, seq.ident()); 
      }
    
      // in case this seq has no variants, ignore this seq
      if (mapChrPosGT->find(string(seq.ident())) == mapChrPosGT->end()) {
        fprintf(stderr, "No variants found in vcf for contig \'%s\'. Skipping.\n", seq.ident());
        continue;
      }

      //  for each seqId
      fprintf(stderr, "Processing \'%s\'\n", seq.ident());
      
      //  get chr specific posGTs
      vector<posGT*>  *posGTlist = mapChrPosGT->at(seq.ident());

      //  temporary output files
      string tmpDebugName = tmpPrefix + to_string(seqId) + ".debug";
      compressedFileWriter  *tmpDebug  = new compressedFileWriter(tmpDebugName.c_str());

      string tmpVcfName = tmpPrefix + to_string(seqId) + ".vcf";
      compressedFileWriter  *tmpVcf    = new compressedFileWriter(tmpVcfName.c_str());

      //  get sequence combinations on each posGT list
      for (uint64 posGtIdx = 0; posGtIdx < posGTlist->size(); posGtIdx++) {

        // initialize variables
        posGT  *posGt = posGTlist->at(posGtIdx);
        uint32 rStart = posGt->_rStart;  // 0-based
        if (rStart > K_PADD) { rStart -= K_PADD; }
        else { rStart = 0; }

        uint32 rEnd   = posGt->_rEnd;             // 1-based
        if (rEnd < seq.length() - K_PADD) {  rEnd += K_PADD;  }
        else { rEnd = seq.length(); }

        //  fprintf(stderr, "\n[ DEBUG ] :: %s : %u - %u\n", seq.ident(), rStart, rEnd);

        vector<gtAllele*> *gts = posGt->_gts;
        vector<uint32> refIdxList;
        vector<uint32> refLenList;
        map<int, vector<char*> > mapPosHap;
        vector<int>     path;

         
        //  load mapPosHap
        //  fprintf(stderr, "[ DEBUG ] :: gts->size = %lu | ", gts->size());
        for (uint32 i = 0; i < gts->size(); i++) {
           gtAllele *gt = gts->at(i);
           refIdxList.push_back(gt->_pos - rStart);
           refLenList.push_back(gt->_refLen);

           //  fprintf(stderr, "gt->_pos = %u ",  gt->_pos);
           //  add alleles. alleles.at(0) is always the ref allele
           mapPosHap.insert(pair<int, vector<char*> >(i, *(gt->alleles)));
        }
        //  fprintf(stderr, "\n");
        
        // temporary sequence to hold ref bases
        char* refTemplate = new char[(ksize*2+rEnd-rStart)];      
        
        //  load original sequence from rStart to rEnd
        if ( ! seq.copy(refTemplate, rStart, rEnd, true )) {
          fprintf(stderr, "Invalid region specified: %s : %u - %u\n", seq.ident(), rStart, rEnd);
          continue;
        }
        // DEBUG			fprintf(stderr, "%s\n", refTemplate);

        if ( refIdxList.size() > comb ) {
          fprintf(stderr, "PANIC : Combination %s:%u-%u has too many variants ( found %lu > %u ) to evaluate. Consider filtering the vcf upfront. Skipping...\n", seq.ident(), rStart, rEnd, gts->size(), comb);
          continue;
        }

        varMer* seqMer = new varMer(posGt);

        //  traverse through each gt combination
        //  fprintf(stderr, "[ DEBUG ] :: traverse begin\n");
        traverse(0, refIdxList, refLenList, mapPosHap, refTemplate, path, seqMer);
        //  fprintf(stderr, "[ DEBUG ] :: traverse done\n");

        //  score each combination
        //  fprintf(stderr, "[ DEBUG ] :: score begin\n");
        seqMer->score(rlookup, alookup, copyKmerDict);
        //  fprintf(stderr, "[ DEBUG ] :: score completed\n");
        
        // store the avgK of the reference to compute the delta
        //  double refAvgK = seqMer->getAvgAbsK(0);

        //  print to debug
        for (uint64 idx = 0; idx < seqMer->seqs.size(); idx++) {
          fprintf(tmpDebug->file(), "%lu\t%s:%u-%u\t%s\t%u\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t",
            varMerId++,
            seq.ident(),
            rStart,
            rEnd,
            seqMer->seqs.at(idx).c_str(), //  seq
            seqMer->numMs.at(idx),	//  missing
            seqMer->getMinAbsK(idx),
            seqMer->getMaxAbsK(idx),
            seqMer->getMedAbsK(idx),
            seqMer->getAvgAbsK(idx),
            //seqMer->getAvgAbsdK(idx, refAvgK),
            seqMer->getTotdK(idx)
          );

          //  new vcf records
          //  fprintf(stderr, "%s:%u-%u seqMer->gtPaths.at(idx).size() %d\n", seq.ident(), rStart, rEnd, seqMer->gtPaths.at(idx).size());
          if  ( seqMer->gtPaths.at(idx).size() > 0 ) {
            for (uint64 i = 0; i < seqMer->gtPaths.at(idx).size(); i++) {
              // Ignore the ref-allele (0/0) GTs
              // print only the non-ref allele variants for fixing
              int altIdx = seqMer->gtPaths.at(idx).at(i);
              if (altIdx > 0) {
                fprintf(tmpDebug->file(), "%s %u . %s %s . PASS . GT 1/1  ",
                  seq.ident(),
                  (gts->at(i)->_pos+1),
                  gts->at(i)->alleles->at(0),
                  gts->at(i)->alleles->at(altIdx)
                );
              }
            }
          }
          fprintf(tmpDebug->file(), "\n");
          fflush(tmpDebug->file());
        }

        // generate vcfs
        if (bykstar) {
          // Experimental: output vcf according to k*
          fprintf(tmpVcf->file(), "%s", seqMer->bestVariant().c_str());
          fflush(tmpVcf->file());
        } else {
          // Filter vcf and print as it was in the original vcf, conservatively
          vector<vcfRecord*> records = seqMer->bestVariantOriginalVCF();
          if (records.size() > 0) {
            for (uint64 i = 0; i < records.size(); i++) {
              records.at(i)->save(tmpVcf);
              fflush(tmpVcf->file());
            }
          }
        }
     
        delete seqMer;
        delete[] refTemplate;

      } // end of for loop per posGtIdx
 
      delete tmpDebug;
      delete tmpVcf;

#pragma omp critical
      {
        //  fprintf(stderr, "[DEBUG] :: Writing output from %d - %s : %s\n", seqId, seq.ident(), tmpDebugName.c_str());
        ifstream tmpFileI(tmpDebugName, std::ios_base::binary);
        if (tmpFileI.is_open()) {
          oDebug << tmpFileI.rdbuf();
          tmpFileI.close();
          remove(tmpDebugName.c_str());
        } else {
          fprintf(stderr, "Unable to open file %s\n", tmpDebugName.c_str());
        }

        tmpFileI.open(tmpVcfName, std::ios_base::binary);
        if (tmpFileI.is_open()) {
          oVcf << tmpFileI.rdbuf();
          tmpFileI.close();
          remove(tmpVcfName.c_str());
        } else {
          fprintf(stderr, "Unable to open file %s\n", tmpVcfName.c_str());
        }
      }
    } // end of parallel for loop per sequence
  }
}

int
main(int argc, char **argv) {
  char           *seqName       = NULL;
  char           *vcfName       = NULL;
  char           *outName       = NULL;
  char           *seqDBname     = NULL;
  char           *readDBname    = NULL;
  char           *pLookupTable  = NULL;

  uint64          minV        = 0;
  uint64          maxV        = UINT64_MAX;
  static uint64   ipeak       = 0;
  bool            skipMissing = false;
  bool            nosplit     = false;
  bool            bykstar     = true;
  uint32          threads     = omp_get_max_threads();
  uint32          memory1     = 0;
  uint32          memory2     = 0;
  uint32          reportType  = OP_NONE;
  uint32		  comb = 15;

  vector<const char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-sequence") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-seqmers") == 0) {
      seqDBname = argv[++arg];

    } else if (strcmp(argv[arg], "-readmers") == 0) {
      readDBname = argv[++arg];

    } else if (strcmp(argv[arg], "-peak") == 0) {
      ipeak = strtouint64(argv[++arg]);
      
    } else if (strcmp(argv[arg], "-lookup") == 0) {

      pLookupTable = argv[++arg];

    } else if (strcmp(argv[arg], "-vcf") == 0) {
      vcfName = argv[++arg];

    } else if (strcmp(argv[arg], "-output") == 0) {
      outName = argv[++arg];

    } else if (strcmp(argv[arg], "-min") == 0) {
      minV = strtouint64(argv[++arg]);

    } else if (strcmp(argv[arg], "-max") == 0) {
      maxV = strtouint64(argv[++arg]);

    } else if (strcmp(argv[arg], "-threads") == 0) {
      threads = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-memory1") == 0) {
      memory1 = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-memory2") == 0) {
      memory2 = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-nosplit") == 0) {
      nosplit = true;

    } else if (strcmp(argv[arg], "-disable-kstar") == 0) {
      bykstar = false;

    } else if (strcmp(argv[arg], "-hist") == 0) {
      reportType = OP_HIST;

    } else if (strcmp(argv[arg], "-dump") == 0) {
      reportType = OP_DUMP;

    } else if (strcmp(argv[arg], "-skipMissing") == 0) {
      skipMissing = true;

    } else if (strcmp(argv[arg], "-vmer") == 0) {
      reportType = OP_VAR_MER;

    } else if (strcmp(argv[arg], "-comb") == 0) {
      comb = strtouint32(argv[++arg]);

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (seqName == NULL)
    err.push_back("No input sequences (-sequence) supplied.\n");
  if (seqDBname == NULL)
    err.push_back("No sequence meryl database (-seqmers) supplied.\n");
  if (readDBname == NULL)
    err.push_back("No read meryl database (-readmers) supplied.\n");
  if (ipeak == 0)
    err.push_back("Peak=0 or no haploid peak (-peak) supplied.\n");
  if (outName == NULL)
    err.push_back("No output (-output) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s <report-type>            \\\n", argv[0]);
    fprintf(stderr, "         -sequence <seq.fasta>     \\\n");
    fprintf(stderr, "         -seqmers  <seq.meryl>     \\\n");
    fprintf(stderr, "         -readmers <read.meryl>    \\\n");
    fprintf(stderr, "         -peak     <haploid_peak>  \\\n");
    fprintf(stderr, "         -lookup   <lookup_table>  \\\n");
    fprintf(stderr, "         -vcf      <input.vcf>     \\\n");
    fprintf(stderr, "         -output   <output>        \n\n");
    fprintf(stderr, "  Predict the kmer consequences of variant calls <input.vcf> given the consensus sequence <seq.fasta>\n");
    fprintf(stderr, "  and lookup the k-mer multiplicity in the consensus sequence <seq.meryl> and in the reads <read.meryl>.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Input -sequence and -vcf files can be FASTA or FASTQ; uncompressed, gz, bz2 or xz compressed\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Each input database can be filtered by value.  More advanced filtering\n");
    fprintf(stderr, "  requires a new database to be constructed using meryl.\n");
    fprintf(stderr, "    -min   m    Ignore kmers with value below m\n");
    fprintf(stderr, "    -max   m    Ignore kmers with value above m\n");
    fprintf(stderr, "    -threads t  Multithreading for meryl lookup table construction, dump and hist.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Memory usage can be limited, within reason, by sacrificing kmer lookup\n");
    fprintf(stderr, "  speed.  If the lookup table requires more memory than allowed, the program\n");
    fprintf(stderr, "  exits with an error.\n");
    fprintf(stderr, "    -memory1 m   Don't use more than m GB memory for loading seqmers\n");
    fprintf(stderr, "    -memory2 m   Don't use more than m GB memory for loading readmers\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -lookup file   Optional input vector of probabilities.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Exactly one report type must be specified.\n");
    fprintf(stderr, "\n\n");
    fprintf(stderr, "  -hist\n");
    fprintf(stderr, "   Generate a 0-centered k* histogram for sequences in <input.fasta>.\n");
    fprintf(stderr, "   Positive k* values are expected collapsed copies.\n");
    fprintf(stderr, "   Negative k* values are expected expanded  copies.\n");
    fprintf(stderr, "   Closer to 0 means the expected and found k-mers are well balenced, 1:1.\n");
    fprintf(stderr, "   Reports QV at the end, in stderr.\n");
    fprintf(stderr, "   Required: -sequence, -seqmers, -readmers, -peak, and -output.\n");
    fprintf(stderr, "   Optional: -lookup <probabilities> use probabilities to adjust multiplicity to copy number\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   Output: k* <tab> frequency\n");
    fprintf(stderr, "\n\n");
    fprintf(stderr, "  -dump\n");
    fprintf(stderr, "   Dump readK, asmK, and k* per bases (k-mers) in <input.fasta>.\n");
    fprintf(stderr, "   Required: -sequence, -seqmers, -readmers, -peak, and -output\n");
    fprintf(stderr, "   Optional: -skipMissing will skip the missing kmer sites to be printed\n");
    fprintf(stderr, "             -lookup <probabilities> use probabilities to adjust multiplicity to copy number\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   Output: seqName <tab> seqPos <tab> readK <tab> asmK <tab> k*\n");
    fprintf(stderr, "      seqName    - name of the sequence this kmer is from\n");
    fprintf(stderr, "      seqPos     - start position (0-based) of the kmer in the sequence\n");
    fprintf(stderr, "      readK      - normalized read copies (read multiplicity / peak)\n");
    fprintf(stderr, "      asmK       - assembly copies as found in <seq.meryl>\n");
    fprintf(stderr, "      k*         - 0-centered k* value\n");
    fprintf(stderr, "\n\n");
    fprintf(stderr, "  -vmer\n");
    fprintf(stderr, "   Score each variant, or variants within distance k and their combinations by k*.\n");
    fprintf(stderr, "   Required: -sequence, -seqmers, -readmers, -peak, -vcf, and -output\n");
    fprintf(stderr, "   Optional: -comb <N>  set the max N of combinations of variants to be evaluated (default: 15)\n"); 
    fprintf(stderr, "             -nosplit   without this options combinations larger than N are split\n");   
    fprintf(stderr, "             -disable-kstar  output variants by kstar. *experimental*\n");   
    //fprintf(stderr, "                        if chosen, use bcftools to compress and index, and consensus -H 1 -f <seq.fata> to polish.\n");
    //fprintf(stderr, "                        first ALT in heterozygous alleles are better supported by avg. |k*|.\n");
    fprintf(stderr, "             -lookup <probabilities> use probabilities to adjust multiplicity to copy number\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   Output files: <output>.debug and <output>.polish.vcf\n");
    fprintf(stderr, "    <output>.debug : some useful info for debugging.\n");
    fprintf(stderr, "                     seqName <tab> varMerStart <tab> varMerEnd <tab> varMerSeq <tab> score <tab> path\n");
    fprintf(stderr, "      varMerID                - unique numbering, starting from 0\n");
    fprintf(stderr, "      varMerRange             - seqName:start-end. position (0-based) of the variant (s), including sequences upstream and downstream of k-1 bp\n");
    fprintf(stderr, "      varMerSeq               - combination of variant sequence to evalute\n");
    fprintf(stderr, "      numMissings             - total number of missing kmers\n");
    fprintf(stderr, "      min k*                  - minimum of all |k*| for non-missing kmers. -1 when all kmers are missing.\n");
    fprintf(stderr, "      max k*                  - maximum of all |k*| for non-missing kmers. -1 when all kmers are missing.\n");
    fprintf(stderr, "      median k*               - median  of all |k*| for non-missing kmers. -1 when all kmers are missing.\n");
    fprintf(stderr, "      avg k*                  - average of all |k*| for non-missing kmers. -1 when all kmers are missing.\n");
    fprintf(stderr, "      avg ref-alt k*          - difference between reference and alternate average k*.\n");
    fprintf(stderr, "      delta kmer multiplicity - cumulative sum of kmer multiplicity variation. Positive values imply recovered kmers. Negative values overrepresented kmers introduced.\n");
    fprintf(stderr, "      record          - vcf record with <tab> replaced to <space>. only non-reference alleles are printed with GT being 1/1.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    <output>.polish.vcf : variants chosen.\n");
    //fprintf(stderr, "     use bcftools view -Oz <output>.polish.vcf and bcftools consensus -H 2 -f <seq.fata> to polish.\n");
    //fprintf(stderr, "     ALT alleles are favored with more support compared to the REF allele.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  omp_set_num_threads(threads);

  vector<string> copyKmerDict;

  if (!(pLookupTable == NULL)) {

    //  Read probabilities lookup table for 1-4 copy kmers.

    int it = 0;
    int r;
    double p;

    ifstream inputFile(pLookupTable, std::ios::in | std::ios::binary);

    if(inputFile.fail()){
      fprintf(stderr, "Error: failed to locate lookup table!\n");
      exit (EXIT_FAILURE);
    }

    // test file open   
    if (inputFile) {        
      string prob;

      // read the elements in the file into a vector  
      while ( inputFile >> prob ) {
        prob.erase(std::remove(prob.begin(), prob.end(), '\n'), prob.end());
        copyKmerDict.push_back(prob);
      }

      fprintf(stderr, "-- Loading copy-number lookup table '%s' (size '%lu').\n\n", pLookupTable, copyKmerDict.size());

      while (it < copyKmerDict.size()) {
        string s = copyKmerDict[it];
        std::string delimiter = ",";
        r = (int) stod(s.substr(0, s.find(delimiter)));
        p = (double) stod(s.erase(0, s.find(delimiter) + delimiter.length()));

        fprintf(stderr, "Copy-number: %u\t\tReadK: %u\tProbability: %f\n", it+1, r, p);

        it++;
      }

      fprintf(stderr, "\n");

      if (!inputFile.eof()) {
        cout << "Error: couldn't read lookup table!\n";
        exit (-1);
      }
    }
  }

  // pick readK function based on dictionary availability

  getreadK = &getreadKdef;

  if (!copyKmerDict.empty()) {
  	getreadK = &getreadKprob;	
  }	

  //  Open read kmers, build a lookup table.

  fprintf(stderr, "-- Loading kmers from '%s' into lookup table.\n", readDBname);

  merylFileReader*  readDB    = new merylFileReader(readDBname);
  merylExactLookup* readLookup = new merylExactLookup();

  readLookup->load(readDB, memory2, true, false, minV, maxV);

  // Open asm kmers, build a lookup table.
  
  fprintf(stderr, "\n-- Loading kmers from '%s' into lookup table.\n", seqDBname);
  merylFileReader*  asmDB      = new merylFileReader(seqDBname);
  merylExactLookup* asmLookup = new merylExactLookup();

  asmLookup->load(asmDB, memory1, true, false);

  delete readDB;   //  Not needed anymore.
  delete asmDB;    //  Not needed anymore.

  //  Open input sequences.
  dnaSeqFile  *seqFile = NULL;
  fprintf(stderr, "-- Opening sequences in '%s'.\n", seqName);
  seqFile = new dnaSeqFile(seqName);

  peak = ipeak;

  //  Check report type
  if (reportType == OP_HIST) {
    fprintf(stderr, "-- Generate histogram of the k* metric to '%s'.\n", outName);
    histKmetric(outName, seqName, seqFile, readLookup, asmLookup, copyKmerDict, threads);
  }
  if (reportType == OP_DUMP) {
    fprintf(stderr, "-- Dump per-base k* metric to '%s'.\n", outName);
    dumpKmetric(outName, seqName, seqFile, readLookup, asmLookup, skipMissing, copyKmerDict, threads);
  }
  if (reportType == OP_VAR_MER) {

    //  Open vcf file
    if (vcfName == NULL) {
      fprintf(stderr, "No variant call (-vcf) supplied.\n");
      exit (-1);
    }
    fprintf(stderr, "-- Opening vcf file '%s'.\n", vcfName);
    vcfFile* inVcf = new vcfFile(vcfName);

    fprintf(stderr, "-- Generate variant mers and score them.\n");
    varMers(seqName, seqFile, inVcf, readLookup, asmLookup, outName, comb, nosplit, copyKmerDict, bykstar, threads);

    delete inVcf;
  }
  
  delete seqFile;
  delete readLookup;
  delete asmLookup;

  fprintf(stderr, "Bye!\n");

  exit(0);
}
