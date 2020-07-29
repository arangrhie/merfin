
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
#include "types.H"

#include <vector>
#include <map>

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
traverse(int     idx,
  vector<uint32>    refIdxList,
  vector<uint32>    refLenList,
  map<int, vector<char*> >  posHaps,
  string         candidate,
  vector<int>    path,
  varMer*        seqMer) {

  if (idx < refIdxList.size()) {

    vector<char*> haps = posHaps.at(idx);
    for (int j = 0; j < haps.size(); j++) {

      path.push_back(j);

      char* hap = haps.at(j);

      //  Make the replaced ver.
      string replaced = replace(candidate, refIdxList.at(idx), refLenList.at(idx), hap);

      //  Shift refIdx
      int delta = strlen(hap) - refLenList.at(idx);

      //  Apply to the rest of the positions
      for (int i = idx + 1; i < refIdxList.size(); i++) {
        refIdxList.at(i) += delta;
      }

      //  Change refLen
      refLenList.at(idx) = strlen(hap);

      //  Traverse
      replaced = traverse(idx + 1, refIdxList, refLenList, posHaps, replaced, path, seqMer);

      //  Only do something when all nodes are visited
      if (idx == refIdxList.size() - 1) {
        seqMer->addSeqPath(replaced, path);
      }

      path.pop_back();

    }

  }
  idx--;
  return candidate;
}

void
histKmetric(char               *outName,
            dnaSeqFile         *sfile,
            merylExactLookup   *rlookup,
            merylExactLookup   *alookup,
            uint64              peak) {

  dnaSeq seq;
  uint64 fValue = 0;
  uint64 rValue = 0;
  double readK;
  double asmK;
  double kMetric;



  //  compressedFileWriter *k_values = new compressedFileWriter(concat(outName, ".gz"));
  compressedFileWriter *k_hist   = new compressedFileWriter(concat(outName, ".hist"));


  //  variables for generating histogram
  uint64   histMax  = 32 * 1024 * 1024;
  uint64 * overHist = new uint64[histMax];	// positive k* values, overHist[0] = bin 0.0 ~ 0.2
  uint64 * undrHist = new uint64[histMax];	// negative k* values
  uint64   missing  = 0;			// missing kmers (0)
  uint64   histBin  = 0;

  for (uint64 ii = 0; ii < histMax; ii++) {
    overHist[ii] = 0;
    undrHist[ii] = 0;
  }

  for (uint32 seqId=0; sfile->loadSequence(seq); seqId++) {
    kmerIterator kiter(seq.bases(), seq.length());

    while (kiter.nextBase()) {
      if (kiter.isValid() == true) {
        fValue = 0;
        rValue = 0;
        rlookup->exists(kiter.fmer(), fValue);
        rlookup->exists(kiter.rmer(), rValue);

        readK = (double) (fValue + rValue) / peak;

        fValue = 0;
        rValue = 0;
        alookup->exists(kiter.fmer(), fValue);
        alookup->exists(kiter.rmer(), rValue);

        asmK  = (double) (fValue + rValue);

        if ( readK == 0 ) {
          kMetric = 0;
          missing++;
        } else if ( asmK > readK ) {
          kMetric = asmK / readK - 1;
          undrHist[(uint64) ((kMetric + 0.1) / 0.2)]++;
          kMetric = kMetric * -1;
        } else { // readK > asmK
          kMetric = readK / asmK - 1;
          overHist[(uint64) ((kMetric + 0.1 ) / 0.2)]++;
        }

/*** Need to put this out as a separate function
        fprintf(k_values->file(), "%s\t%lu\t%.2f\t%.2f\t%.2f\n",
                 seq.name(),
                 kiter.position(),
                 readK,
                 asmK,
                 kMetric
                 );
***/
      }
    }
  }

  for (uint64 ii = histMax - 1; ii > 0; ii--) {
     if (undrHist[ii] > 0)  fprintf(k_hist->file(), "%.1f\t%lu\n", ((double) ii * -0.2), undrHist[ii]);
  }
  fprintf(k_hist->file(), "%.1f\t%lu\n", 0.0, (undrHist[0]+overHist[0]));
  for (uint64 ii = 1; ii < histMax; ii++) {
     if (overHist[ii] > 0)  fprintf(k_hist->file(), "%.1f\t%lu\n", ((double) ii * 0.2), overHist[ii]);
  }
}


void
generateVarMers(dnaSeqFile       *sfile,
                vcfFile          *vfile,
                merylExactLookup *rlookup,
                merylExactLookup *alookup,
                uint64            peak,
                char             *out) {

  //  output file
  compressedFileWriter *oFasta  = new compressedFileWriter(concat(out, ".fasta"));
  compressedFileWriter *oKAfter = new compressedFileWriter(concat(out, ".k_after"));
  compressedFileWriter *oDebug  = new compressedFileWriter(concat(out, ".debug"));

  //  What is the kmer size?
  uint32 ksize = kmer::merSize();

  //  Merge posGTlist for each chr within ksize
  fprintf(stderr, "Merge variants within %u-mer bases\n", ksize);
  vfile->mergeChrPosGT(ksize);

  // print CHR rStart rEnd POS HAP1 HAP2 minHAP1 minHAP2 to out.debug
  
  map<string, vector<posGT*>*> *mapChrPosGT = vfile->_mapChrPosGT;

  vector<posGT*>  *posGTlist;
  vector<gtAllele*>     *gts;

  string   seqName;
  char     fString[65];
  char     rString[65];
  string   kmer;
  dnaSeq   seq;
  uint64   fValue = 0;
  uint64   rValue = 0;
  
  posGT    *posGt;
  gtAllele *gt;
  uint32   pos;
  uint32   rStart;
  uint32   rEnd;
  uint32   K_PADD = ksize - 1;

  uint32   hapIdx = 0;

  vector<uint32> refIdxList;
  vector<uint32> refLenList;
  map<int, vector<char*> > mapPosHap;
  vector<char*>   haps;
  vector<int>     path;
 
  // temporal sequence to hold ref bases
  char      refTemplate[1024] = "";

  for (uint32 seqId=0; sfile->loadSequence(seq); seqId++) {
    //  for each seqId
    seqName = string(seq.name());
    fprintf(stderr, "\nProcessing \'%s\'\n", seq.name());

    //  in case no seq.name() available, ignore this seqName
    if (mapChrPosGT->find(seqName) == mapChrPosGT->end())
      continue;

    //  get chr specific posGTs
    posGTlist   = mapChrPosGT->at(seq.name());

    //  get sequence combinations on each posGT list
    for (int posGtIdx = 0; posGtIdx < posGTlist->size(); posGtIdx++) {

      // initialize variables
      posGt       = posGTlist->at(posGtIdx);
      rStart      = posGt->_rStart;  // 0-based
      if (rStart > K_PADD) { rStart -= K_PADD; }
      else { rStart = 0; }

      rEnd        = posGt->_rEnd;             // 1-based
      if (rEnd < seq.length() - K_PADD) {  rEnd += K_PADD;  }
      else { rEnd = seq.length(); }

      gts         = posGt->_gts;
      refIdxList.clear();
      refLenList.clear();
      path.clear();
      haps.clear();
      mapPosHap.clear();

      //  load original sequence from rStart to rEnd
      if ( ! seq.copy(refTemplate, rStart, rEnd, true )) {
        fprintf(stderr, "Invalid region specified: %s : %u - %u\n", seq.name(), rStart, rEnd);
        return;
      }
      // DEBUG  fprintf(stderr, "%s\n", refTemplate);
       
      //  load mapPosHap
      for (int i = 0; i < gts->size(); i++) {
         gt = gts->at(i);
         refIdxList.push_back(gt->_pos - rStart);
         refLenList.push_back(gt->_refLen);
         haps.push_back(gt->_hap1);
         haps.push_back(gt->_hap2);

         //  include ref allele in case it was not in the hap1 hap2
         if (!gt->_hasRef) {
           haps.push_back(gt->_ref);
         }
         mapPosHap.insert(pair<int, vector<char*> >(i, haps));
         fprintf(stderr, "%s %u %u %u %s %s\n",
                 seq.name(), rStart, rEnd, gt->_pos, gt->_hap1, gt->_hap2);
         haps.clear();
      }

      varMer* seqMer = new varMer(posGt);

      //  traverse through each gt combination
      traverse(0, refIdxList, refLenList, mapPosHap, refTemplate, path, seqMer);

      //  score each combination
      seqMer->score(rlookup, alookup, peak);

      //  print output sorted haps
      fprintf(oDebug->file(), "Sorted by minFreq\n");
      multimap<uint64, int> minFreqs = seqMer->getMinFreqs();
      for ( pair<uint64, int> p : minFreqs ) {
        if ( p.first > 0 ) {
          int idx = p.second;

          fprintf(oDebug->file(), "%s\t%u\t%u\t%s\t%lu\tpath =",
          seq.name(),
          rStart,
          rEnd,
          seqMer->seqs.at(idx).c_str(),
          p.first);

          for ( int i = 0; i < seqMer->paths.at(idx).size(); i++) {
            fprintf(oDebug->file(), "\t%u ( %d )", gts->at(i)->_pos, seqMer->paths.at(idx).at(i));
          }

          fprintf(oDebug->file(), "\n");
        }
      }

      fprintf(oDebug->file(), "Sorted by k*\n");
      multimap<double, int> minKs = seqMer->getMinKs();
      for ( pair<double, int> p : minKs ) {
        if ( p.first > 0.0 ) {
          int idx = p.second;

          fprintf(oDebug->file(), "%s\t%u\t%u\t%s\t%.1f\tpath =",
          seq.name(),
          rStart,
          rEnd,
          seqMer->seqs.at(idx).c_str(),
          p.first);

          for ( int i = 0; i < seqMer->paths.at(idx).size(); i++) {
            fprintf(oDebug->file(), "\t%u ( %d )", gts->at(i)->_pos, seqMer->paths.at(idx).at(i));
          }

          fprintf(oDebug->file(), "\n");
        }
      }

      fprintf(oDebug->file(), "\n");
      //  polish seq
    }
  }
}


int
main(int argc, char **argv) {
  char           *seqName    = NULL;
  char           *vcfName    = NULL;
  char           *outName    = NULL;
  char           *seqDBname  = NULL;
  char           *readDBname  = NULL;

  uint64          minV       = 0;
  uint64          maxV       = UINT64_MAX;
  uint64          peak       = 0;

  uint32          threads    = omp_get_max_threads();
  uint32          memory     = 0;

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-sequence") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-seqmers") == 0) {
      seqDBname = argv[++arg];

    } else if (strcmp(argv[arg], "-readmers") == 0) {
      readDBname = argv[++arg];

    } else if (strcmp(argv[arg], "-peak") == 0) {
      peak = strtouint64(argv[++arg]);

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

    } else if (strcmp(argv[arg], "-memory") == 0) {
      memory = strtouint32(argv[++arg]);

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
  if (vcfName == NULL)
    err.push_back("No variant call (-vcf) supplied.\n");
  if (peak == 0)
    err.push_back("Peak=0 or no haploid peak (-peak) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s \\\n", argv[0]);
    fprintf(stderr, "         -sequence <input.fasta>   \\\n");
    fprintf(stderr, "         -seqmers  <input.meryl>   \\\n");
    fprintf(stderr, "         -readmers <read.meryl>    \\\n");
    fprintf(stderr, "         -peak     <haploid_peak>  \\\n");
    fprintf(stderr, "         -vcf      <input.vcf>     \\\n");
    fprintf(stderr, "         -output   <output>        \n\n");
    fprintf(stderr, "  Predict the kmer from <input.vcf> given sequence <input.fasta>\n");
    fprintf(stderr, "  and lookup the k-mer multiplicity from <input.meryl>.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Input files can be FASTA or FASTQ; uncompressed, gz, bz2 or xz compressed\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Output:\n");
    fprintf(stderr, "    <output>.fasta    - polished sequence\n");
    fprintf(stderr, "    <output>.alt.vcf  - alternative variants (not chosen due to lower multiplicity support\n");
    fprintf(stderr, "    <output>.k_before - k* hist before, 0-centered\n");
    fprintf(stderr, "    <output>.k_after  - k* hist after,  0-centered\n");
    fprintf(stderr, "    <output>.debug    - kmers generated from each variant, multiplicity, K* before and after (est)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Each input database can be filtered by value.  More advanced filtering\n");
    fprintf(stderr, "  requires a new database to be constructed using meryl.\n");
    fprintf(stderr, "    -min   m    Ignore kmers with value below m\n");
    fprintf(stderr, "    -max   m    Ignore kmers with value above m\n");
    fprintf(stderr, "    -threads t  Number of threads to use when constructing lookup table.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Memory usage can be limited, within reason, by sacrificing kmer lookup\n");
    fprintf(stderr, "  speed.  If the lookup table requires more memory than allowed, the program\n");
    fprintf(stderr, "  exits with an error.\n");
    fprintf(stderr, "    -memory m   Don't use more than m GB memory\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  omp_set_num_threads(threads);

  //  Open read kmers, build a lookup table.

  fprintf(stderr, "-- Loading kmers from '%s' into lookup table.\n", readDBname);

  merylFileReader*  merylDB    = new merylFileReader(readDBname);
  merylExactLookup* readLookup = new merylExactLookup(merylDB, memory, minV, maxV);

  if (readLookup->configure() == false)
    exit(1);

  readLookup->load();

  // Open asm kmers, build a lookup table.
  
  fprintf(stderr, "-- Loading kmers from '%s' into lookup table.\n", seqDBname);
  merylDB = new merylFileReader(seqDBname);
  merylExactLookup* asmLookup = new merylExactLookup(merylDB, memory, 0, UINT64_MAX);

  if (asmLookup->configure() == false)
    exit(1);

  asmLookup->load();

  delete merylDB;   //  Not needed anymore.



  //  Open input sequences.

  dnaSeqFile  *seqFile = NULL;
  if (seqName != NULL) {
    fprintf(stderr, "-- Opening sequences in '%s'.\n", seqName);
    seqFile = new dnaSeqFile(seqName);
  }

  //  Open vcf file
  fprintf(stderr, "-- Opening vcf file '%s'.\n", vcfName);
  vcfFile* inVcf = new vcfFile(vcfName);

  //  Dump before kmetric

/*** Comment when testing takes too long ***
  fprintf(stderr, "-- Dump k* metrics to '%s.k_before.gz' and '%s.k_before.hist'.\n", outName, outName);
  histKmetric(concat(outName, ".k_before"), seqFile, readLookup, asmLookup, peak);
  fprintf(stderr, "\n");
*******************************************/

  //  seqFile->reopen();
  delete seqFile;
  seqFile = new dnaSeqFile(seqName);

  fprintf(stderr, "-- Generate variant mers and score them.\n");
  generateVarMers(seqFile, inVcf, readLookup, asmLookup, peak, outName);

  //  Dump after kmetric
  
  delete seqFile;
  delete inVcf;
  delete readLookup;
  delete asmLookup;

  fprintf(stderr, "Bye!\n");

  exit(0);
}
