
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
#include "system.H"
#include "sequence.H"
#include "kmers.H"
#include "files.H"

#include "vcf.H"

#include <string>
#include <vector>

#ifndef MERFIN_GLOBAL
#define MERFIN_GLOBAL

#define OP_NONE       0
#define OP_HIST       1
#define OP_DUMP       2
#define OP_VAR_MER    3
#define OP_COMPL      4



////////////////////////////////////////
//
//  This object holds the input sequence and output results.  A new object is
//  created for each input, passed to the compute function, then passed to
//  the output function and deleted there.
//
class merfinInput {
public:
public:
  merfinInput() {
  };
  ~merfinInput() {
    delete [] undr;
    delete [] over;

    delete [] dumpReadK;
    delete [] dumpAsmK;
    delete [] dumpKMetric;

    //delete [] vmers;
  };

  //  Input data.

  dnaSeq        seq;
  kmerIterator  kiter;

  //  Histogram results.

  uint64     kasm        = 0;   //  number of kmers in the sequence processed
  uint64     kmissing    = 0;   //  number of kmers missing from the db
  double     koverCpy    = 0;   //  double, because it's summing fractions

  uint32     undrMax     = 0;
  uint64    *undr        = nullptr;

  uint32     overMax     = 0;
  uint64    *over        = nullptr;

  //  Dump results (also uses kasm and kmissing from above).

  double    *dumpReadK   = nullptr;
  double    *dumpAsmK    = nullptr;
  double    *dumpKMetric = nullptr;

  //  Variant results.

  //uint64     vmersLen    = 0;
  //uint64     vmersMax    = 0;
  //char      *vmers       = nullptr;
};



////////////////////////////////////////
//
//  This object holds per-thread data for -vmer work.
//   - a handful of arrays that we can recycle for each input
//   - a per-thread debug logging file
//
class
merfinThrData {
public:
  merfinThrData() {
  }
  ~merfinThrData() {
    delete oDebug;
    delete oVcf;
  }

  //  Our ident.

  uint32                  threadID;

  //  For debug logging.

  compressedFileWriter   *oDebug   = nullptr;
  uint64                  varMerId = 0;

  compressedFileWriter   *oVcf     = nullptr;
};



////////////////////////////////////////
//
//  This object holds command line parameters and input data.  The input data
//  should be treated as read-only; multiple threads will access it with no
//  locking.
//
//  This object also holds any summary data for output.  This data should be
//  modified only in the sweatshop 'output' function - no locking is needed as
//  the sweatshop guarantees only one thread is active at a time here.
//
class merfinGlobal {
public:
  ~merfinGlobal() {
    delete seqFile;
    delete inVcf;
    delete asmLookup;
    delete readLookup;
  }

public:
  void      load_Kmetric(void);
  void      load_Kmers(void);
  void      open_Inputs(void);

  void      getK(kmvalu   seqValue,
                 kmvalu   asmValue,
                 double  &readK, double  &asmK, double  &prob);
  void      getK(kmer     fmer,
                 kmer     rmer,
                 double  &readK, double  &asmK, double  &prob);
  double    getKmetric(double readK, double asmK);

  //  Parameters for histogram mode.
public:

  //  Parameters for dump mode.
  //   - dumpOutFile is created on demand in the output thread.
public:
  bool                    skipMissing  = false;
  compressedFileWriter   *dumpOutFile  = nullptr;

  //  Parameters for vmer mode.
  //   - dumpOutFile is created on demand in the output thread.
public:
  compressedFileWriter   *oVCF         = nullptr;





  //  Output data for histogram and data modes.  Data mode only uses histKasm
  //  and histKmissing.
public:
  uint64                  histKasm     = 0;   //  number of kmers in all seqs processed
  uint64                  histKmissing = 0;   //  number of kmers missing from the db
  double                  histKoverCpy = 0;   //  double, because it's summing fractions

  uint32                  histUndrMax  = 0;
  uint64                 *histUndr     = nullptr;

  uint32                  histOverMax  = 0;
  uint64                 *histOver     = nullptr;

  void                    reportHistogram(void);
  void                    reportCompleteness(void);

  //  Output data for variant mode.
public:

  //  Output data for completeness
public:
  double                  complTotal[64] = {0.0};
  double                  complUndrC[64] = {0.0};


  //  Input data.
public:
  char                   *seqName       = nullptr;
  dnaSeqFile             *seqFile       = nullptr;

  char                   *vcfName       = nullptr;
  vcfFile                *inVcf         = nullptr;

  char                   *outName       = nullptr;

  char                   *seqDBname     = nullptr;
  merylExactLookup       *asmLookup;

  char                   *readDBname    = nullptr;
  merylExactLookup       *readLookup;

  uint64                  minV          = 0;           //  Kmer value limits for the
  uint64                  maxV          = uint64max;   //  above two lookup tables.

  char                   *pLookupTable  = nullptr;
  std::vector<uint32>     copyKmerK;
  std::vector<double>     copyKmerP;

  //  Run-time parameters.
public:
  uint32                  threads       = omp_get_max_threads();
  double                  maxMemory     = getMaxMemoryAllowed() / 1024.0 / 1024.0 / 1024.0;

  uint32                  reportType    = OP_NONE;





public:
  double                  peak          = 0;
  bool                    nosplit       = false;
  bool                    bykstar       = true;
  uint32                  comb          = 15;
};



inline
double
merfinGlobal::getKmetric(double readK, double asmK) {
  if (readK == 0)
    return(0);

  if (asmK > readK)
    return((asmK / readK - 1) * -1);

  if (asmK < readK)
    return(readK / asmK - 1);

  return(0);  //  Otherwise equal, return zero.
}


#endif  //  MERFIN_GLOBAL
