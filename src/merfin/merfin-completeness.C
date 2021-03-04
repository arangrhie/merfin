
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
#include "kmers.H"


class dbInput {
public:
  dbInput(char *path, uint32 piece) {
    _stream = new merylFileReader(path);
    _stream->enableThreads(piece);

    nextMer();
  }

  ~dbInput() {
    delete _stream;
  }

  bool      nextMer(void) {
    _valid = _stream->nextMer();
    _kmer  = _stream->theFMer();
    _value = _stream->theValue();

    return(_valid);
  }

  merylFileReader *_stream;
  bool             _valid;
  kmer             _kmer;
  kmvalu           _value;
};


void
computeCompleteness(merfinGlobal *G) {
  dbInput  *readDB[64];
  dbInput  *asmDB[64];

  //  Open all 64 meryl pieces and clear intermediate storage (that should
  //  already be cleared).

  for (uint32 ii=0; ii<64; ii++) {
    readDB[ii] = new dbInput(G->readDBname, ii);
    asmDB[ii]  = new dbInput(G->seqDBname, ii);

    G->complTotal[ii] = 0.0;
    G->complUndrC[ii] = 0.0;
  }


#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 ii=0; ii<64; ii++) {
    dbInput *S = readDB[ii];
    dbInput *A = asmDB[ii];

    while ((S->_valid == true) ||
           (A->_valid == true)) {
      double  readK = 0.0;
      double  asmK  = 0.0;
      double  prob  = 0.0;


      //  If the kmers are equal, compute readK, asmK and prob as usual.
      //
      if        ((S->_valid == true) &&
                 (A->_valid == true) &&
                 (S->_kmer == A->_kmer)) {
        G->getK(S->_value, A->_value, readK, asmK, prob);

        S->nextMer();
        A->nextMer();
      }

      //  If there is no assmebly kmer for this read kmer,
      //   - either A is not a valid kmer or
      //   -        A is larger than S
      //  compute readK, asmK and prob with 0 for the assembly value.
      //
      else if ((S->_valid == true) && ((A->_valid == false) ||
                                       (S->_kmer < A->_kmer))) {
        G->getK(S->_value, 0, readK, asmK, prob);

        S->nextMer();
      }

      //  Otherwise there is no read kmer for this assembly kmer
      //   - either S is not a valid kmer or
      //   -        S is larger than A
      //  do nothing, in particular, skip the summations after this
      //  if-else-else block.
      //
      else {
        A->nextMer();
        continue;
      }

      //  Now with valid readK, asmK and prob, sum our metric.

      G->complTotal[ii] += readK;

      if (readK > asmK)
        G->complUndrC[ii] += readK - asmK;
    }

    fprintf(stderr, "thread %2u total %12.2f underc %15.5f completeness %0.8f\n",
            ii, G->complTotal[ii], G->complUndrC[ii], 1.0 - G->complUndrC[ii] / G->complTotal[ii]);
  }

  //  Cleanup inputs.

  for (uint32 ii=0; ii<64; ii++) {
    delete  readDB[ii];
    delete  asmDB[ii];
  }

  //  Report.

  double  total = 0;
  double  undrc = 0;

  for (uint32 ii=0; ii<64; ii++) {
    total += G->complTotal[ii];
    undrc += G->complUndrC[ii];
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "TOTAL readK:   %15.2f\n", total);
  fprintf(stderr, "TOTAL undrcpy:    %15.5f\n", undrc);
  fprintf(stderr, "COMPLETENESS:             %0.5f\n", 1.0 - undrc / total);
}



void
merfinGlobal::reportCompleteness(void) {
}
