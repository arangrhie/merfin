
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



void
processDump(void *G, void *T, void *S) {
  merfinGlobal  *g = (merfinGlobal *)G;
  //merfinHisto   *t = (merfinHisto *)T;
  merfinInput   *s = (merfinInput *)S;

  fprintf(stderr, "Processing sequence %s for dumping\n", s->seq.ident());

  //  Output for each base in the input sequence, if 'skipMissing' == false:
  //    %s\t%lu\t%.2f\t%.2f\t%.2f\n  - seq.ident(), kiter.position(), readK, asmK, kMetric
  //
  //  Output for each sequence:
  //    "%s\t%lu\t%lu\t%lu\n", seq.ident(), seq_missing, total_missing, total_kmers

  if (g->skipMissing == false) {
    allocateArray(s->dumpReadK,   s->seq.length() + 1, _raAct::clearNew);
    allocateArray(s->dumpAsmK,    s->seq.length() + 1, _raAct::clearNew);
    allocateArray(s->dumpKMetric, s->seq.length() + 1, _raAct::clearNew);
  }

  double  readK   = 0;
  double  asmK    = 0;
  double  prob    = 0;

  while (s->kiter.nextBase()) {
    if (s->kiter.isValid() == false)
      continue;

    s->kasm++;

    //  This is mirroring merfinGlobal::getKmetric(), but directly computing
    //  a histogram index instead of reversing the compute of Kmetric.

    g->getK(s->kiter.fmer(), s->kiter.rmer(),
            readK, asmK, prob);

    if (readK == 0) {
      s->kmissing++;
    }

    if (g->skipMissing == false) {
      uint64 pp = s->kiter.position();

      s->dumpReadK[pp]   = readK;
      s->dumpAsmK[pp]    = asmK;
      s->dumpKMetric[pp] = g->getKmetric(readK, asmK);
    }
  }
}



void
outputDump(void *G, void *S) {
  merfinGlobal  *g = (merfinGlobal *)G;
  merfinInput   *s = (merfinInput *)S;

  fprintf(stderr, "Output sequence %s\n", s->seq.ident());

  //  Open the output if we're writing output and don't have it open already.

  if ((g->skipMissing == false) &&
      (g->dumpOutFile == nullptr))
    g->dumpOutFile = new compressedFileWriter(g->outName);

  //  Write the output if we're writing output.

  if (g->skipMissing == false)
    for (uint64 pp=0; pp<s->seq.length(); pp++)
      if ((s->dumpReadK[pp]   != 0.0) ||
          (s->dumpAsmK[pp]    != 0.0) ||
          (s->dumpKMetric[pp] != 0.0))
        fprintf(g->dumpOutFile->file(), "%s\t%lu\t%.2f\t%.2f\t%.2f\n",
                s->seq.ident(), pp, s->dumpReadK[pp], s->dumpAsmK[pp], s->dumpKMetric[pp]);

  //  Add to the global missing/asm counts, write a message for each sequence.

  g->histKmissing += s->kmissing;
  g->histKasm     += s->kasm;

  fprintf(stderr, "%s\t%lu\t%lu\t%lu\n",
          s->seq.ident(), s->kmissing, g->histKmissing, g->histKasm);

  delete s;
}
