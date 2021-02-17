
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



//  Compute a QV given the number of kmers in some set (kval) and the
//  total number of kmers in the universe (ktot).  
double
histoQV(double kval, double ktot) {
  double base   = kval / ktot;
  double kinv   = 1.0  / kmer::merSize();
  double qv     = -10.0 * log10( 1.0 - pow(1.0 - base, kinv));

  //fprintf(stderr, "histoQV(%f, %f) - 1.0-%f ^ %f - %f\n", kval, ktot, base, kinv, qv);

  return(qv);
}



void
processHistogram(void *G, void *T, void *S) {
  merfinGlobal  *g = (merfinGlobal *)G;
  //merfinHisto   *t = (merfinHisto *)T;
  merfinInput   *s = (merfinInput *)S;

  fprintf(stderr, "Processing sequence %s for histogram\n", s->seq.ident());

  //  Allocate space for our result.

  allocateArray(s->undr, s->undrMax, 1024);
  allocateArray(s->over, s->overMax, 1024);

  //  Iterate over every kmer in the sequence and ....

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
      continue;
    }

    if (asmK > readK) {
      uint32  idx = ((asmK / readK - 1) + 0.1) / 0.2;

      increaseArray(s->undr, idx, s->undrMax, 1024, _raAct::copyData | _raAct::clearNew);
      s->undr[idx]++;

      //  TODO: Check if this kmer was already counted. Only if not,
      //  overcpy += (asmK - readK)
      //  overcpy += (asmK - readK) / asmK
      s->koverCpy += (1.0 - readK / asmK) * prob;
    }

    else {
      uint32  idx = ((readK / asmK - 1) + 0.1) / 0.2;

      increaseArray(s->over, idx, s->overMax, 1024, _raAct::copyData | _raAct::clearNew);
      s->over[idx]++;
    }
  }
}



void
outputHistogram(void *G, void *S) {
  merfinGlobal  *g = (merfinGlobal *)G;
  merfinInput   *s = (merfinInput *)S;

  fprintf(stderr, "Output sequence %s\n", s->seq.ident());

  //  Allocate space if needed.

  if (g->histUndr == nullptr) {
    allocateArray(g->histUndr, g->histUndrMax, 2048);
    allocateArray(g->histOver, g->histOverMax, 2048);
  }

  //  Copy the per-sequence data into the global data

  g->histKmissing += s->kmissing;
  g->histKasm     += s->kasm;
  g->histKoverCpy += s->koverCpy;

  increaseArray(g->histUndr, s->undrMax, g->histUndrMax, 1024, _raAct::copyData | _raAct::clearNew);
  for (uint32 ii=0; ii<s->undrMax; ii++)
    g->histUndr[ii] += s->undr[ii];

  increaseArray(g->histOver, s->overMax, g->histOverMax, 1024, _raAct::copyData | _raAct::clearNew);
  for (uint32 ii=0; ii<s->overMax; ii++)
    g->histOver[ii] += s->over[ii];

  //  Emit a message for each read.

  fprintf(stderr, "%s\t%lu\t%lu\t%lu\t%.2f\n",
          s->seq.ident(),
          s->kmissing,
          g->histKmissing,
          s->kasm,
          histoQV(s->kmissing, s->kasm));

  delete s;
}



void
merfinGlobal::reportHistogram(void) {

  //  If no data to report, don't report.

  if ((histUndr == nullptr) ||
      (histOver == nullptr))
    return;

  //  Dump the histogram to a file.

  compressedFileWriter *k_hist = new compressedFileWriter(outName);

  for (uint64 ii = histUndrMax-1; ii > 0; ii--)
    if (histUndr[ii] > 0)
      fprintf(k_hist->file(), "%.1f\t%lu\n", ((double) ii * -0.2), histUndr[ii]);

  fprintf(k_hist->file(), "%.1f\t%lu\n", 0.0, histUndr[0] + histOver[0]);

  for (uint64 ii = 1; ii < histOverMax; ii++)
    if (histOver[ii] > 0)
      fprintf(k_hist->file(), "%.1f\t%lu\n", ((double) ii *  0.2), histOver[ii]);

  delete k_hist;

  //  Write a summary.

  fprintf(stderr, "\n");
  fprintf(stderr, "K-mers not found in reads (missing) : %lu\n", histKmissing);
  fprintf(stderr, "K-mers overly represented in assembly: %.2f\n", histKoverCpy);
  fprintf(stderr, "K-mers found in the assembly: %lu\n", histKasm);
  fprintf(stderr, "Missing QV: %.2f\n", histoQV(histKmissing,                histKasm));
  fprintf(stderr, "Merfin QV*: %.2f\n", histoQV(histKmissing + histKoverCpy, histKasm));
  fprintf(stderr, "*** Note this QV is valid only if -seqmer was generated with -sequence ***\n\n");
  fprintf(stderr, "*** Missing QV only considers missing kmers as errors. Merfin QV* includes overrepresented kmers. ***\n\n");
  fprintf(stderr, "*** When the lookup table is provided, missing QV includes weighted low frequency kmers, otherwise it is identical to Merqury QV. ***\n\n");
}
