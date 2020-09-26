
/******************************************************************************
 *
 *  This file is part of meryl-utility, a collection of miscellaneous code
 *  used by Meryl, Canu and others.
 *
 *  This software is based on:
 *    'Canu' v2.0              (https://github.com/marbl/canu)
 *  which is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */
 
#ifndef KMETRIC_H
#define KMETRIC_H

#include "runtime.H"
#include "types.H"
#include "kmers.H"

extern uint64 peak;

typedef double (*func_t)(
  		   merylExactLookup   *rlookup,
           merylExactLookup   *alookup,
           kmer                fmer,
           kmer                rmer,
           vector<string>	  &copyKmerDict,
           double             &readK,
           double             &asmK,
           double			  &prob);
           
double
getKmetricDef(
  		   merylExactLookup   *rlookup,
           merylExactLookup   *alookup,
           kmer                fmer,
           kmer                rmer,
           vector<string>     &copyKmerDict,
           double             &readK,
           double             &asmK,
           double			  &prob);
           
double
getKmetricProb(
  		   merylExactLookup   *rlookup,
           merylExactLookup   *alookup,
           kmer                fmer,
           kmer                rmer,
           vector<string>     &copyKmerDict,
           double             &readK,
           double             &asmK,
           double			  &prob);
           
#endif  //  KMETRIC_H
