/*
   Copyright (c) 2019, 2020 by Jianhua Yang <yangjh7@mail.sysu.edu.cn>
   endSeeker: A computational software for identifying 2’-O-Methylation sites from Nm-REP-seq data.
   Date: 2019/11/18 @ Sun Yat-sen University
*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include <getopt.h>
#include "BamReader.h"
#include "BamAux.h"
using namespace BamTools;
#include <map>
#include <algorithm>
#include <ios>
#include <iostream>
#include <string>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <vector>
using namespace std;

#include "bioUtils.h"
#include "faiFile.h"
#include "bedFile.h"
#include "samFile.h"

#include "endSeeker.h"

char version[] = "endSeeker version 0.1";
void usage(void);

int main(int argc, char *argv[])
{
  char *outFile   = NULL;
  char *ctrlFile  = NULL;
  char *treatFile = NULL;
  char *faFile    = NULL;
  char *faiFile   = NULL;
  char *geneFile  = NULL;
  char *typeStr   = NULL;

  FILE *genomefp  = NULL;
  FILE *faifp     = NULL;
  FILE *outfp     = NULL;
  FILE *bedfp     = NULL;

  int showVersion = 0;
  int showHelp    = 0;
  int c           = 0;
  int i           = 0;

  struct parameterInfo paraInfo;
  /* parse commmand line parameters */

  if (argc == 1)
  {
    usage();
  }

  const char *shortOptions = "vhVcUsno:p:m:f:l:L:t:F:t:a:A:T:I:w:r:g:M:y:";

  const struct option longOptions[] =
  {
    { "verbose" , no_argument , NULL, 'v' },
    { "help" , no_argument , NULL, 'h' },
    { "version" , no_argument , NULL, 'V' },
    { "collapser" , no_argument , NULL, 'c' },
    { "psi" , no_argument , NULL, 'U' },
    { "strand" , no_argument , NULL, 's' },
    { "norm" , no_argument , NULL, 'n' },
    { "outfile" , required_argument , NULL, 'o' },
    { "pvalue" , required_argument, NULL, 'p' },
    { "fold" , required_argument, NULL, 'f' },
    { "mfold" , required_argument, NULL, 'm' },
    { "min-len" , required_argument, NULL, 'l' },
    { "max-len" , required_argument, NULL, 'L' },
    { "min-tag" , required_argument, NULL, 't' },
    { "fdr" , required_argument , NULL, 'F' },
    { "rpm" , required_argument , NULL, 'r' },
    { "fa" , required_argument , NULL, 'a' },
    { "fai" , required_argument , NULL, 'A' },
    { "treat" , required_argument , NULL, 'T' },
    { "input" , required_argument , NULL, 'I' },
    { "gene" , required_argument , NULL, 'g' },
    { "window" , required_argument , NULL, 'w' },
    { "model" , required_argument , NULL, 'M' },
    { "type" , required_argument , NULL, 'y' },
    {NULL, 0, NULL, 0} ,  /* Required at end of array. */
  };

  paraInfo.verbose    = 0;
  paraInfo.genomeSize = 0;
  paraInfo.pvalue     = 1.0;
  paraInfo.fdr        = 5.0;
  paraInfo.mfold      = 0.1;
  paraInfo.fold       = 1.0;
  paraInfo.minTag     = 5;
  paraInfo.minLen     = 15;
  paraInfo.maxLen     = 1000;
  paraInfo.strand     = 1;
  paraInfo.psi        = 0;
  paraInfo.windowSize = 20;
  paraInfo.rpm        = 0.001;
  paraInfo.norm       = 0;
  paraInfo.collapser  = 0;
  paraInfo.geneModel  = 1; // Note for input genes
  paraInfo.type       = 1;

  while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0)
  {
    switch (c)
    {
    case 'v':
      paraInfo.verbose = 1;
      break;
    case 'h':
      showHelp = 1;
      break;
    case 'V':
      showVersion = 1;
      break;
    case 'c':
      paraInfo.collapser = 1;
      break;
    case 'U':
      paraInfo.psi = 1;
      break;
    case 's':
      paraInfo.strand = 1;
      break;
    case 'n':
      paraInfo.norm = 1;
      break;
    case 'o':
      outFile = optarg;
      break;
    case 'a':
      faFile  = optarg;
      break;
    case 'A':
      faiFile = optarg;
      break;
    case 'T':
      treatFile = optarg;
      break;
    case 'I':
      ctrlFile  = optarg;
      break;
    case 'g':
      geneFile  = optarg;
      break;
    case 'y':
      typeStr  = optarg;
      break;
    case 'p':
      paraInfo.pvalue = atof(optarg);
      break;
    case 'f':
      paraInfo.fold = atof(optarg);
      break;
    case 'm':
      paraInfo.mfold = atof(optarg);
      break;
    case 't':
      paraInfo.minTag = atof(optarg);
      break;
    case 'r':
      paraInfo.rpm = atof(optarg);
      break;
    case 'F':
      paraInfo.fdr = atof(optarg);
      break;
    case 'l':
      paraInfo.minLen = atoi(optarg);
      break;
    case 'L':
      paraInfo.maxLen = atoi(optarg);
      break;
    case 'w':
      paraInfo.windowSize = atoi(optarg);
      break;
    case 'M':
      paraInfo.geneModel = atoi(optarg);
      break;
    case '?':
      showHelp = 1;
      break;
    default:
      usage();
    }
  }

  // help for version
  if (showVersion)
  {
    fprintf(stderr, "%s", version);
    exit(1);
  }

  if (showHelp)
  {
    usage();
    exit(1);
  }

  if (treatFile == NULL)
  {
    fprintf(stderr, "ERROR: please set the option: --treat <treat alignments>\n");
    usage();
  }

  if (ctrlFile == NULL)
  {
    fprintf(stderr, "ERROR: please set the option: --input <rna-seq alignments>\n");
    usage();
  }

  if (faFile != NULL)
  {
    genomefp = (FILE *) fopen(faFile, "r");
    if (genomefp == NULL)
    {
      fprintf(stderr, "ERROR: Can't open genome file: %s\n", faFile);
      usage();
    }
  }
  else
  {
    fprintf(stderr, "ERROR: please set the option: --fa <genome file>\n");
    usage();
  }

  if (faiFile != NULL)
  {
    faifp = (FILE *) fopen(faiFile, "r");
    if (faifp == NULL)
    {
      fprintf(stderr, "ERROR: Can't open fai file: %s\nYou can use samtools faidx to generate fai file\n", faiFile);
      usage();
    }
  }
  else
  {
    fprintf(stderr, "ERROR: please set the option: --fai <fai file>\n");
    usage();
  }

  if (outFile == NULL)
  {
    outfp = stdout;
  }
  else
  {
    outfp = (FILE *) fopen(outFile, "w");
    if (outfp == NULL)
    {
      fprintf(stderr, "ERROR: Can't open %s\n", outFile);
      usage();
    }
  }

  if (geneFile != NULL)
  {
    bedfp = (FILE *) fopen(geneFile, "r");
    if (bedfp == NULL)
    {
      fprintf(stderr, "ERROR: Can't open gene file: %s\n", geneFile);
      usage();
    }
    if (paraInfo.geneModel != 2)
      paraInfo.geneModel = 1;
  }
  else
  {
    fprintf(stderr, "ERROR: please set the option: --gene <gene file, BED12 format>\n");
  }

  if (typeStr != NULL)
  {
    if (strcmp(typeStr, "mgr") == 0) paraInfo.type = 2;
    if (strcmp(typeStr, "rep") == 0) paraInfo.type = 1;
  }
  if (paraInfo.windowSize < 10)
  {
    fprintf(stderr, "the window size must be large than 10 with -w option\n");
    usage();
  }

  runEndSeeker(&paraInfo, genomefp, faifp, bedfp, treatFile, ctrlFile, outfp);

  fclose(genomefp);
  fclose(faifp);
  fclose(outfp);
  if (bedfp != NULL) fclose(bedfp);
  return 0;
}

void usage(void)
{
  fprintf(stderr, "%s", "Usage:  endSeeker [options] --fa <genome seq> --fai <fai file> --gene <bed12 file> --treat <alignments> --input <input alignments>\n\
endSeeker: A computational software for identifying 2’-O-Methylation sites from Nm-REP-seq data.\n\
[options]\n\
--fa <string>                  : genome sequence<fasta format>[Required]\n\
--fai <string>                 : genome fai file<fai format>[Required].\n\
                                 using \"samtools faidx\" to generate fai file\n\
--gene <string>                : gene file <BED12 format>[Required]\n\
--treat <string>               : file treated by mgR/mgR+OED file<BAM format>[Required]\n\
--input <string>               : input file<BAM format>[Required]\n\
-v/--verbose                   : verbose information\n\
-V/--version                   : endSeeker version\n\
-h/--help                      : help informations\n\
-s/--strand                    : strand-specific sequencing data\n\
-n/--norm                      : normalized reads to the locus number\n\
-c/--collapser                 : keep duplication, deault is false\n\
-o/--outfile <string>          : output file\n\
-t/--min-tag <double>          : minimum tag number for each end site, default>=5.0 read\n\
-r/--rpm <double>              : minimum rpm value for each end site, default>=0.001\n\
-f/--fold <int>                : minimum fold-change[default>=1.0]\n\
-w/--window <int>              : window size around the end position[default=20]\n\
-l/--min-len <int>             : minimum length of reads, default=15\n\
-L/--max-len <int>             : maximum length of reads, default=1000\n\
");
  exit(1);
}
