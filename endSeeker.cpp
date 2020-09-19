/*
   Copyright (c) 2019, 2020 by Jianhua Yang <yangjh7@mail.sysu.edu.cn>
   endSeeker: A computational software for identifying 2’-O-Methylation sites from Nm-REP-seq data.
   Date: 2019/11/18 @ Sun Yat-sen University
*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<assert.h>
#include<math.h>
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

rtsSiteVector rtsSiteList;

int cmpScorePval(const rtsSiteInfo *x, const rtsSiteInfo *y)
{
  if (x->score != y->score)
    return x->score < y->score;
  else
    return x->siteModel < y->siteModel;
}

template<typename T>
string NumberToString(T Number)
{
  ostringstream ss;
  ss << Number;
  return ss.str();
}

void runEndSeeker(struct parameterInfo *paraInfo, FILE *genomefp,
                  FILE *faifp, FILE *bedfp,
                  char *treatBamFile, char *ctrlBamFile, FILE *outfp)
{
  bool skipDeletion = 0;
  bool skipSplice = 0;
  double treatTotalNum = 0;
  double ctrlTotalNum  = 0;
  long   genomeSize    = 1;

  BamReader treatReader;
  BamReader ctrlReader;

  chromBed12Map treatBedHash;
  chromBed12Map ctrlBedHash;

  chromStrandBed12Map bed12GeneHash;

  if (!treatReader.Open(treatBamFile))
  {
    cerr << "Failed to open BAM file " << treatBamFile << endl;
    exit(1);
  }
  if (!ctrlReader.Open(ctrlBamFile))
  {
    cerr << "Failed to open BAM file " << ctrlBamFile << endl;
    exit(1);
  }

  fprintf(stderr, "#read genome fai file\n");
  faidxMap faiHash;
  paraInfo->genomeSize = readFai(faifp, faiHash);

  fprintf(stderr, "#get chromosome size\n");
  chromSizeMap mapSize;
  paraInfo->genomeSize = getGenomeSize(treatReader, mapSize);

  if (paraInfo->geneModel)
  {
    fprintf(stderr, "#get transcript annotation information\n");
    int geneNum = getBedToMap(bedfp, bed12GeneHash);
    fprintf(stderr, " get %d transcript from annotation file\n", geneNum);
  }

  fprintf(stderr, "#read treatment bam file to bed list\n");
  treatTotalNum = readSEBamToBed12Map(treatReader, treatBedHash, skipDeletion, skipSplice, paraInfo->collapser);
  fprintf(stderr, "#read control bam file to bed list\n");
  ctrlTotalNum = readSEBamToBed12Map(ctrlReader, ctrlBedHash, skipDeletion, skipSplice, paraInfo->collapser);
  fprintf(stderr, " treatment readNum=%.f control readNum=%.f\n", treatTotalNum, ctrlTotalNum);

  if (paraInfo->norm)
  {
    fprintf(stderr, "#normalizing reads in treatment sample\n");
    treatTotalNum = normalizedBed12Reads(treatBedHash);
    fprintf(stderr, "#normalizing reads in control sample\n");
    ctrlTotalNum = normalizedBed12Reads(ctrlBedHash);
  }

  paraInfo->tTotalNum = treatTotalNum;
  paraInfo->cTotalNum = ctrlTotalNum;

  outputHeader(paraInfo, outfp);

  callEndSites(paraInfo, mapSize, bed12GeneHash, genomefp, faiHash, treatBedHash, ctrlBedHash);

  int rtsNum = outputRtsSites(paraInfo, rtsSiteList, outfp);
  fprintf(stderr, "#find %d Nm sites\n", rtsNum);

  treatReader.Close();
  ctrlReader.Close();
  if (paraInfo->verbose)  fprintf(stderr, "#free memory...\n");

  freeChromBed12Map(treatBedHash);
  freeChromBed12Map(ctrlBedHash);
  freeRtsSiteVector(rtsSiteList);
  freeChromStrandMap(bed12GeneHash);
}

void outputHeader(struct parameterInfo *paraInfo, FILE *outfp)
{
  // output the head line
  fprintf(outfp, "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\t");
  if (paraInfo->geneModel) fprintf(outfp, "geneName\tgeneStart\tgeneEnd\t");
  fprintf(outfp, "modifiedBase\t");
  fprintf(outfp, "endReadNum\tendRPM\t");
  fprintf(outfp, "upFC\tupCtrlFC\t");
  fprintf(outfp, "downFC\tdownCtrlFC\t");
  fprintf(outfp, "extendSeq\n");
}

int outputRtsSites(struct parameterInfo *paraInfo, rtsSiteVector &rtsSiteList, FILE *outfp)
{
  int siteNum = 0;
  // sort the list
  sort(rtsSiteList.begin(), rtsSiteList.end(), cmpScorePval);

  // remove sites with same position
  rtsSiteMap rtsSiteHash;
  for (rtsSiteVector::iterator vecItr = rtsSiteList.begin(); vecItr != rtsSiteList.end(); vecItr++)
  {
    rtsSiteInfo *rts = *vecItr;
    string key = rts->chrom + NumberToString(rts->chromStart) + rts->strand;
    rtsSiteHash[key] = rts;
  }

  // push sites into new siteVector
  rtsSiteVector newSiteList;
  for (rtsSiteMap::iterator it = rtsSiteHash.begin(); it != rtsSiteHash.end(); ++it)
  {
    rtsSiteInfo *rts = it->second;
    newSiteList.push_back(rts);
  }
  // sort the newSiteList
  sort(newSiteList.begin(), newSiteList.end(), cmpScorePval);

  /* don't remove the sites with same sequences
  // remove sites with same sequence
  rtsSiteMap rtsUniqSeqHash;
  for (rtsSiteVector::iterator vecItr = newSiteList.begin(); vecItr != newSiteList.end(); vecItr++)
  {
    rtsSiteInfo *rts = *vecItr;
    string key = rts->extSeq;
    rtsUniqSeqHash[key] = rts;
  }*/

  for (rtsSiteVector::iterator vecItr = newSiteList.begin(); vecItr != newSiteList.end(); vecItr++)
  {
    rtsSiteInfo *rts = *vecItr;
    if (rts->modPos <= 1) continue;
    siteNum += 1;
    fprintf(outfp, "%s\t%d\t%d\tendSeeker-%d\t%.5f\t%c\t", rts->chrom, rts->chromStart, rts->chromEnd, siteNum, rts->score, rts->strand);
    if (paraInfo->geneModel) fprintf(outfp, "%s\t%d\t%d\t", rts->geneName, rts->geneStart, rts->geneEnd);
    fprintf(outfp, "%c\t", rts->stopBase);
    fprintf(outfp, "%.5f\t%.5f\t", rts->tsNum, rts->tsRpm);
    fprintf(outfp, "%.5f\t%.5f\t", rts->tUpFC, rts->upCtrlFC);
    fprintf(outfp, "%.5f\t%.5f\t", rts->tDownFC, rts->downCtrlFC);
    fprintf(outfp, "%s\n", rts->extSeq);
    fflush(outfp);
  }

  return siteNum;
}

void callEndSites(struct parameterInfo *paraInfo, chromSizeMap &mapSize,
                  chromStrandBed12Map &bed12GeneHash,
                  FILE *genomefp, faidxMap &faiHash,
                  chromBed12Map &treatBedHash, chromBed12Map &ctrlBedHash)
{
  chromBed12Map::iterator it;
  char strand  = '.';
  int tReadNum = 0;
  int cReadNum = 0;
  int rtsNum   = 1;
  for (it = treatBedHash.begin(); it != treatBedHash.end(); ++it)
  {
    bed12Vector tBedList = it->second;
    bed12Vector cBedList = ctrlBedHash[it->first];
    char *chrom = (char *)it->first.c_str();
    string chromStr(chrom);
    if (faiHash.find(chromStr) == faiHash.end())
    {
      fprintf(stderr, "can't not find the chromosome %s, skip it.\n", chrom);
      continue;
    }
    faidx *fai   = faiHash[chromStr];
    int chromLen = mapSize[it->first];
    if (tBedList.size() >= 1 && cBedList.size() >= 1)
    {
      if (paraInfo->strand)
      {
        strand = '+';
        if (strand == '+')
        {
          profileInfo *treatProfile = (profileInfo *)safeMalloc(sizeof(struct profileInfo));
          profileInfo *ctrlProfile = (profileInfo *)safeMalloc(sizeof(struct profileInfo));
          tReadNum = getChromReadVals(paraInfo, tBedList, treatProfile, chromLen, strand);
          fprintf(stderr, "#treatment sample: get %d reads from chrom=%s and strand=%c\n", tReadNum, chrom, strand);
          cReadNum = getChromReadVals(paraInfo, cBedList, ctrlProfile, chromLen, strand);
          fprintf(stderr, "#control sample: get %d reads from chrom=%s and strand=%c\n", cReadNum, chrom, strand);

          if (paraInfo->geneModel == 1 || paraInfo->geneModel == 2)
          {
            fprintf(stderr, "#identifying modification sites in genes from chrom %s and strand: %c\n", chrom, strand);
            rtsNum  = getGeneEndSites(paraInfo, genomefp, fai, chrom, strand, chromLen, treatProfile, ctrlProfile, bed12GeneHash);
          }
          if (paraInfo->geneModel == 0 || paraInfo->geneModel == 2)
          {
            fprintf(stderr, "#identifying modification sites in chromosomes from chrom %s and strand: %c\n", chrom, strand);
            rtsNum = getChromEndSites(paraInfo, genomefp, fai, chrom, strand, chromLen, treatProfile, ctrlProfile);
          }
          freeProfiles(treatProfile);
          freeProfiles(ctrlProfile);
        }
        strand = '-';
        if (strand == '-')
        {
          profileInfo *treatProfile = (profileInfo *)safeMalloc(sizeof(struct profileInfo));
          profileInfo *ctrlProfile = (profileInfo *)safeMalloc(sizeof(struct profileInfo));
          tReadNum = getChromReadVals(paraInfo, tBedList, treatProfile, chromLen, strand);
          fprintf(stderr, "#treatment sample: get %d reads from chrom=%s and strand=%c\n", tReadNum, chrom, strand);
          cReadNum = getChromReadVals(paraInfo, cBedList, ctrlProfile, chromLen, strand);
          fprintf(stderr, "#control sample: get %d reads from chrom=%s and strand=%c\n", cReadNum, chrom, strand);
          if (paraInfo->geneModel == 1 || paraInfo->geneModel == 2)
          {
            fprintf(stderr, "#identifying modification sites in genes from chrom=%s and strand=%c\n", chrom, strand);
            rtsNum  = getGeneEndSites(paraInfo, genomefp, fai, chrom, strand, chromLen, treatProfile, ctrlProfile, bed12GeneHash);
          }
          if (paraInfo->geneModel == 0 || paraInfo->geneModel == 2)
          {
            fprintf(stderr, "#identifying modification sites in chromosomes from chrom=%s and strand=%c\n", chrom, strand);
            rtsNum = getChromEndSites(paraInfo, genomefp, fai, chrom, strand, chromLen, treatProfile, ctrlProfile);
          }
          freeProfiles(treatProfile);
          freeProfiles(ctrlProfile);
        }
      } // if strand
    }// if sizes
  } // for bed hash
}

int getGeneEndSites(struct parameterInfo *paraInfo,
                    FILE *gfp, faidx *fai,
                    char *chrom, char strand, int chromLen,
                    profileInfo *treatProfile, profileInfo *ctrlProfile,
                    chromStrandBed12Map &bed12GeneHash)
{
  int i = 0;
  int extLen = 15;
  double *taProfile  = treatProfile->heightProfile;
  double *caProfile  = ctrlProfile->heightProfile;
  double *teProfile  = treatProfile->endProfile;
  double *ceProfile  = ctrlProfile->endProfile;

  string chromStrand = chrom;
  chromStrand = chromStrand + strand;
  if (bed12GeneHash.find(chromStrand) == bed12GeneHash.end())
  {
    fprintf(stderr, "#can't found the chromosome %s and strand %c in gene file, skip it.\n", chrom, strand);
    return 0;
  }
  bed12Vector bed12List = bed12GeneHash[chromStrand];
  for (bed12Vector::iterator vecItr = bed12List.begin(); vecItr != bed12List.end(); vecItr++)
  {
    CBed12 *bed12 = *vecItr;

    regionLen *regLen = getRegionLen(bed12);
    int geneLen = regLen->geneLen;
    if (geneLen < 15)
    {
      fprintf(stderr, "#the length(%d) of gene: %s less than 15 nt, skip it\n", geneLen, bed12->name);
      safeFree(regLen);
      continue;
    }

    profileInfo *treatGeneProfile = (profileInfo *)safeMalloc(sizeof(struct profileInfo));
    profileInfo *ctrlGeneProfile  = (profileInfo *)safeMalloc(sizeof(struct profileInfo));

    double tReadVal  = getGeneVals(paraInfo, treatProfile, treatGeneProfile, bed12, geneLen);
    double cReadVal  = getGeneVals(paraInfo, ctrlProfile, ctrlGeneProfile, bed12, geneLen);

    if (tReadVal < paraInfo->minTag || cReadVal < paraInfo->minTag)
    {
      freeProfiles(treatGeneProfile);
      freeProfiles(ctrlGeneProfile);
      safeFree(regLen);
      continue;
    }

    char *geneSeq = getGeneSeq(gfp, fai, bed12, geneLen);
    findStopSites(paraInfo, chrom, geneLen, geneSeq, bed12, treatGeneProfile, ctrlGeneProfile);

    // free memory
    safeFree(regLen);
    safeFree(geneSeq);
    freeProfiles(treatGeneProfile);
    freeProfiles(ctrlGeneProfile);
  } // for end
  return 1;
}

int findStopSites(struct parameterInfo *paraInfo, char *chrom,
                  int geneLen, char *geneSeq, CBed12 *bed12,
                  profileInfo *treatProfile, profileInfo *ctrlProfile)
{
  int i = 0;
  int extLen = EXTEND_LEN;
  int rtsNum = 0;
  double *taProfile  = treatProfile->heightProfile; // treatment sample
  double *caProfile  = ctrlProfile->heightProfile; // control sample
  double *teProfile  = treatProfile->endProfile; // treatment sample
  double *ceProfile  = ctrlProfile->endProfile; // control sample
  geneReadInfo *geneInfo = getGeneReadInfo(treatProfile, ctrlProfile, 0, geneLen);
  if (geneInfo->tsSiteNum < 1)
  {
    safeFree(geneInfo);
    return 0;
  }
  for (i = 1; i < geneLen - 1; i++)
  {
    int endPos   = i;
    double tsNum = teProfile[endPos] + pseudCount; //reads whose 3’-ends are located at surveyed terminated sites in treatment sample
    if (tsNum < paraInfo->minTag) continue;
    double csNum = ceProfile[endPos] + pseudCount; //reads whose 3’-ends are located at surveyed terminated sites in control sample
    int prePos   = endPos - 1;
    int afterPos = endPos + 1;

    double taNum = taProfile[endPos] + pseudCount; // read coverage at surveyed terminated sites in treatment sample
    double caNum = caProfile[endPos] + pseudCount; // read coverage at surveyed terminated sites in control sample

    double tspNum = teProfile[prePos] + pseudCount; // reads whose 3’-ends are located one nucleotide upstream at surveyed terminated sites
    double cspNum = ceProfile[prePos] + pseudCount; // reads whose 3’-ends are located one nucleotide upstream at surveyed terminated sites

    double tapNum = teProfile[afterPos] + pseudCount; // reads whose 3’-ends are located one nucleotide downstream at surveyed terminated sites
    double capNum = ceProfile[afterPos] + pseudCount; // reads whose 3’-ends are located one nucleotide downstream at surveyed terminated sites

    double tsRpm = getRpmVals(tsNum, paraInfo->tTotalNum);
    double csRpm = getRpmVals(csNum, paraInfo->cTotalNum);

    int idxStart = endPos - paraInfo->windowSize;
    int idxEnd   = endPos + paraInfo->windowSize;
    if (idxStart < 0) idxStart = 0;
    if (idxEnd > geneLen) idxEnd = geneLen;

    geneReadInfo *geneReionInfo = getGeneReadInfo(treatProfile, ctrlProfile, idxStart, idxEnd);
    double rcsm     = geneReionInfo->csTotalNum / geneReionInfo->csSiteNum; // control start mean number in region

    double tUpFC    = tsNum / tspNum;
    double cUpFC    = csNum / rcsm;
    double upCtrlFC = tUpFC / cUpFC;

    double tDownFC  = tsNum / tapNum;
    double cDownFC  = csNum / rcsm;
    double downCtrlFC  = tDownFC / cDownFC;

    safeFree(geneReionInfo);

    if (tUpFC >= paraInfo->fold
        && tDownFC >= paraInfo->fold
        && upCtrlFC >= paraInfo->fold
        && downCtrlFC >= paraInfo->fold)
    {
      int modPos = endPos;
      if (paraInfo->type == 2) modPos = prePos;
      rtsSiteInfo *stopInfo = (rtsSiteInfo *)safeMalloc(sizeof(rtsSiteInfo));
      stopInfo->siteModel  = 1;
      stopInfo->chrom      = strClone(chrom);
      stopInfo->chromStart = getGenomicPos(modPos, bed12);
      stopInfo->chromEnd   = stopInfo->chromStart + 1;
      stopInfo->geneStart  = modPos;
      stopInfo->geneEnd    = stopInfo->geneStart + 1;
      stopInfo->geneName   = strClone(bed12->name);
      stopInfo->score      = tUpFC;
      stopInfo->strand     = bed12->strand;
      stopInfo->tsNum      = tsNum;
      stopInfo->tsRpm      = tsRpm;
      stopInfo->tUpFC      = tUpFC;
      stopInfo->cUpFC      = cUpFC;
      stopInfo->upCtrlFC   = upCtrlFC;
      stopInfo->tDownFC    = tDownFC;
      stopInfo->cDownFC    = cDownFC;
      stopInfo->downCtrlFC = downCtrlFC;

      stopInfo->stopBase = geneSeq[modPos];

      int extStart = modPos - extLen;
      int extEnd   = modPos + extLen + 1;
      if (extStart < 0) extStart = 0;
      if (extEnd > geneLen) extEnd = geneLen;
      int seqLen = extEnd - extStart;
      char *extSeq = (char*)safeMalloc(sizeof(char) * (seqLen + 1));
      strncpy(extSeq, geneSeq + extStart, seqLen);
      int endIdx = modPos - extStart;
      extSeq[endIdx] = 'm';
      stopInfo->extSeq = extSeq;
      stopInfo->modPos = endIdx;

      if (filterCandidate(paraInfo, stopInfo))
      {
        freeRtsSite(stopInfo);
        continue;
      }
      else
      {
        rtsSiteList.push_back(stopInfo);
        rtsNum++;
      }
    }// pfold
  } // for loop for gene length
  safeFree(geneInfo);
  return rtsNum;
}

int filterCandidate(struct parameterInfo *paraInfo, rtsSiteInfo *rts)
{
  int filterTag = 1;

  if (rts->tUpFC   < paraInfo->fold) return filterTag;
  if (rts->tDownFC < paraInfo->fold) return filterTag;
  if (rts->upCtrlFC < paraInfo->fold) return filterTag;
  if (rts->downCtrlFC < paraInfo->fold) return filterTag;

  filterTag = 0;
  return filterTag;
}

/* for identifying Nm sites at genome-wide scale */
int getChromEndSites(struct parameterInfo *paraInfo,
                     FILE *genomefp, faidx *fai,
                     char *chrom, char strand, int chromLen,
                     profileInfo *treatProfile, profileInfo *ctrlProfile)
{
  int i = 0;
  int rtsNum = 0;
  int extLen = EXTEND_LEN;
  double *taProfile  = treatProfile->heightProfile;
  double *caProfile  = ctrlProfile->heightProfile;
  double *teProfile  = treatProfile->endProfile;
  double *ceProfile  = ctrlProfile->endProfile;
  for (i = 1; i < chromLen; i++)
  {
    int endPos   = i;
    double tsNum = teProfile[endPos] + pseudCount; // treatment start
    double tsRpm = getRpmVals(tsNum, paraInfo->tTotalNum);
    if (tsNum < paraInfo->minTag || tsRpm < paraInfo->rpm) continue;
    double csNum = ceProfile[endPos] + pseudCount; // control start
    int prePos   = endPos - 1;
    int afterPos = endPos + 1;

    if (strand == '-')
    {
      prePos = endPos + 1;
      afterPos = endPos - 1;
    }

    double taNum = taProfile[endPos] + pseudCount;
    double caNum = caProfile[endPos] + pseudCount;

    double tspNum = teProfile[prePos] + pseudCount;
    double cspNum = ceProfile[prePos] + pseudCount;

    double tapNum = teProfile[afterPos] + pseudCount;
    double capNum = ceProfile[afterPos] + pseudCount;

    double csRpm = getRpmVals(csNum, paraInfo->cTotalNum);

    double tUpFC = tsNum / tspNum;
    double cUpFC = csNum / cspNum;
    double upCtrlFC = tUpFC / cUpFC;

    double tDownFC = tsNum / tapNum;
    double cDownFC = csNum / capNum;
    double downCtrlFC = tDownFC / cDownFC;

    if (tUpFC >= paraInfo->fold
        && tDownFC >= paraInfo->fold
        && upCtrlFC >= paraInfo->fold
        && downCtrlFC >= paraInfo->fold)
    {
      int modPos = endPos;
      if (paraInfo->type == 2) modPos = prePos;
      rtsSiteInfo *stopInfo = (rtsSiteInfo *)safeMalloc(sizeof(rtsSiteInfo));
      stopInfo->siteModel  = 1;
      stopInfo->chrom      = strClone(chrom);
      stopInfo->chromStart = modPos;
      stopInfo->chromEnd   = stopInfo->chromStart + 1;
      stopInfo->geneStart  = modPos;
      stopInfo->geneEnd    = stopInfo->geneStart + 1;
      stopInfo->geneName   = strClone(chrom);
      stopInfo->score      = tUpFC;
      stopInfo->strand     = strand;
      stopInfo->tsNum      = tsNum;
      stopInfo->tsRpm      = tsRpm;
      stopInfo->tUpFC      = tUpFC;
      stopInfo->cUpFC      = cUpFC;
      stopInfo->upCtrlFC   = upCtrlFC;
      stopInfo->tDownFC    = tDownFC;
      stopInfo->cDownFC    = cDownFC;
      stopInfo->downCtrlFC = downCtrlFC;

      char *rtsSeq = faidxFetchSeq(genomefp, fai, modPos, modPos + 1, strand);
      stopInfo->stopBase = rtsSeq[0];
      safeFree(rtsSeq);

      int extStart = modPos - extLen;
      int extEnd   = modPos + extLen + 1;
      if (extStart < 0) extStart = 0;
      if (extEnd > fai->len) extEnd = fai->len;
      char *extSeq = faidxFetchSeq(genomefp, fai, extStart, extEnd, strand);
      int endIdx = modPos - extStart;
      extSeq[endIdx] = 'm';
      stopInfo->extSeq = extSeq;
      stopInfo->modPos = endIdx;

      if (filterCandidate(paraInfo, stopInfo))
      {
        freeRtsSite(stopInfo);
        continue;
      }
      else
      {
        rtsSiteList.push_back(stopInfo);
        rtsNum++;
      }

    }// pfold
  } // for loop for chromosome length
  return rtsNum;
}

double getGeneVals(struct parameterInfo *paraInfo, profileInfo *profile,
                   profileInfo *geneProfile, CBed12 *bed12, int geneLen)
{
  int i = 0;
  int j = 0;
  int k = 0;
  int idx = 0;
  double geneVal = 0;
  geneProfile->endProfile    = (double *)safeMalloc(sizeof(double) * geneLen);
  geneProfile->heightProfile   = (double *)safeMalloc(sizeof(double) * geneLen);
  for (i = 0; i < bed12->blockCount; i++)
  {
    if (bed12->strand == '+')
    {
      idx = i;
    }
    else  // reverse the all items
    {
      idx = bed12->blockCount - i - 1;
    }
    int chromStart = bed12->chromStart + bed12->chromStarts[idx];
    int chromEnd = chromStart + bed12->blockSizes[idx];
    if (bed12->strand == '+')
    {
      for (j = chromStart; j < chromEnd; j++)
      {
        geneProfile->endProfile[k]  = profile->endProfile[j];
        geneProfile->heightProfile[k] = profile->heightProfile[j];
        geneVal += profile->heightProfile[j];
        k++;
      }
    }// if + strand
    else
    {
      for (j = chromEnd - 1; j >= chromStart; j--)
      {
        geneProfile->endProfile[k]  = profile->endProfile[j];
        geneProfile->heightProfile[k] = profile->heightProfile[j];
        geneVal += profile->heightProfile[j];
        k++;
      }
    } // else - strand
  } // for count
  return geneVal;
}

char *getGeneSeq(FILE *gfp, faidx *fai, CBed12 *bed12, int geneLen)
{
  int i = 0;
  int idx = 0;
  char *seq = NULL;
  seq = (char *)safeMalloc(sizeof(char) * (geneLen + 1));
  for (i = 0; i < bed12->blockCount; i++)
  {
    if (bed12->strand == '+')
    {
      idx = i;
    }
    else  // reverse the all items
    {
      idx = bed12->blockCount - i - 1;
    }
    int chromStart = bed12->chromStart + bed12->chromStarts[idx];
    int chromEnd = chromStart + bed12->blockSizes[idx];
    char *subSeq = faidxFetchSeq(gfp, fai, chromStart, chromEnd, bed12->strand);
    strcat(seq, subSeq);
    safeFree(subSeq);
  }
  seq[geneLen] = '\0';
  return seq;
}

regionLen *getRegionLen(CBed12 *bed12)
{
  int blockCount, i;
  int geneLen  = 0;
  int utr5Len  = 0;
  int utr3Len  = 0;
  int cdsLen   = 0;
  regionLen *reg = (regionLen *)safeMalloc(sizeof(struct regionLenInfo));
  int thickStart = bed12->thickStart;
  int thickEnd   = bed12->thickEnd;
  int orgStart   = bed12->chromStart;
  int orgEnd     = bed12->chromEnd;
  char strand    = bed12->strand;
  blockCount = bed12->blockCount;
  for (i = 0; i < blockCount; i++)
  {
    int chromStart = bed12->chromStart + bed12->chromStarts[i];
    int chromEnd = chromStart + bed12->blockSizes[i];
    geneLen += bed12->blockSizes[i];
    int utr5Ovl = overlapLength(orgStart, thickStart, chromStart, chromEnd);
    if (utr5Ovl > 0)
    {
      utr5Len += utr5Ovl;
    }
    int utr3Ovl = overlapLength(thickEnd, orgEnd, chromStart, chromEnd);
    if (utr3Ovl > 0)
    {
      utr3Len += utr3Ovl;
    }
    int cdsOvl = overlapLength(thickStart, thickEnd, chromStart, chromEnd);
    if (cdsOvl > 0)
    {
      cdsLen += cdsOvl;
    }
  } // for blockCount
  if (strand == '-')
  {
    int tmpUtr5Len = utr5Len;
    utr5Len = utr3Len;
    utr3Len = tmpUtr5Len;
  }
  reg->geneLen = geneLen;
  reg->utr5Len = utr5Len;
  reg->utr3Len = utr3Len;
  reg->cdsLen  = cdsLen;
  return reg;
}

geneReadInfo *getGeneReadInfo(profileInfo *treatProfile, profileInfo *ctrlProfile, int start, int end)
{
  int i = 0;
  double *taProfile  = treatProfile->heightProfile;
  double *caProfile  = ctrlProfile->heightProfile;
  double *teProfile  = treatProfile->endProfile;
  double *ceProfile  = ctrlProfile->endProfile;
  geneReadInfo *geneInfo = (geneReadInfo*) safeMalloc(sizeof(geneReadInfo));

  geneInfo->tsTotalNum = pseudCount;
  geneInfo->csTotalNum = pseudCount;
  geneInfo->taTotalNum = pseudCount;
  geneInfo->caTotalNum = pseudCount;
  geneInfo->tsSiteNum = pseudCount;
  geneInfo->csSiteNum = pseudCount;

  for (i = start; i < end; i++)
  {
    int endPos   = i;
    double tsNum = teProfile[endPos];
    double csNum = ceProfile[endPos];
    if (tsNum <= 0 && csNum <= 0) continue;

    double taNum = taProfile[endPos];
    double caNum = caProfile[endPos];

    geneInfo->tsTotalNum += tsNum;
    geneInfo->csTotalNum += csNum;
    geneInfo->taTotalNum += taNum;
    geneInfo->caTotalNum += caNum;
    if (tsNum > 0) geneInfo->tsSiteNum += 1;
    if (csNum > 0) geneInfo->csSiteNum += 1;
  }
  return geneInfo;
}

int getGenomicPos(int pos, CBed12 *bed12)
{
  int blockCount, i;
  int geneLen  = 0;
  int blockStart  = 0;
  int blockEnd = 0;
  int genomePos = 0;
  char strand = bed12->strand;
  blockCount = bed12->blockCount;
  if (bed12->strand == '+') // plus strand
  {
    for (i = 0; i < blockCount; i++)
    {
      int chromStart = bed12->chromStart + bed12->chromStarts[i];
      int chromEnd = chromStart + bed12->blockSizes[i];
      blockEnd += bed12->blockSizes[i];
      blockStart = blockEnd - bed12->blockSizes[i];
      if (pos >= blockStart && pos < blockEnd)
      {
        genomePos = chromStart + (pos - blockStart);
        break;
      }
    }// for blockCount
  } // if end
  else // minus strand
  {
    for (i = blockCount - 1; i >= 0; i--)
    {
      int chromStart = bed12->chromStart + bed12->chromStarts[i];
      int chromEnd = chromStart + bed12->blockSizes[i];
      blockEnd += bed12->blockSizes[i];
      blockStart = blockEnd - bed12->blockSizes[i];
      if (pos >= blockStart && pos < blockEnd)
      {
        genomePos = chromStart + (blockEnd - pos) - 1;
        break;
      }
    }// for blockCount
  } // if end
  return genomePos;
}

int getChromReadVals(struct parameterInfo *paraInfo, bed12Vector &bedList,
                     profileInfo *profile, int chromLen, char strand)
{
  int i = 0;
  int tagNum = 0;
  profile->endProfile = (double *)safeMalloc(sizeof(double) * chromLen);
  profile->heightProfile = (double *)safeMalloc(sizeof(double) * chromLen);
  for (i = 0; i < chromLen; i++)
  {
    profile->endProfile[i] = 0;
    profile->heightProfile[i] = 0;
  }

  //fprintf(stderr, "# safeZeroedMalloc for chrom length %d\n", chromLen);
  for (bed12Vector::iterator vecItr = bedList.begin(); vecItr != bedList.end(); vecItr++)
  {
    CBed12 *bed = *vecItr;
    if (paraInfo->strand) // skip other strand
    {
      if (bed->strand != strand)
      {
        continue;
      }
    }
    int start = bed->chromStart;
    int end = bed->chromEnd;

    int readLen = 0;
    for (i = 0; i < bed->blockCount; i++)
    {
      readLen += bed->blockSizes[i];
    }
    if (readLen < paraInfo->minLen || readLen > paraInfo->maxLen) continue;

    int startIdx = start;
    int endIdx = end - 1;
    if (bed->strand == '-')
    {

      startIdx = end - 1;
      endIdx = start;
    }
    profile->endProfile[endIdx] += bed->score;

    // for coverage profilings
    for (i = 0; i < bed->blockCount; i++) // modified and remove the 1
    {
      int bStart = bed->chromStart + bed->chromStarts[i];
      int bEnd   = bStart + bed->blockSizes[i];
      int j = 0;
      for (j = bStart; j < bEnd; j++)
      {
        profile->heightProfile[j] += bed->score;
      }
    }
    tagNum++;
  } // for end
  return tagNum;
}

double getRpmVals(double readNum, double totalNum)
{
  double rpm = readNum * 1000000 / totalNum;
  return rpm;
}

void freeProfiles(profileInfo *profile)
{
  safeFree(profile->endProfile);
  safeFree(profile->heightProfile);
  safeFree(profile);
}

double round(double r)
{
  return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

void freeRtsSite(rtsSiteInfo *rts)
{
  safeFree(rts->chrom);
  safeFree(rts->geneName);
  safeFree(rts->extSeq);
  safeFree(rts);
}

void freeRtsSiteVector(rtsSiteVector &rtsSiteList)
{
  for (rtsSiteVector::iterator vecItr = rtsSiteList.begin(); vecItr != rtsSiteList.end(); vecItr++)
  {
    rtsSiteInfo *rts = *vecItr;
    freeRtsSite(rts);
  }
  rtsSiteList.clear();
}
