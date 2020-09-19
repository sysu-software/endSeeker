/*
   Copyright (c) 2019, 2020 by Jianhua Yang <yangjh7@mail.sysu.edu.cn>
   endSeeker: A computational software for identifying 2’-O-Methylation sites from Nm-REP-seq data.
   Date: 2019/11/18 @ Sun Yat-sen University
*/
#ifndef endSeeker_HEAD_H
#define endSeeker_HEAD_H

#define pseudCount 1.0
#define EXTEND_LEN 15

typedef struct parameterInfo
{
  int    verbose;
  int    collapser;
  int    psi;
  int    minLen;
  int    maxLen;
  int    windowSize;
  int    strand;
  int    norm;
  int    geneModel;
  int    type;
  double minTag;
  double pvalue;
  double fdr;
  double mfold;
  double fold;
  double tTotalNum;
  double cTotalNum;
  double rpm;
  long   genomeSize;
} parameterInfo;

typedef struct profileInfo
{
  double *endProfile;
  double *heightProfile;
} profileInfo;

struct geneReadInfo
{
  double tsTotalNum;
  double csTotalNum;
  double taTotalNum;
  double caTotalNum;
  double tsSiteNum;
  double csSiteNum;
};
typedef struct geneReadInfo geneReadInfo;

typedef struct regionLenInfo
{
  int geneLen;
  int utr5Len;
  int utr3Len;
  int cdsLen;
} regionLen;

struct rtsSiteInfo
{
  int siteModel; /* sites from genes or chromosomes */
  char *chrom;
  int chromStart;
  int chromEnd;
  int geneStart;
  int geneEnd;
  int modPos;
  char *geneName;
  char *extSeq;
  double score;
  char strand;
  char stopBase;

  double tsNum;
  double tsRpm;

  double tUpFC;
  double cUpFC;

  double tDownFC;
  double cDownFC;

  double upCtrlFC;
  double downCtrlFC;
};
typedef struct rtsSiteInfo rtsSiteInfo;

typedef vector<rtsSiteInfo *> rtsSiteVector;

typedef map<string, rtsSiteInfo *> rtsSiteMap;

extern rtsSiteVector rtsSiteList;

void runEndSeeker(struct parameterInfo *paraInfo, FILE *genomefp,
                  FILE *faifp, FILE *bedfp,
                  char *treatBamFile, char *ctrlBamFile, FILE *outfp);

void outputHeader(struct parameterInfo *paraInfo, FILE *outfp);

int outputRtsSites(struct parameterInfo *paraInfo, rtsSiteVector &rtsSiteList, FILE *outfp);

void callEndSites(struct parameterInfo *paraInfo, chromSizeMap &mapSize,
                  chromStrandBed12Map &bed12GeneHash,
                  FILE *genomefp, faidxMap &faiHash,
                  chromBed12Map &treatBedHash, chromBed12Map &ctrlBedHash);

int getChromReadVals(struct parameterInfo *paraInfo, bed12Vector &bedList,
                     profileInfo *profile, int chromLen, char strand);

regionLen *getRegionLen(CBed12 *bed12);

double getGeneVals(struct parameterInfo *paraInfo, profileInfo *profile,
                   profileInfo *geneProfile, CBed12 *bed12, int geneLen);

geneReadInfo *getGeneReadInfo(profileInfo *treatProfile,
                              profileInfo *ctrlProfile,
                              int start, int end);

char *getGeneSeq(FILE *gfp, faidx *fai, CBed12 *bed12, int geneLen);

int findStopSites(struct parameterInfo *paraInfo, char *chrom,
                  int geneLen, char *geneSeq, CBed12 *bed12,
                  profileInfo *treatProfile, profileInfo *ctrlProfile);

int getGeneEndSites(struct parameterInfo *paraInfo,
                    FILE *gfp, faidx *fai,
                    char *chrom, char strand, int chromLen,
                    profileInfo *treatProfile, profileInfo *ctrlProfile,
                    chromStrandBed12Map &bed12GeneHash);

int getChromEndSites(struct parameterInfo *paraInfo,
                     FILE *genomefp, faidx *fai,
                     char *chrom, char strand, int chromLen,
                     profileInfo *treatProfile, profileInfo *ctrlProfile);

double getRpmVals(double readNum, double totalNum);

int getGenomicPos(int pos, CBed12 *bed12);

void freeProfiles(profileInfo *profile);

double round(double r);

void freeRtsSite(rtsSiteInfo *rts);

void freeRtsSiteVector(rtsSiteVector &rtsSiteList);

int filterCandidate(struct parameterInfo *paraInfo, rtsSiteInfo *rts);

int cmpScorePval(const rtsSiteInfo *x, const rtsSiteInfo *y);

#endif /* End endSeeker_HEAD_H */
