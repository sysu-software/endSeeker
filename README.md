# endSeeker
endSeeker: A computational software for identifying 2'-O-Methylation sites from Nm-REP-seq data.

Overview:
---------
endSeeker is an software to identify novel Nm sites by calculating the 3'-end coverage of the candidate Nm sites from Nm-REP-seq data. 

Usage:
---------
Usage:  endSeeker [options] --fa <genome seq> --fai <fai file> --gene <bed12 file> --treat <alignments> --input <input alignments><BR>
endSeeker: A computational software for identifying 2'-O-Methylation sites from Nm-REP-seq data.<BR>
[options]<BR>
--fa <string>          : genome sequence<fasta format>[Required]<BR>
--fai <string>         : genome fai file<fai format>[Required].<BR>
                         using "samtools faidx" to generate fai file<BR>
--gene <string>        : gene file <BED12 format>[Required]<BR>
--treat <string>       : file treated by mgR/mgR+OED file<BAM format>[Required]<BR>
--input <string>       : input file<BAM format>[Required]<BR>
-v/--verbose           : verbose information<BR>
-V/--version           : endSeeker version<BR>
-h/--help              : help informations<BR>
-s/--strand            : strand-specific sequencing data<BR>
-n/--norm              : normalized reads to the locus number<BR>
-c/--collapser         : keep duplication, deault is false<BR>
-o/--outfile <string>  : output file<BR>
-t/--min-tag <double>  : minimum tag number for each end site, default>=5.0 read<BR>
-r/--rpm <double>      : minimum rpm value for each end site, default>=0.001<BR>
-f/--fold <int>        : minimum fold-change[default>=1.0]<BR>
-w/--window <int>      : window size around the end position[default=20]<BR>
-l/--min-len <int>     : minimum length of reads, default=15<BR>
-L/--max-len <int>     : maximum length of reads, default=1000<BR>


Installation:<BR>
---------
Download endSeeker-0.1.tar.gz from https://github.com/sysu-software/endSeeker/releases ; unpack it, and make:<BR>
tar -xzvf endSeeker-0.1.tar.gz<BR>
cd endSeeker-0.1<BR>
make<BR>

System requirements:
---------
Operating system: endSeeker is designed to run on POSIX-compatible platforms, including UNIX, Linux and Mac OS/X. We have tested  most extensively on Linux and MacOS/X because these are the machines we develop on.<BR>
Compiler: The source code is compiled with  the C++ compiler g++. We test the code using the g++ compilers.<BR>
Libraries and other installation requirements: endSeeker includes one software library: the BamTools library package. All will automatically compile during endSeeker installation process.<BR>
By default, endSeeker does not require any additional libraries to be installed by you.<BR>

Prerequisites:<BR>
---------
Dependencies: The input of endSeeker is BAM file. So you need the short read mapper STAR or other mappers<BR>
You can get the most fresh versions:<BR>
(1)	STAR: https://github.com/alexdobin/STAR<BR>
(2)	Samtools: http://www.htslib.org/<BR>
You need to have the reference genome, fai file, annotation file, and  STAR indexes for genome and annotation.<BR>
You can constructed these datasets by yourself using following steps:<BR>
As an example, let's assume you use human genome (version hg38).<BR>
(1)	Genome:<BR>
mkdir genome<BR>
wget -c 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'<BR>
gzip -d hg38.fa.gz<BR>
cd ..<BR>
(2)	Annotation:<BR>
You can use Table Browser to get the BED12 for genome annotation(e.g. GENCODE)<BR>
http://genome.ucsc.edu/cgi-bin/hgTables<BR>
e.g. You can save the output files in genome directory: hg38.gencode.bed12<BR>
(3) Build the genome index and align reads to genome<BR>
STAR --runMode genomeGenerate --genomeDir ./starIndex --genomeFastaFiles ./genome/hg38.fa --sjdbGTFfile gencode.v30.annotation.gtf --sjdbOverhang 100<BR>
(4)build the fai index:<BR>
samtools faidx hg38.fa<BR><BR>
(5)Align reads to genome using STAR<BR>
STAR parameters as follows: --alignEndsType EndToEnd --outFilterType BySJout --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 30 --outFilterMismatchNmax 15 --outFilterMismatchNoverLmax 0.1 --outFilterScoreMin 0 --outFilterScoreMinOverLread 0 --outFilterMatchNmin 15 --outFilterMatchNminOverLread 0.8 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --seedSearchStartLmax 15 --seedSearchStartLmaxOverLread 1 --seedSearchLmax 0 --seedMultimapNmax 20000 --seedPerReadNmax 1000 --seedPerWindowNmax 100 --seedNoneLociPerWindow 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts<BR>

run endSeeker:
---------
cd example;
../endSeeker --norm -t 5 -f 3 --fa chr21.fa --fai chr21.fa.fai --gene human_rRNA_genes.bed12 --treat MgR_treatment_sample.bam --input MgR_control_sample.bam \>endSeeker_candidate_Nm_sites.txt<BR>

Output:
---------
#chrom	chromStart	chromEnd	name	score	strand	geneName	geneStart	geneEnd	modifiedBase	endReadNum	endRPM	upFC	upCtrlFC	downFC	downCtrlFC	extendSeq<BR>
chr21	8217780	8217781	endSeeker-1	3.06000	+	NR_146148|28S|rRNA	3894	3895	C	306.00000	135.37146	3.06000	3.99876	3.97403	5.19320	AGCGGGGAAAGAAGAmCCTGTTGAGCTTGAC<BR>

Note: # is comment line<BR>

Acknowledgements:
---------
Thanks a lot to everyone who contributed to the public code (e.g. BamTools, Samtools) used by endSeeker.<BR>

Contact :
---------
*****************************************************************************************<BR>
 \*	endSeeker - A computational software for identifying 2'-O-Methylation sites from Nm-REP-seq data.<BR>
 \*<BR>
 \*	Author : Jian-Hua Yang <yangjh7@mail.sysu.edu.cn><BR>
 \* <BR>
 \*	RNA Information Center, School of Life Sciences, Sun Yat-Sen University<BR>
 \*	<BR>
 \*  Create date: 11/18/2019<BR>
 \*  <BR>
 \*  last modified time: 09/01/2020<BR>
 ****************************************************************************************<BR>
