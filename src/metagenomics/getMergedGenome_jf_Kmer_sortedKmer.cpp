#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>
#include "../mps3Lib/read_block_test.h"
#include "../mps3Lib/otherFunc.h"
#include "../mps3Lib/index_info.h"
#include "general/chromosomeSeq_info_vec.h"
#include "general/species_info.h"
#include "general/genus_info.h"
#include "general/phylum_info.h"
time_t nowtime;
struct tm *local;

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 9)
	{
		cout << "Executable inputChromosomeSeq2taxoInfoFile outputInfo2taxoInfoFilePrefix_jf_Kmer_sortedKmer inputDir_fa JellyFishBinPath Kmer_length threads_num jellyFishBitMapSize outputDirPrefix_jf_Kmer_sortedKmer" << endl;
		exit(1);
	}
	string inputChromosomeSeq2taxoInfoFile = argv[1];

	string outputInfo2taxoInfoFilePrefix_jf_Kmer_sortedKmer = argv[2];
	string outputJf2taxoInfoFile = outputInfo2taxoInfoFilePrefix_jf_Kmer_sortedKmer + ".jf";
	string outputKmer2taxoInfoFile = outputInfo2taxoInfoFilePrefix_jf_Kmer_sortedKmer + ".Kmer";
	string outputSortedKmer2taxoInfoFile = outputInfo2taxoInfoFilePrefix_jf_Kmer_sortedKmer + ".sortedKmer";
	
	string inputDir_fa = argv[3];
	string JellyFishBinPath = argv[4];
	JellyFishBinPath += "/";
	string Kmer_length_str = argv[5];
	int Kmer_length = atoi(Kmer_length_str.c_str());
	string threads_num_str = argv[6];
	int threads_num = atoi(threads_num_str.c_str());
	string jellyFishBitMapSizeStr = argv[7];
	
	string outDirPrefix = argv[8];
	string outputDir_jf = outDirPrefix + "_jf/";
	string outputDir_Kmer = outDirPrefix + "_Kmer/";
	string outputDir_sortedKmer = outDirPrefix + "_sortedKmer/";
	string outputDir_tmpSortDir = outDirPrefix + "_tmpSortDir";
	string cmd_mkdir_tmpSortDir = "mkdir " + outputDir_tmpSortDir;
	system(cmd_mkdir_tmpSortDir.c_str());

	ChromosomeSeq_Info_Vec tmpChrSeqInfoVec;
	cout << "start to initiate inputChromosomeSeq2taxoInfoFile" << endl;
	tmpChrSeqInfoVec.initiate(inputChromosomeSeq2taxoInfoFile);
	cout << "start to initiate_fa_fileVec_keepDir" << endl;
	tmpChrSeqInfoVec.initiate_fa_fileVec_keepDir(inputDir_fa);
	cout << "start to initiate_jf_fileVec_creatDir" << endl;
	tmpChrSeqInfoVec.initiate_jf_fileVec_creatDir(outputDir_jf, Kmer_length);
	cout << "start to initiate_Kmer_fileVec_creatDir" << endl;
	tmpChrSeqInfoVec.initiate_Kmer_fileVec_creatDir(outputDir_Kmer, Kmer_length);
	cout << "start to initiate_sortedKmer_fileVec_creatDir" << endl;
	tmpChrSeqInfoVec.initiate_sortedKmer_fileVec_creatDir(outputDir_sortedKmer, Kmer_length);
	cout << "start to printMergedGenome_jf_info" << endl;
	tmpChrSeqInfoVec.printMergedGenome_jf_info(outputJf2taxoInfoFile, Kmer_length);
	cout << "start to printMergedGenome_Kmer_info" << endl;
	tmpChrSeqInfoVec.printMergedGenome_Kmer_info(outputKmer2taxoInfoFile, Kmer_length);
	cout << "start to printMergedGenome_sortedKmer_info" << endl;
	tmpChrSeqInfoVec.printMergedGenome_sortedKmer_info(outputSortedKmer2taxoInfoFile, Kmer_length);

	cout << "start to generateMergedGenome_jf" << endl;
	tmpChrSeqInfoVec.generateMergedGenome_jf(JellyFishBinPath, Kmer_length, threads_num, jellyFishBitMapSizeStr);
	cout << "start to generateMergedGenome_Kmer" << endl;
	tmpChrSeqInfoVec.generateMergedGenome_Kmer(JellyFishBinPath, Kmer_length);
	cout << "start to generateMergedGenome_sortedKmer" << endl;
	tmpChrSeqInfoVec.generateMergedGenome_sortedKmer(outputDir_tmpSortDir, Kmer_length);

	return 0;
}