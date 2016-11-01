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
time_t nowtime;
struct tm *local;

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 9)
	{
		cout << "Executable inputChromosomeSeq2taxoInfoFile outputJf2taxoInfoFile inputDir_fa JellyFishBinPath Kmer_length threads_num jellyFishBitMapSize outputDir_jf" << endl;
		exit(1);
	}
	string inputChromosomeSeq2taxoInfoFile = argv[1];
	string outputMergedJf2taxoInfoFile = argv[2];
	string inputDir_fa = argv[3];
	string JellyFishBinPath = argv[4];
	JellyFishBinPath += "/";
	string Kmer_length_str = argv[5];
	int Kmer_length = atoi(Kmer_length_str.c_str());
	string threads_num_str = argv[6];
	int threads_num = atoi(threads_num_str.c_str());
	string jellyFishBitMapSizeStr = argv[7];
	string outputDir_jf = argv[8];
	outputDir_jf += "/";

	ChromosomeSeq_Info_Vec tmpChrSeqInfoVec;
	tmpChrSeqInfoVec.initiate(inputChromosomeSeq2taxoInfoFile);
	tmpChrSeqInfoVec.initiate_fa_fileVec_keepDir(inputDir_fa);
	tmpChrSeqInfoVec.initiate_jf_fileVec_creatDir(outputDir_jf, Kmer_length);
	tmpChrSeqInfoVec.printMergedGenome_jf_info(outputMergedJf2taxoInfoFile, Kmer_length);
	tmpChrSeqInfoVec.generateMergedGenome_jf(JellyFishBinPath, Kmer_length, threads_num, jellyFishBitMapSizeStr);
	return 0;
}