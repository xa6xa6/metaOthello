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
	if(argc != 4)
	{
		cout << "Executable inputChromosomeSeq2taxoInfoFile outputFa2taxoInfoFile outputDir_fa" << endl;
		exit(1);
	}
	string inputChromosomeSeq2taxoInfoFile = argv[1];
	string outputMergedFa2taxoInfoFile = argv[2];
	string outputDir_fa = argv[3];

	ChromosomeSeq_Info_Vec tmpChrSeqInfoVec;
	tmpChrSeqInfoVec.initiate_mergedSpeciesFa(inputChromosomeSeq2taxoInfoFile);
	tmpChrSeqInfoVec.initiate_fa_fileVec_creatDir_mergedSpeciesFa(outputDir_fa);
	tmpChrSeqInfoVec.printMergedSpecies_fa_info(outputMergedFa2taxoInfoFile);
	tmpChrSeqInfoVec.generateMergedSpecies_fa();
	return 0;
}