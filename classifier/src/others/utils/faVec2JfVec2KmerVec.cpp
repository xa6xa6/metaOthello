// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
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
#include "../general/fa2jf_info.h"
#include "../general/jf2Kmer_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 10)
	{
		cout << "Executable#0" << endl;
		cout << "jellyFishiBin#1" << endl;
		cout << "threads_num#2" << endl;
		cout << "count_min#3" << endl;
		cout << "Kmer_length#4" << endl;
		cout << "bf_size#5" << endl;
		cout << "inputIdListFile#6" << endl;
		cout << "inputDir_fa#7" << endl;
		cout << "outputDir_jf#8" << endl;
		cout << "outputDir_Kmer#9" << endl;
		exit(1);
	}	

	string jellyFishiBin = argv[1];
	string threads_num_str = argv[2];
	string count_min_str = argv[3];
	string Kmer_length_str = argv[4];
	string bf_size_str = argv[5];
	string inputIdListFile = argv[6];
	string inputDir_fa = argv[7]; inputDir_fa += "/";
	string outputDir_jf = argv[8]; outputDir_jf += "/";
	string mkdir_jf = "mkdir " + outputDir_jf;
	system(mkdir_jf.c_str());
	string outputDir_Kmer = argv[9]; outputDir_Kmer += "/";
	string mkdir_Kmer = "mkdir " + outputDir_Kmer;
	system(mkdir_Kmer.c_str());

	vector<string> idVec;
	ifstream idList_ifs(inputIdListFile.c_str());
	while(!idList_ifs.eof())
	{
		string tmpStr;
		getline(idList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		idVec.push_back(tmpStr);
	}
	idList_ifs.close();
	cout << "start to do faVec 2 jfVec and jfVec 2 KmerVec" << endl;
	// faVec 2 jfVec
	for(int tmp = 0; tmp < idVec.size(); tmp++)
	{
		string tmpId = idVec[tmp];
		cout << endl << "tmpFileIndex: " << tmp << "; tmpId: " << tmpId << endl;
		cout << "start to do fa 2 jf" << endl;
		string tmpInputFaFile = inputDir_fa + tmpId + ".fa";
		string tmpOutputJfFile = outputDir_jf + tmpId + ".jf"; 
		cout << "tmpInputFaFile: " << tmpInputFaFile << endl;
		cout << "tmpOutputJfFile: " << tmpOutputJfFile << endl;
		Fa2jf_Info tmpFa2jfInfo;
		tmpFa2jfInfo.initaite_withMinCount(tmpInputFaFile, tmpOutputJfFile, jellyFishiBin,
			threads_num_str, count_min_str, Kmer_length_str, bf_size_str);
		tmpFa2jfInfo.fa2jf();
		cout << "start to do jf 2 Kmer" << endl;
		string tmpInputJfFile = outputDir_jf + tmpId + ".jf";
		string tmpOutputKmerFile = outputDir_Kmer + tmpId + ".Kmer";
		cout << "tmpInputJfFile: " << tmpInputJfFile << endl;
		cout << "tmpOutputKmerFile: " << tmpOutputKmerFile << endl;
		Jf2Kmer_Info jf2KmerInfo;
		jf2KmerInfo.initiate(tmpInputJfFile, tmpOutputKmerFile, jellyFishiBin);
		jf2KmerInfo.jf2Kmer();
		cout << "end of fa 2 jf 2 Kmer" << endl;
	}
	cout << "All jobs done!" << endl;
	return 0;
}