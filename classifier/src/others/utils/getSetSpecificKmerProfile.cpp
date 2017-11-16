// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
// MEMORY: 1.5*N*M/256; metagenomics: 9*M/256
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <errno.h>
#include "../Ye_implementation/othello.h"
#include "../Ye_implementation/mulothindex.h"
#include <chrono>
#include <inttypes.h>
#include "../Ye_implementation/io_helper.h"
#include "../mps3Lib/read_block_test.h"
#include "../mps3Lib/otherFunc.h"
#include "../mps3Lib/index_info.h"
//#include "../query/general/queryConstantDef.h"
//#include "../query/general/querySeq_info.h"
//#include "../../Ye_implementation/othelloindex.h"
//#include "general/regionAssignment_info.h"

using namespace std;

typedef unsigned long long keyT;
typedef uint64_t valueT;
typedef uint16_t setIdT;
typedef uint8_t freqT;

vector<keyT> keys;
vector<valueT> values;

keyT convertKmer2Key(char*s, int kmerLength)
{
    int tmpLength = 0;
    switch (*s) 
    {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            keyT ret = 0;
            while (*s == 'A' || *s == 'C' || *s =='T' || *s =='G')
            {
                ret <<=2;
                switch (*s) 
                {
                    case 'T':
                        ret ++;
                    case 'G':
                        ret ++;
                    case 'C':
                        ret ++;
                }
                s++;
                tmpLength ++;
                if(tmpLength == kmerLength)
                    return ret;
            }
    }
}

int main(int argc, char** argv) 
{
	if(argc != 8)
	{
		cout << "Executable#0 inputKmerIndexBitArrayFile#1" << endl;
		cout << "inputKmerFreqUnionFile#2 inputKmerSetFileListFile#3" << endl;
		cout << "outputDir#4 KmerLength#5 splitBit#6 KmerNumMax#7" << endl;
		exit(1);
	}
	string inputKmerSetFileListFile = argv[3];
	string outputDir = argv[4];
	outputDir += "/";
	string cmd_mkdir = "mkdir " + outputDir;
	system(cmd_mkdir.c_str());
	string output_stats = outputDir + "stats.txt";
	string output_log = outputDir + "log.txt";
	ofstream stats_ofs(output_stats.c_str());
	ofstream log_ofs(output_log.c_str());

	vector<string> allKmerFilePathVec;
	ifstream KmerSetFileList_ifs(inputKmerSetFileListFile.c_str());
	while(!KmerSetFileList_ifs.eof())
	{
		string tmpStr;
		getline(KmerSetFileList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		allKmerFilePathVec.push_back(tmpStr);
		log_ofs << "tmpFile: " << tmpStr << endl;
	}
	KmerSetFileList_ifs.close();
	int setNum = allKmerFilePathVec.size();
	int repetitiveKmerSetId = setNum + 1;

	string KmerLengthStr = argv[5];
	int Kmer_length = atoi(KmerLengthStr.c_str());
    string splitbitStr = argv[6];
    int splitbit = atoi(splitbitStr.c_str());
    string KmerNumMaxStr = argv[7];
    unsigned long long int KmerNumMax = atoll(KmerNumMaxStr.c_str());

    IOHelper<keyT,valueT> *helper;
    helper = new ConstantLengthKmerHelper<keyT,valueT>(Kmer_length,splitbit);

    cout << "start to load union kmer index bit array file" << endl;
    MulOthIndex<keyT> * moth;
    moth = new MulOthIndex<keyT>(argv[1], helper);

    cout << "start to load union kmer frequency file and output repetitive Kmer" << endl;
    freqT* freqArray = new freqT[KmerNumMax];
    string inputKmerFreqUnionFile = argv[2];
    string output_repetitive_Kmer_file = outputDir + "repetitive.Kmer";
    ofstream repetitiveKmer_ofs(output_repetitive_Kmer_file.c_str());
    string output_alien_Kmer_file = outputDir + "alien.Kmer";
    ofstream alienKmer_ofs(output_alien_Kmer_file.c_str());    
    ifstream KmerFreqUnion_ifs(inputKmerFreqUnionFile.c_str());
    unsigned long long tmpLineNO = 0; 
    while(!KmerFreqUnion_ifs.eof())
    {
    	string tmpStr;
    	getline(KmerFreqUnion_ifs, tmpStr);
    	if(tmpStr == "")
    		break;
        tmpLineNO ++;
        unsigned long long tmpThousandIndex = tmpLineNO / 1000000;
        if(tmpLineNO == tmpThousandIndex * 1000000)          
            cout << "Processed Line #: " << tmpLineNO << endl; 
    	string tmpKmerStr = tmpStr.substr(0, Kmer_length);
	    char tmpKmerCharArray[Kmer_length];
		std::copy(tmpKmerStr.begin(), tmpKmerStr.end(), tmpKmerCharArray);
		keyT tmpKeyT = convertKmer2Key(tmpKmerCharArray, Kmer_length);
		valueT tmpValueT = moth->query(tmpKeyT);
    	string tmpFreqStr = tmpStr.substr(Kmer_length + 1);
    	unsigned long long tmpFreq = atoll(tmpFreqStr.c_str());
        if(tmpFreq <= 255)
 		   	freqArray[tmpValueT] = tmpFreq;
    	else
    		freqArray[tmpValueT] = 255;
    	if(tmpFreq >= 2)
            repetitiveKmer_ofs << tmpKmerStr << "\t" << repetitiveKmerSetId << endl;
        else if(tmpFreq == 0)
            alienKmer_ofs << tmpKmerStr << "\t0" << endl;
        else // set-specific Kmer
        {}
    }
    KmerFreqUnion_ifs.close();
    repetitiveKmer_ofs.close();
    alienKmer_ofs.close();

    cout << "start to initiate setIdArray" << endl;
   	setIdT* setIdArray = new setIdT[KmerNumMax];

   	cout << "start to get set-specific Kmer file for each set" << endl;
   	string outputDir_setSpecific = outputDir + "set_specific/";
   	string mkdir_setSpecific = "mkdir " + outputDir_setSpecific;
   	system(mkdir_setSpecific.c_str());
    int kmerSetNum = allKmerFilePathVec.size();
    for(int tmpSetIndex = 0; tmpSetIndex < kmerSetNum; tmpSetIndex++)
    {
    	string output_tmpSetSpecific = outputDir_setSpecific + int_to_str(tmpSetIndex+1) + ".Kmer";
    	ofstream tmpSetSpecific_ofs(output_tmpSetSpecific.c_str());
    	string tmpKmerSetFile = allKmerFilePathVec[tmpSetIndex];
    	ifstream tmpKmerSet_ifs(tmpKmerSetFile.c_str());
		while(!tmpKmerSet_ifs.eof())
		{
			string tmpStr;
			getline(tmpKmerSet_ifs, tmpStr);
			if(tmpStr == "")
				break;
			string tmpKmerStr = tmpStr.substr(0, Kmer_length);
			char tmpSeqChar[Kmer_length];
			std::copy(tmpKmerStr.begin(), tmpKmerStr.end(), tmpSeqChar);
			keyT tmpKeyT = convertKmer2Key(tmpSeqChar, Kmer_length);
			valueT tmpValueT = moth->query(tmpKeyT);
			if((freqArray[tmpValueT] == 1)&&(setIdArray[tmpValueT] == 0))
			{
				int tmpSetId = tmpSetIndex + 1;
				tmpSetSpecific_ofs << tmpKmerStr << "\t" << tmpSetId << endl;
				setIdArray[tmpValueT] = tmpSetId;
			}
		}
    	tmpKmerSet_ifs.close();
    	tmpSetSpecific_ofs.close();
    }
	cout << "All jobs done!" << endl;
    

    delete freqArray;
    delete helper;
    delete setIdArray;
	delete moth;
	stats_ofs.close();
	log_ofs.close();
	return 0;
}