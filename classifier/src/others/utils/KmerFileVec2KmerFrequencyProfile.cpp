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

using namespace std;

typedef unsigned long long keyT;
typedef uint64_t valueT;
typedef uint8_t freqT;

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
		    cout << "Executable#0 inputUnionKmerIndex#1 inputUnionKmerTwoColumnFile#2" << endl;
		    cout << "inputKmerSetTwoColumnFileListFile#3" << endl;
		    cout << "outputDir#4 KmerLength#5 splitBit#6 KmerNumMax#7" << endl;
		    exit(1);
	  }
	  cout << "start to creat dir" << endl;
	  string outputDir = argv[4];
	  outputDir += "/";
	  string cmd_mkdir = "mkdir " + outputDir;
	  system(cmd_mkdir.c_str());
	  string output_stats = outputDir + "stats.txt";
	  string output_log = outputDir + "log.txt";
	  ofstream stats_ofs(output_stats.c_str());
	  ofstream log_ofs(output_log.c_str());	

	  cout << "start to creat KmerSetFileListFileVec" << endl;
	  string inputKmerSetFileListFile = argv[3];
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

	  string KmerLengthStr = argv[5];
	  int Kmer_length = atoi(KmerLengthStr.c_str());
    string splitbitStr = argv[6];
    int splitbit = atoi(splitbitStr.c_str());
    
    cout << "start to initiate IOhelper" << endl;
    IOHelper<keyT,valueT> *helper;
    helper = new ConstantLengthKmerHelper<keyT,valueT>(Kmer_length, splitbit);
    
    cout << "start to load union kmer with count file" << endl;
    MulOthIndex<keyT> * moth;
    moth = new MulOthIndex<keyT>(argv[1], helper);
    
    cout << "start to initiate freqArray" << endl;
    string KmerNumMaxStr = argv[7];
   	unsigned long long int KmerNumMax = atoll(KmerNumMaxStr.c_str());
   	freqT* freqArray = new freqT[KmerNumMax];
    
    cout << "start to load each Kmer set file" << endl;
    int kmerSetNum = allKmerFilePathVec.size();
    for(int tmp = 0; tmp < kmerSetNum; tmp++)
    {
    	  string tmpKmerSetFile = allKmerFilePathVec[tmp];
        cout << "tmpKmerSetFile: " << tmp << endl << tmpKmerSetFile << endl;
    	  ifstream tmpKmerSet_ifs(tmpKmerSetFile.c_str());
		    while(!tmpKmerSet_ifs.eof())
		    {
			      string tmpStr;
		      	getline(tmpKmerSet_ifs, tmpStr);
			      if(tmpStr == "")
				        break;
            string tmpKmerStr = tmpStr.substr(0, Kmer_length);
			      char tmpKmerCharArray[Kmer_length];
			      std::copy(tmpKmerStr.begin(), tmpKmerStr.end(), tmpKmerCharArray);
			      keyT tmpKeyT = convertKmer2Key(tmpKmerCharArray, Kmer_length);
			      valueT tmpValueT = moth->query(tmpKeyT);
			      freqT tmpOriFreq = freqArray[tmpValueT];
			      if(tmpOriFreq == UINT8_MAX) // new
			          {}
			      else
				        freqArray[tmpValueT] = tmpOriFreq + 1; // originally set-specific, convert 2 repetitive 
		    }
    	  tmpKmerSet_ifs.close();
    }
   	
    cout << "start to scan whole freqArray and fill specificFreqKmerCountArray and output Kmer freq profile" << endl;
    vector<ofstream*> specificFreqKmerOfsVec;
    for(uint8_t tmpFreq = 0; tmpFreq <= UINT8_MAX; tmpFreq++)
    {
      	string tmpOutputFile = outputDir + "Kmer.freq." + int_to_str((int)tmpFreq);
    	  ofstream* tmp_ofs = new ofstream(tmpOutputFile.c_str());
      	specificFreqKmerOfsVec.push_back(tmp_ofs);
    }

   	uint64_t* specificFreqKmerCountArray = new uint64_t[UINT8_MAX + 1]; 
   	string outputKmerFrequencyProfile = outputDir + "Kmer.freq";;
   	ofstream KmerFreq_ofs(outputKmerFrequencyProfile.c_str());
   	string inputUnionKmerTwoColumnFile = argv[2];
   	ifstream unionKmerTwoColumn_ifs(inputUnionKmerTwoColumnFile.c_str());
    unsigned long long tmpLineNO = 0;
   	while(!unionKmerTwoColumn_ifs.eof())
   	{
	      string tmpStr;
	   	  getline(unionKmerTwoColumn_ifs, tmpStr);
	   	  if(tmpStr == "")
	   		    break;
        tmpLineNO ++;
        int tmpThousandIndex = tmpLineNO / 100000;
        if(tmpLineNO == tmpThousandIndex * 100000)          
            cout << "Processed Line #: " << tmpLineNO << endl;
        string tmpKmerStr = tmpStr.substr(0, Kmer_length);
	      char tmpKmerCharArray[Kmer_length];
	      std::copy(tmpKmerStr.begin(), tmpKmerStr.end(), tmpKmerCharArray);
		    keyT tmpKeyT = convertKmer2Key(tmpKmerCharArray, Kmer_length);
		    valueT tmpValueT = moth->query(tmpKeyT);   		
   		  freqT tmpFreq = freqArray[tmpValueT];
   		  specificFreqKmerCountArray[tmpFreq] ++;
   		  KmerFreq_ofs << tmpKmerStr << "\t" << tmpFreq << endl;
   		  (*specificFreqKmerOfsVec[tmpFreq]) << tmpKmerStr << endl;
   	}

   	cout << "start to output specificFreqKmerCountArray" << endl;
   	string outputSpecificFreqKmerCountFile = outputDir + "Kmer.freq.count";
   	ofstream specificFreqKmerCount_ofs(outputSpecificFreqKmerCountFile.c_str());
   	for(uint8_t tmp = 0; tmp <= UINT8_MAX; tmp++)
   		  specificFreqKmerCount_ofs << "Freq:\t" << tmp << "\tCount:" << specificFreqKmerCountArray[tmp] << endl;

    for(uint8_t tmpFreq = 0; tmpFreq <= UINT8_MAX; tmpFreq++)
    {
      	(*specificFreqKmerOfsVec[tmpFreq]).close();
      	delete specificFreqKmerOfsVec[tmpFreq];
    }   	

    cout << "All jobs done!" << endl;
   	specificFreqKmerCount_ofs.close();
   	unionKmerTwoColumn_ifs.close();
  	delete specificFreqKmerCountArray;
  	delete helper;
  	delete moth;
  	delete freqArray;
  	KmerFreq_ofs.close();
  	stats_ofs.close();
  	log_ofs.close();
  	return 0;
}