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
#include "../Ye_implementation/muloth.h"
#include <chrono>
#include <inttypes.h>
#include "../Ye_implementation/io_helper.h"
#include "../mps3Lib/read_block_test.h"
#include "../mps3Lib/otherFunc.h"
#include "../mps3Lib/index_info.h"
#include "../query/general/queryConstantDef.h"
#include "../query/general/querySeq_info.h"
#include "general/regionAssignment_info.h"

using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 8)
	{
		cout << "Executable inputGenomeIndexPath inputMothBitArray outputFolder windowSize KmerRepetitiveProportionMin class_num Kmer_length" << endl;
		exit(1);
	}
	string inputGenomeIndexPath = argv[1];
	string inputMothBitArray = argv[2];
	string outputFolder = argv[3];
	outputFolder += "/";
	string windowSizeStr = argv[4];
	string KmerRepetitiveProportionMinStr = argv[5];
	string class_num_str = argv[6];
	string Kmer_length_str = argv[7];
	string cmd_mkdir_output = "mkdir " + outputFolder;
	system(cmd_mkdir_output.c_str());
	string log_file = outputFolder + "log";
	ofstream log_ofs(log_file.c_str());
	string nonATGCgenomicRegion_file = outputFolder + "nonATGCgenomicRegion.txt";
	ofstream nonATGCgenomicRegion_ofs(nonATGCgenomicRegion_file.c_str());
	string repetitiveKmerRegion_file = outputFolder + "repetitiveKmerRegion.txt";
	ofstream repetitiveKmerRegion_ofs(repetitiveKmerRegion_file.c_str());
	string regionKmerClassInfo_file = outputFolder + "regionKmerClassInfo.txt";
	ofstream regionKmerClassInfo_ofs(regionKmerClassInfo_file.c_str());

    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "start to load index" << endl;
    MulOth<VALUELENGTH, keyT> * moth = new MulOth<VALUELENGTH, keyT>(argv[2]);
    
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "end of loading index" << endl;	
    cout << endl << "[" << asctime(local) << "start to do query" << endl;

	log_ofs << endl << "[" << asctime(local) << "start to initiate indexInfo" << endl;
	cout << endl << "[" << asctime(local) << "start to initiate indexInfo" << endl;
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
    nowtime = time(NULL);
    local = localtime(&nowtime);	
	cout << endl << "[" << asctime(local) << "end of initiating indexInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "end of initiating indexInfo" << endl;

	int windowSize = atoi(windowSizeStr.c_str());
	double KmerRepetitiveProportionMin = atof(KmerRepetitiveProportionMinStr.c_str());
	int class_num = atoi(class_num_str.c_str());
	int Kmer_length = atoi(Kmer_length_str.c_str());
	log_ofs << endl << "windowSize: " << windowSize << endl << "KmerRepetitiveProportionMin: " << KmerRepetitiveProportionMin 
		<< endl << "class_num: " << class_num << endl << "Kmer_length: " << Kmer_length << endl << endl;
	int jump_length = 1;
	int Kmer_num = windowSize - Kmer_length + 1;
	
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		string tmpChrName = indexInfo->returnChrNameStr(tmpChr);
		cout << "tmpChr: " << tmpChrName << endl;
		log_ofs << "tmpChr: " << tmpChrName << endl;
		int tmpChr_length = indexInfo->returnChromLength(tmpChr);
		cout << "tmpChr_length: " << tmpChr_length << endl;
		log_ofs << "tmpChr_length: " << tmpChr_length << endl;
		for(int tmpBase = 0; tmpBase < tmpChr_length; tmpBase += jump_length)
		{
			int tmpSeq_startPos = tmpBase + 1;
			int tmpSeq_endPos = tmpSeq_startPos + windowSize - 1;
			if(tmpSeq_endPos > tmpChr_length)
				break;
			string tmpChrSubSeq = indexInfo->returnChromStrSubstr(tmpChr, tmpSeq_startPos, windowSize);
			if((tmpChrSubSeq.find_first_not_of("ATGCatgc")) != string::npos)
			{
				nonATGCgenomicRegion_ofs << tmpChrName << "\t" << tmpSeq_startPos << "\t" << tmpSeq_startPos + windowSize - 1 << endl;
				continue;
			}			
			if(tmpChrSubSeq.find_first_of("atgc") != string::npos)
			{
				cout << "Error! a/t/g/c exists ..." << endl;
				exit(1);
			}

			string tmpChrSubSeq_id;
			QuerySeq_Info tmpQuerySeqInfo;
			tmpQuerySeqInfo.initiate(tmpChrSubSeq_id, tmpChrSubSeq);
			vector< pair<valueT,int> > tmpQuerySeq_classIdCountPairVec;
			tmpQuerySeqInfo.querySeq_returnClassIdCountPairVec(tmpQuerySeq_classIdCountPairVec, moth, Kmer_length, jump_length);
			
			RegionAssignment_Info tmpRegionAssignmentInfo;
			tmpRegionAssignmentInfo.initiate(class_num - 1, 0);
			tmpRegionAssignmentInfo.getMostKmerCountRegion_SE_1st_2nd(tmpQuerySeq_classIdCountPairVec);
			int tmpQuerySeq_bestClassId = tmpRegionAssignmentInfo.return_validDiscrimitiveClassIndex_best();
			int tmpQuerySeq_bestClassKmerCount = tmpRegionAssignmentInfo.return_validDiscrimitiveClassKmerCount_best();
			int tmpQuerySeq_repetitiveKmerCount = tmpRegionAssignmentInfo.return_validClassKmerCount_repetitive();

			double tmpQuerySeq_KmerRepetitiveProportion = (double)tmpQuerySeq_repetitiveKmerCount/(double)Kmer_num;
			if(tmpQuerySeq_KmerRepetitiveProportion >= KmerRepetitiveProportionMin)
				repetitiveKmerRegion_ofs << tmpChrName << "_" << tmpSeq_startPos << "_" << tmpSeq_endPos << endl;
			regionKmerClassInfo_ofs << tmpChrName << "_" << tmpSeq_startPos << "_" << tmpSeq_endPos 
				<< "\t" << (int)tmpQuerySeq_bestClassId << "_" << tmpQuerySeq_bestClassKmerCount << "_"
				<< tmpQuerySeq_repetitiveKmerCount << endl;
		}
	}

	delete indexInfo;
	free(chrom);
	nonATGCgenomicRegion_ofs.close();
	repetitiveKmerRegion_ofs.close();
	regionKmerClassInfo_ofs.close();
	chrom_bit_file_ifs.close();
	parameter_ifs.close();
	log_ofs.close();
	return 0;
}