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
//#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
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
#include "general/queryConstantDef.h"
#include "general/querySeq_info.h"

using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputLothelloNodeIndex inputFaFile outputFolder Kmer_length threads_num" << endl;
		exit(1);
	}
	int jump_length = 1;
	string outputDirStr = argv[3];
	outputDirStr += "/";
  	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());
   	string log_file = outputDirStr + "log.txt";
   	ofstream log_ofs(log_file.c_str());

    string kmer_length_str = argv[4];
    int kmer_length = atoi(kmer_length_str.c_str());
    string threads_num_str = argv[5];
    int threads_num = atoi(threads_num_str.c_str());

    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "start to load index" << endl;
    MulOth<VALUELENGTH, keyT> * moth = new MulOth<VALUELENGTH, keyT>(argv[1]);
    
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "end of loading index" << endl;	
    cout << endl << "[" << asctime(local) << "start to do query" << endl;

    string output_query_file = outputDirStr + "query.txt";
    //ofstream query_ofs(output_query_file.c_str());
    FILE* p_queryOutput = fopen(output_query_file.c_str(), "w");
    string InputReadFile = argv[2];
	ifstream inputRead_ifs(InputReadFile.c_str());
	
	int readTotalNum = 0;
	int normalRecordNum = 1000000;//0;
	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;
	string line1, line2;
	//int readPairNum = 0;	
	vector<string> readNameVec(normalRecordNum);
	vector<string> readSeqVec(normalRecordNum);
	vector<valueT> readBestClassIdVec(normalRecordNum);
	vector<int> readBestClassCountVec(normalRecordNum);
	vector<valueT> readSecondBestClassIdVec(normalRecordNum);
	vector<int> readSecondBestClassCountVec(normalRecordNum);	
	for(tmpTurn = 0; ; tmpTurn++)
	{
		if(EndOfRecord)
			break;		
		int recordNum = normalRecordNum;
		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
		realRecordNum = normalRecordNum;

		for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
		{
    		if(inputRead_ifs.eof())
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    			
    		}
    		getline(inputRead_ifs, line1);
    		if((inputRead_ifs.eof())||(line1 == ""))
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    			
    		}
    		readNameVec[recordNumTmp] = line1;
    		getline(inputRead_ifs, line2);	
    		readSeqVec[recordNumTmp] = line2;
		}
		readTotalNum += realRecordNum;
		nowtime = time(NULL);
		local = localtime(&nowtime);		
		cout << endl << "[" << asctime(local) << "finish reading Fasta/Fastq file, turn: " << tmpTurn + 1 << endl;
		cout << endl << "[" << asctime(local) << "start to fix reads in phase1, turn: " << tmpTurn + 1 << endl;
		
		omp_set_num_threads(threads_num);
		#pragma omp parallel for schedule(dynamic)
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{
			//cout << "tmpOpenMP: " << tmpOpenMP << endl;
			int threadNO = omp_get_thread_num();
			//cout << "start to initiate tmpQuerySeqInfo " << endl;
			QuerySeq_Info tmpQuerySeqInfo;
			tmpQuerySeqInfo.initiate(readNameVec[tmpOpenMP], readSeqVec[tmpOpenMP]);
			//cout << "start to do query" << endl;
			vector< pair<valueT,int> > tmpQuerySeq_classIdCountPairVec;
			tmpQuerySeqInfo.querySeq_returnClassIdCountPairVec(tmpQuerySeq_classIdCountPairVec, moth, kmer_length, jump_length);
			//cout << "start to select the best class" << endl;
			valueT tmpQuerySeq_bestClassId;
			int tmpQuerySeq_bestClassKmerCount;
			tmpQuerySeqInfo.select_best_class(tmpQuerySeq_bestClassId, tmpQuerySeq_bestClassKmerCount, tmpQuerySeq_classIdCountPairVec);
			//cout << "tmpQuerySeq_bestClassId: " << (int)tmpQuerySeq_bestClassId << endl;
			//cout << "tmpQuerySeq_bestClassKmerCount: " << tmpQuerySeq_bestClassKmerCount << endl;
			//cout << "start to select the second best class" << endl;
			valueT tmpQuerySeq_secondBestClassId;
			int tmpQuerySeq_secondBestClassKmerCount;
			tmpQuerySeqInfo.select_secondBest_class(tmpQuerySeq_secondBestClassId, tmpQuerySeq_secondBestClassKmerCount, tmpQuerySeq_bestClassId, tmpQuerySeq_classIdCountPairVec);
			//cout << "tmpQuerySeq_secondBestClassId: " << (int)tmpQuerySeq_secondBestClassId << endl;
			//cout << "tmpQuerySeq_secondBestClassKmerCount: " << tmpQuerySeq_secondBestClassKmerCount << endl;
			readBestClassIdVec[tmpOpenMP] = tmpQuerySeq_bestClassId;
			readBestClassCountVec[tmpOpenMP] = tmpQuerySeq_bestClassKmerCount;
			readSecondBestClassIdVec[tmpOpenMP] = tmpQuerySeq_secondBestClassId;
			readSecondBestClassCountVec[tmpOpenMP] = tmpQuerySeq_secondBestClassKmerCount;			
		}
		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		cout << endl << "[" << asctime(local) << "start to output ... turn: " << tmpTurn+1 << endl;		
		// for(int tmp = 0; tmp < realRecordNum; tmp++)
		// 	query_ofs << readNameVec[tmp] << "\t" << (int)readBestClassIdVec[tmp] << "\t" 
		// 		<< readBestClassCountVec[tmp] << "\t" << (int)readSecondBestClassIdVec[tmp] << "\t" 
		// 		<< readSecondBestClassCountVec[tmp] << endl << readSeqVec[tmp] << endl;
		for(int tmp = 0; tmp < realRecordNum; tmp++)
			fprintf(p_queryOutput, "%s\t%d\t%d\t%d\t%d\n%s\n", readNameVec[tmp].c_str(), (int)readBestClassIdVec[tmp],
				readBestClassCountVec[tmp], (int)readSecondBestClassIdVec[tmp], readSecondBestClassCountVec[tmp], readSeqVec[tmp].c_str());
	}
	//query_ofs.close();
	fclose(p_queryOutput);
    nowtime = time(NULL);
    local = localtime(&nowtime);
    cout << endl << "[" << asctime(local) << "end of doing query" << endl;
    log_ofs.close();
    delete moth;
	return 0;
}