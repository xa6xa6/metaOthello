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
	if((argc != 10)&&(argc != 11))
	{
		cout << "Executable inputLothelloNodeIndex#1 outputFolder#2 Kmer_length#3 threads_num#4 " << endl;
		cout << "class_num#5 fa_or_fq#6 SE_or_PE#7 DiscrimitiveKmerProportionMin#8 inputReadFile_1#9 (inputReadFile_2#10)" << endl;
		exit(1);
	}
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "start to load index" << endl;
    string kmer_length_str = argv[3];
    int kmer_length = atoi(kmer_length_str.c_str());
    string threads_num_str = argv[4];
    int threads_num = atoi(threads_num_str.c_str());
    string class_num_str = argv[5];
    int class_num = atoi(class_num_str.c_str());
    int repetitiveClassId = class_num - 1;
    //int kmer_length = 20;
    int splitbit = 3;
    IOHelper<keyT,valueT> *helper;
    helper = new ConstantLengthKmerHelper<keyT,valueT>(kmer_length,splitbit);
    MulOth<keyT> * moth;
    moth = new MulOth<keyT>(argv[1], helper);
    
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "end of loading index" << endl;	

	int jump_length = 1;
	string DiscrimitiveKmerProportionMinStr = argv[8];
	double DiscrimitiveKmerProportionMin = atof(DiscrimitiveKmerProportionMinStr.c_str());
    string fa_or_fq_str = argv[6];
    bool fa_or_fq_bool;
    if((fa_or_fq_str == "Fa")||(fa_or_fq_str == "FA")||(fa_or_fq_str == "fa")
    	||(fa_or_fq_str == "Fasta")||(fa_or_fq_str == "FASTA")||(fa_or_fq_str == "fasta"))
    	fa_or_fq_bool = true;
    else if((fa_or_fq_str == "Fq")||(fa_or_fq_str == "FQ")||(fa_or_fq_str == "fq")
    	||(fa_or_fq_str == "Fastq")||(fa_or_fq_str == "FASTQ")||(fa_or_fq_str == "fastq"))
    	fa_or_fq_bool = false;
    else
    {
    	cout << "invalid parameter for fa_or_fq: " << fa_or_fq_str << endl << "Should be Fa or Fq;" << endl;
    	exit(1);
    }
    bool InputAsFastq = (!fa_or_fq_bool);

    string SE_or_PE_str = argv[7];
    bool SE_or_PE_bool;
    if((SE_or_PE_str == "SE")||(SE_or_PE_str == "se")||(SE_or_PE_str == "Se"))
    	SE_or_PE_bool = true;
    else if((SE_or_PE_str == "PE")||(SE_or_PE_str == "pe")||(SE_or_PE_str == "Pe"))
    	SE_or_PE_bool = false;
    else
    {
    	cout << "invalid parameter for SE_or_PE_str: " << SE_or_PE_str << endl << "Shoule be SE or PE" << endl;
    	exit(1);
    }

    if((SE_or_PE_bool && (argc != 10)) || ((!SE_or_PE_bool) && (argc != 11)))
    {
    	cout << "inconsistent SE_or_PE_bool and argc" << endl;
    	cout << "SE_or_PE_bool: " << SE_or_PE_bool << endl;
    	cout << "argc: " << argc << endl;
    	exit(1);
    }

	string outputDirStr = argv[2];
	outputDirStr += "/";
  	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());
   	string log_file = outputDirStr + "log.txt";
   	ofstream log_ofs(log_file.c_str());

	//vector<ofstream*> regionAssignedFaOfsVec;
	vector<FILE*> regionAssignedReadFpVec_1;
	vector<FILE*> regionAssignedReadFpVec_2;
	string faOrFq_suffix_str = ".fq";
	if(fa_or_fq_bool)
		faOrFq_suffix_str = ".fa";
    string invalidAssignment_1_file = outputDirStr + "invalid.1" + faOrFq_suffix_str;
    FILE* fp_invalidAssignment_1 = fopen(invalidAssignment_1_file.c_str(), "w");
	string invalidAssignment_2_file = outputDirStr + "invalid.2" + faOrFq_suffix_str;
	FILE* fp_invalidAssignment_2 = fopen(invalidAssignment_2_file.c_str(), "w");
	// end_1 file
    for(int tmp = 1; tmp <= class_num - 2; tmp++)
    {
    	string tmpRegionAssignedReadFile = outputDirStr + "region_" + int_to_str(tmp) + ".1";
    	tmpRegionAssignedReadFile += faOrFq_suffix_str;
    	FILE* fp_tmpReigonAssignedRead = fopen(tmpRegionAssignedReadFile.c_str(), "w");
    	regionAssignedReadFpVec_1.push_back(fp_tmpReigonAssignedRead);
    }
    if(!SE_or_PE_bool) // end 2 file, if SE_or_PE_bool = false;
    {
	    for(int tmp = 1; tmp <= class_num - 2; tmp++)
	    {
	    	string tmpRegionAssignedReadFile = outputDirStr + "region_" + int_to_str(tmp) + ".2";
	    	tmpRegionAssignedReadFile += faOrFq_suffix_str;
	    	FILE* fp_tmpReigonAssignedRead = fopen(tmpRegionAssignedReadFile.c_str(), "w");
	    	regionAssignedReadFpVec_2.push_back(fp_tmpReigonAssignedRead);
	    }
    }
    cout << endl << "[" << asctime(local) << "start to do query" << endl;
    //string output_query_file = outputDirStr + "query.txt";
    //ofstream query_ofs(output_query_file.c_str());
    string InputReadFile = argv[9];
	ifstream inputRead_ifs(InputReadFile.c_str());
	string InputReadFile_PE;
	if(!SE_or_PE_bool)
		InputReadFile_PE = argv[10];
	else
		InputReadFile_PE = "NULL";
	ifstream inputRead_PE_ifs(InputReadFile_PE.c_str());

	int readTotalNum = 0;
	int normalRecordNum = 10000000;//0;
	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;
	string line1, line2, line3, line4, line1_PE, line2_PE, line3_PE, line4_PE;
	//int readPairNum = 0;	
	// *************** needed for both SE and PE *************//
	vector<string> readName1Vec(normalRecordNum);
	vector<string> readSeq1Vec(normalRecordNum);
	vector<string> readQualSeq1Vec(normalRecordNum);
	// *******************************************************//
	// ***************  only needed for PE  ******************//
	vector<string> readName2Vec(normalRecordNum);
	vector<string> readSeq2Vec(normalRecordNum);
	vector<string> readQualSeq2Vec(normalRecordNum);

	vector<int> readAssignedRegionIndexVec(normalRecordNum);
	
	#ifdef ASSIGN_INFO
	string output_assign_info = outputDirStr + "assignInfo.txt";
	FILE* fp_assignInfo = fopen(output_assign_info.c_str(), "w");
	string output_readAssignment = outputDirStr + "readAssignment.txt";
	FILE* fp_readAssignment = fopen(output_readAssignment.c_str(), "w");
	vector<double> readAssignedRegionConfidenceScoreVec(normalRecordNum);
	vector<double> readAssignedRegionUniquenessScoreVec(normalRecordNum);	
	vector<int> readAssignmentClassIndex_best(normalRecordNum);
	vector<int> readAssignmentClassIndex_secondBest(normalRecordNum);
	vector<int> readAssignmentClassKmerCount_best(normalRecordNum);
	vector<int> readAssignmentClassKmerCount_secondBest(normalRecordNum);
	vector<int> readAssignmentClassKmerCount_repetitive(normalRecordNum);
	vector<int> readAssignmentClassKmerCount_invalid(normalRecordNum);
	vector<int> readAssignmentClassKmerCount_totalValid(normalRecordNum);
	#endif
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
    		readName1Vec[recordNumTmp] = line1;
    		getline(inputRead_ifs, line2);	
    		readSeq1Vec[recordNumTmp] = line2;
    		// read quality sequences for end 1 reads
    		if(InputAsFastq)
    		{
    			getline(inputRead_ifs, line3);
    			getline(inputRead_ifs, line4);
    			readQualSeq1Vec[recordNumTmp] = line4;   
    		}
    		if(!SE_or_PE_bool) // if paired end reads are provied
    		{	
	    		getline(inputRead_PE_ifs, line1_PE); // readName_2
	    		readName2Vec[recordNumTmp] = line1_PE;
	    		getline(inputRead_PE_ifs, line2_PE); // readSeq_2
	    		readSeq2Vec[recordNumTmp] = line2_PE;
	    		if(InputAsFastq)
	    		{
	    			getline(inputRead_PE_ifs, line3_PE);
	    			getline(inputRead_PE_ifs, line4_PE);
	    			readQualSeq2Vec[recordNumTmp] = line4_PE;
	    		}
    		}
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
			QuerySeq_Info tmpQuerySeqInfo_1;
			tmpQuerySeqInfo_1.initiate(readName1Vec[tmpOpenMP], readSeq1Vec[tmpOpenMP]);
			vector< pair<valueT,int> > tmpQuerySeq_classIdCountPairVec_1;
			tmpQuerySeqInfo_1.querySeq_returnClassIdCountPairVec(
				tmpQuerySeq_classIdCountPairVec_1, moth, kmer_length, jump_length);
			
			RegionAssignment_Info tmpRegionAssignmentInfo;
			tmpRegionAssignmentInfo.initiate(class_num-1, 0);
			if(SE_or_PE_bool)
				tmpRegionAssignmentInfo.getMostKmerCountRegion_SE_1st_2nd(
					tmpQuerySeq_classIdCountPairVec_1);
			else
			{
				QuerySeq_Info tmpQuerySeqInfo_2;
				tmpQuerySeqInfo_2.initiate(readName2Vec[tmpOpenMP], readSeq2Vec[tmpOpenMP]);
				vector< pair<valueT,int> > tmpQuerySeq_classIdCountPairVec_2;
				tmpQuerySeqInfo_2.querySeq_returnClassIdCountPairVec(
					tmpQuerySeq_classIdCountPairVec_2, moth, kmer_length, jump_length);
				tmpRegionAssignmentInfo.getMostKmerCountRegion_PE_1st_2nd(
					tmpQuerySeq_classIdCountPairVec_1, tmpQuerySeq_classIdCountPairVec_2);
			}

			double tmpDiscrimitiveKmerProportion = (double)(tmpRegionAssignmentInfo.return_validDiscrimitiveClassKmerCount_best())/(double)(tmpRegionAssignmentInfo.return_validDiscrimitiveClassKmerCount_total());
			if((tmpDiscrimitiveKmerProportion >= DiscrimitiveKmerProportionMin)&&(tmpRegionAssignmentInfo.return_validDiscrimitiveClassKmerCount_best() > 2))
				readAssignedRegionIndexVec[tmpOpenMP] = tmpRegionAssignmentInfo.return_validDiscrimitiveClassIndex_best();
			else
				readAssignedRegionIndexVec[tmpOpenMP] = 0;
			// if(tmpRegionAssignmentInfo.return_validDiscrimitiveClassKmerCount_best() 
			// 	== tmpRegionAssignmentInfo.return_validDiscrimitiveClassKmerCount_secondBest())
			// 	readAssignedRegionIndexVec[tmpOpenMP] = 0;
			#ifdef ASSIGN_INFO
			tmpRegionAssignmentInfo.get_confidence_score_best();
			tmpRegionAssignmentInfo.get_uniqueness_score_best();
			readAssignedRegionConfidenceScoreVec[tmpOpenMP] 
				= tmpRegionAssignmentInfo.return_confidence_score_best();
			readAssignedRegionUniquenessScoreVec[tmpOpenMP] 
				= tmpRegionAssignmentInfo.return_uniqueness_score_best();			
			readAssignmentClassIndex_best[tmpOpenMP]
				= tmpRegionAssignmentInfo.return_validDiscrimitiveClassIndex_best();
			readAssignmentClassIndex_secondBest[tmpOpenMP]
				= tmpRegionAssignmentInfo.return_validDiscrimitiveClassIndex_secondBest();
			readAssignmentClassKmerCount_best[tmpOpenMP]
				= tmpRegionAssignmentInfo.return_validDiscrimitiveClassKmerCount_best();
			readAssignmentClassKmerCount_secondBest[tmpOpenMP]
				= tmpRegionAssignmentInfo.return_validDiscrimitiveClassKmerCount_secondBest();
			readAssignmentClassKmerCount_repetitive[tmpOpenMP]
				= tmpRegionAssignmentInfo.return_validClassKmerCount_repetitive();
			readAssignmentClassKmerCount_invalid[tmpOpenMP]
				= tmpRegionAssignmentInfo.return_invalidClassKmerCount();
			readAssignmentClassKmerCount_totalValid[tmpOpenMP]
				= tmpRegionAssignmentInfo.return_validDiscrimitiveClassKmerCount_total();
			#endif
		}
		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		cout << endl << "[" << asctime(local) << "start to output ... turn: " << tmpTurn+1 << endl;		
		if(fa_or_fq_bool)
		{			
			for(int tmp = 0; tmp < realRecordNum; tmp++)
			{
				int tmpRegionIndex = readAssignedRegionIndexVec[tmp];
				if((tmpRegionIndex < 1)||(tmpRegionIndex > class_num - 2))
					fprintf(fp_invalidAssignment_1, "%s\n%s\n", 
						readName1Vec[tmp].c_str(), readSeq1Vec[tmp].c_str());
				else	
					fprintf(regionAssignedReadFpVec_1[tmpRegionIndex - 1], "%s\n%s\n", 
						readName1Vec[tmp].c_str(), readSeq1Vec[tmp].c_str());
			}
			if(!SE_or_PE_bool) // PE, fa
			{
				for(int tmp = 0; tmp < realRecordNum; tmp++)
				{
					int tmpRegionIndex = readAssignedRegionIndexVec[tmp];
					if((tmpRegionIndex < 1)||(tmpRegionIndex > class_num - 2))
						fprintf(fp_invalidAssignment_2, "%s\n%s\n", 
							readName2Vec[tmp].c_str(), readSeq2Vec[tmp].c_str());
					else
						fprintf(regionAssignedReadFpVec_2[tmpRegionIndex - 1], "%s\n%s\n", 
							readName2Vec[tmp].c_str(), readSeq2Vec[tmp].c_str());				
				}
			}
		}
		else // fq
		{	
			for(int tmp = 0; tmp < realRecordNum; tmp++)
			{
				int tmpRegionIndex = readAssignedRegionIndexVec[tmp];
				if((tmpRegionIndex < 1)||(tmpRegionIndex > class_num - 2))
					fprintf(fp_invalidAssignment_1, "%s\n%s\n+%s\n%s\n", 
						readName1Vec[tmp].c_str(), readSeq1Vec[tmp].c_str(), 
						readName1Vec[tmp].c_str(), readQualSeq1Vec[tmp].c_str());
				else	
					fprintf(regionAssignedReadFpVec_1[tmpRegionIndex - 1], "%s\n%s\n+%s\n%s\n", 
						readName1Vec[tmp].c_str(), readSeq1Vec[tmp].c_str(), 
						readName1Vec[tmp].c_str(), readQualSeq1Vec[tmp].c_str());
			}
			if(!SE_or_PE_bool) // PE, fq
			{
				for(int tmp = 0; tmp < realRecordNum; tmp++)
				{
					int tmpRegionIndex = readAssignedRegionIndexVec[tmp];
					if((tmpRegionIndex < 1)||(tmpRegionIndex > class_num - 2))
						fprintf(fp_invalidAssignment_2, "%s\n%s\n+%s\n%s\n", 
							readName2Vec[tmp].c_str(), readSeq2Vec[tmp].c_str(),
							readName2Vec[tmp].c_str(), readQualSeq2Vec[tmp].c_str());
					else
						fprintf(regionAssignedReadFpVec_2[tmpRegionIndex - 1], "%s\n%s\n+%s\n%s\n", 
							readName2Vec[tmp].c_str(), readSeq2Vec[tmp].c_str(),
							readName2Vec[tmp].c_str(), readQualSeq2Vec[tmp].c_str());					
				}
			}			
		}
		#ifdef ASSIGN_INFO
		for(int tmp = 0; tmp < realRecordNum; tmp++)
		{	
			fprintf(fp_assignInfo, "%s\t%s\n%f\t%f\n%d\t%d\n%d\t%d\n%d\t%d\t%d\n", 
				readName1Vec[tmp].c_str(), readName2Vec[tmp].c_str(), 
				readAssignedRegionConfidenceScoreVec[tmp], readAssignedRegionUniquenessScoreVec[tmp],
				readAssignmentClassIndex_best[tmp], readAssignmentClassKmerCount_best[tmp],
				readAssignmentClassIndex_secondBest[tmp], readAssignmentClassKmerCount_secondBest[tmp],
				readAssignmentClassKmerCount_totalValid[tmp],
				readAssignmentClassKmerCount_repetitive[tmp], readAssignmentClassKmerCount_invalid[tmp]);		
			fprintf(fp_readAssignment, "%s\t%d\n", readName1Vec[tmp].c_str(), readAssignedRegionIndexVec[tmp]);
		}
		#endif
	}

	#ifdef ASSIGN_INFO
	fclose(fp_assignInfo);
	fclose(fp_readAssignment);
	#endif
    for(int tmp = 1; tmp <= class_num - 2; tmp++)
    	fclose(regionAssignedReadFpVec_1[tmp - 1]);
    if(!SE_or_PE_bool) // end 2 file, if SE_or_PE_bool = false;
    {
	    for(int tmp = 1; tmp <= class_num - 2; tmp++)
	    	fclose(regionAssignedReadFpVec_2[tmp - 1]); 	
    }
    fclose(fp_invalidAssignment_1);
    fclose(fp_invalidAssignment_2);
	inputRead_ifs.close();
	inputRead_PE_ifs.close();
    nowtime = time(NULL);
    local = localtime(&nowtime);    
    cout << endl << "[" << asctime(local) << "end of doing query" << endl;
    log_ofs.close();
    delete helper;
    delete moth;
	return 0;
}