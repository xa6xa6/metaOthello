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
#include "../../Ye_implementation/othello.h"
#include "../../Ye_implementation/muloth.h"
#include <chrono>
#include <inttypes.h>
#include "../../Ye_implementation/io_helper.h"
#include "../../mps3Lib/read_block_test.h"
#include "../../mps3Lib/otherFunc.h"
#include "../../mps3Lib/index_info.h"
#include "../../query/general/queryConstantDef.h"
#include "../../query/general/querySeq_info.h"
#include "../general/geneCount_info.h"
#include "../general/geneAssignment_info.h"
#include "../general/fusionCount_info.h"
#include "../general/fusionAssignment_info.h"
using namespace std;

time_t nowtime;
struct tm *local;

void get_Fa_or_Fq_SE_or_PE_bool(bool& fa_or_fq_bool, bool& SE_or_PE_bool, string& fa_or_fq_str, string& SE_or_PE_str, int argc)
{
    if((fa_or_fq_str == "Fa")||(fa_or_fq_str == "FA")||(fa_or_fq_str == "fa")||(fa_or_fq_str == "Fasta")||(fa_or_fq_str == "FASTA")||(fa_or_fq_str == "fasta"))
    	fa_or_fq_bool = true;
    else if((fa_or_fq_str == "Fq")||(fa_or_fq_str == "FQ")||(fa_or_fq_str == "fq")||(fa_or_fq_str == "Fastq")||(fa_or_fq_str == "FASTQ")||(fa_or_fq_str == "fastq"))
    	fa_or_fq_bool = false;
    else
    {
    	cout << "invalid parameter for fa_or_fq: " << fa_or_fq_str << endl << "Should be Fa or Fq;" << endl;
    	exit(1);
    }
    // bool InputAsFastq = (!fa_or_fq_bool);
    // bool SE_or_PE_bool;
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
}

void reorder_removeIsolatedKmerValue(vector<valueT>& tmpKmerValueVec_filtered, vector<valueT>& tmpKmerValueVec_raw, int tmpQuerySeqKmerNum_1)
{
	int tmpKmerValueVec_raw_size = tmpKmerValueVec_raw.size();
	if(tmpKmerValueVec_raw_size == 1)
		tmpKmerValueVec_filtered.push_back(0);
	else if(tmpKmerValueVec_raw_size == 2)
	{
		if(tmpKmerValueVec_raw[0] == tmpKmerValueVec_raw[1])
			tmpKmerValueVec_filtered.push_back(tmpKmerValueVec_raw[0]);
		else
			tmpKmerValueVec_filtered.push_back(0);
	}
	else // tmpKmerValueVec_raw_size > 2
	{
		int tmpQuerySeqKmerNum_2 = tmpKmerValueVec_raw_size - tmpQuerySeqKmerNum_1;
		// the 1st set of KmerValueVec
		if(tmpQuerySeqKmerNum_1 == 1)
			tmpKmerValueVec_filtered.push_back(0);
		else if(tmpQuerySeqKmerNum_1 == 2)
		{
			if(tmpKmerValueVec_raw[0] == tmpKmerValueVec_raw[1])
			{
				tmpKmerValueVec_filtered.push_back(tmpKmerValueVec_raw[0]);
				tmpKmerValueVec_filtered.push_back(tmpKmerValueVec_raw[1]);
			}
			else
			{
				tmpKmerValueVec_filtered.push_back(0);		
				tmpKmerValueVec_filtered.push_back(0);
			}
		}
		else // tmpQuerySeqKmerNum_1 > 2
		{
			// the 1st
			if(tmpKmerValueVec_raw[0] == tmpKmerValueVec_raw[1])
				tmpKmerValueVec_filtered.push_back(tmpKmerValueVec_raw[0]);
			else
				tmpKmerValueVec_filtered.push_back(0);
			// inter 
			for(int tmp = 1; tmp < tmpQuerySeqKmerNum_1 - 1; tmp++)
			{
				if((tmpKmerValueVec_raw[tmp-1] == tmpKmerValueVec_raw[tmp])
					||(tmpKmerValueVec_raw[tmp] == tmpKmerValueVec_raw[tmp+1]))
					tmpKmerValueVec_filtered.push_back(tmpKmerValueVec_raw[tmp]);
				else
					tmpKmerValueVec_filtered.push_back(0);
			}
			// the last
			if(tmpKmerValueVec_raw[tmpQuerySeqKmerNum_1 - 2] == tmpKmerValueVec_raw[tmpQuerySeqKmerNum_1 - 1])
				tmpKmerValueVec_filtered.push_back(tmpKmerValueVec_raw[tmpQuerySeqKmerNum_1 - 1]);
			else
				tmpKmerValueVec_filtered.push_back(0);			
		}

		// the 2nd set of KmerValueVec
		if(tmpQuerySeqKmerNum_2 == 1)
			tmpKmerValueVec_filtered.push_back(0);
		else if(tmpQuerySeqKmerNum_2 == 2)
		{
			if(tmpKmerValueVec_raw[tmpQuerySeqKmerNum_1] == tmpKmerValueVec_raw[tmpQuerySeqKmerNum_1+1])
			{
				tmpKmerValueVec_filtered.push_back(tmpKmerValueVec_raw[tmpQuerySeqKmerNum_1]);
				tmpKmerValueVec_filtered.push_back(tmpKmerValueVec_raw[tmpQuerySeqKmerNum_1 + 1]);
			}
			else
			{
				tmpKmerValueVec_filtered.push_back(0);		
				tmpKmerValueVec_filtered.push_back(0);
			}
		}
		else // tmpQuerySeqKmerNum_2 > 2
		{
			// the last
			if(tmpKmerValueVec_raw[tmpKmerValueVec_raw_size - 2] == tmpKmerValueVec_raw[tmpKmerValueVec_raw_size - 1])
				tmpKmerValueVec_filtered.push_back(tmpKmerValueVec_raw[tmpKmerValueVec_raw_size - 1]);
			else
				tmpKmerValueVec_filtered.push_back(0);
			// inter 
			for(int tmp = tmpKmerValueVec_raw_size - 2; tmp >= tmpQuerySeqKmerNum_1 + 1; tmp--)
			{
				if((tmpKmerValueVec_raw[tmp - 1] == tmpKmerValueVec_raw[tmp])
					||(tmpKmerValueVec_raw[tmp] == tmpKmerValueVec_raw[tmp + 1]))
					tmpKmerValueVec_filtered.push_back(tmpKmerValueVec_raw[tmp]);
				else
					tmpKmerValueVec_filtered.push_back(0);
			}
			// the 1st
			if(tmpKmerValueVec_raw[tmpQuerySeqKmerNum_1] == tmpKmerValueVec_raw[tmpQuerySeqKmerNum_1 + 1])
				tmpKmerValueVec_filtered.push_back(tmpKmerValueVec_raw[tmpQuerySeqKmerNum_1]);
			else
				tmpKmerValueVec_filtered.push_back(0);			
		}
	}
}

int main(int argc, char** argv)
{
	if((argc != 10)&&(argc != 11))
	{
		cout << "#1 inputLothelloNodeIndex" << endl;
		cout << "#2 outputFolder" << endl; 
		cout << "#3 Kmer_length" << endl;
		cout << "#4 threads_num " << endl;
		cout << "#5 fa_or_fq" << endl; 
		cout << "#6 SE_or_PE" << endl;
		cout << "#7 inputGene2indexFile" << endl;
		cout << "#8 jump_length" << endl;
		cout << "#9 inputReadFile_1" << endl;
		cout << "(#10 inputReadFile_2)" << endl;
		exit(1);
	}
	//int jump_length = 1;
	int normalRecordNum = 1000000;
	string inputLothelloNodeIndexStr = argv[1];
	string outputDirStr = argv[2];
	string kmer_length_str = argv[3];
    string threads_num_str = argv[4];
    string fa_or_fq_str = argv[5];
    string SE_or_PE_str = argv[6];
	string inputGeneIdListFile = argv[7];
	string jump_length_str = argv[8];
	string inputReadFile_1 = argv[9];
	string inputReadFile_2;
	if(argc == 11)
		inputReadFile_2 = argv[10];

	int jump_length = atoi(jump_length_str.c_str());
	outputDirStr += "/";
  	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());
   	string log_file = outputDirStr + "log.txt";
   	ofstream log_ofs(log_file.c_str());

   	string command_line_str = "Command_Line:\n#1 inputLothelloNodeIndex: " + inputLothelloNodeIndexStr
   		+ "\n#2 outputFolder: " + outputDirStr + "\n#3 Kmer_length: " + kmer_length_str + "\n#4 threads_num: "
   		+ threads_num_str + "\n#5 fa_or_fq: " + fa_or_fq_str + "\n#6 SE_or_PE: " + SE_or_PE_str 
   		+ "\n#7 inputSpeciesTaxoIdFile: " + inputGeneIdListFile + " \n#8 jump_length: " + jump_length_str
   		+ "\n#9 inputReadFile_1: " + inputReadFile_1;
   	if(argc == 11)
   		command_line_str = command_line_str + "\n#10 inputReadFile_2: " + inputReadFile_2;
   	log_ofs << command_line_str << endl;
   	cout << command_line_str << endl;

    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "start to load index" << endl;
    int kmer_length = atoi(kmer_length_str.c_str());
    int threads_num = atoi(threads_num_str.c_str());
    int splitbit = 6;
    IOHelper<keyT,valueT> *helper;
    helper = new ConstantLengthKmerHelper<keyT,valueT>(kmer_length,splitbit);
    MulOth<keyT> * moth;
    moth = new MulOth<keyT>(argv[1], helper);  	

    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "end of loading index" << endl;

    #ifdef FUSION_GENE_READ_DETECTION
    FusionCount_Info fusionCountInfo;
    fusionCountInfo.initaite_gene2reissuedFile(inputGeneIdListFile);
    int geneTotalNum = fusionCountInfo.return_geneTotalNum();
    #else
    GeneCount_Info geneCountInfo;
    geneCountInfo.initaite_gene2reissuedFile(inputGeneIdListFile); 
    int geneTotalNum = geneCountInfo.return_geneTotalNum();
    #endif

    //int geneTotalNum = 57363;
    int repetitiveClassId = geneTotalNum + 1;

    bool fa_or_fq_bool, SE_or_PE_bool;
    get_Fa_or_Fq_SE_or_PE_bool(fa_or_fq_bool, SE_or_PE_bool, fa_or_fq_str, SE_or_PE_str, argc);
	bool InputAsFastq = (!fa_or_fq_bool);

    cout << endl << "[" << asctime(local) << "start to do query" << endl;
    string InputReadFile = inputReadFile_1;
	ifstream inputRead_ifs(InputReadFile.c_str());
	string InputReadFile_PE;
	if(!SE_or_PE_bool)
		InputReadFile_PE = inputReadFile_2;
	else
		InputReadFile_PE = "NULL";
	ifstream inputRead_PE_ifs(InputReadFile_PE.c_str());
	
	int readTotalNum = 0;
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

	#ifdef FUSION_GENE_READ_DETECTION
	vector< pair<int,int> > assignment_fusionReissuedIdPairVec(normalRecordNum);
	vector< pair<int,int> > kmerCountPairVec(normalRecordNum);
	vector<int> fusionSiteVec(normalRecordNum);
	vector<int> kmerClusterDistanceVec(normalRecordNum);
	vector< pair<int,int> > kmerCountPairVec_compatible(normalRecordNum);
	vector<int> fusionSiteVec_compatible(normalRecordNum);
	vector<int> kmerClusterDistanceVec_compatible(normalRecordNum);	
	string output_fusion_assignment_file = outputDirStr + "fusion_assignment_all.txt";
	FILE * fusion_assignment_pOut = fopen(output_fusion_assignment_file.c_str(), "w");
	string output_fusion_assignment_file_valid = outputDirStr + "fusion_assignment_valid.txt";
	FILE * fusion_assignment_pOut_valid = fopen(output_fusion_assignment_file_valid.c_str(), "w");
	string output_fusion_assignment_file_invalid = outputDirStr + "fusion_assignment_invalid.txt";
	FILE * fusion_assignment_pOut_invalid = fopen(output_fusion_assignment_file_invalid.c_str(), "w");		
	FusionAssignment_Info fusionAssignmentInfo;
	fusionAssignmentInfo.initiate(1, geneTotalNum, repetitiveClassId);
	#else
	vector<int> assignment_geneReissuedIdVec(normalRecordNum);
	string output_gene_assignment_file = outputDirStr + "gene_assignment_all.txt";
	FILE * gene_assignment_pOut = fopen(output_gene_assignment_file.c_str(), "w");
	string output_gene_assignment_file_valid = outputDirStr + "gene_assignment_all_valid.txt";
	FILE * gene_assignment_pOut_valid = fopen(output_gene_assignment_file_valid.c_str(), "w");	
	string output_gene_assignment_file_invalid = outputDirStr + "gene_assignment_all_invalid.txt";
	FILE * gene_assignment_pOut_invalid = fopen(output_gene_assignment_file_invalid.c_str(), "w");	
	GeneAssignment_Info geneAssignmentInfo;
	geneAssignmentInfo.initiate(1, geneTotalNum, repetitiveClassId);
	#endif
	for(tmpTurn = 0; ; tmpTurn++)
	{
		if(EndOfRecord)
			break;		
		int recordNum = normalRecordNum;
		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
		log_ofs << endl << "[" << asctime(local) << "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
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
		log_ofs << endl << "[" << asctime(local) << "finish reading Fasta/Fastq file, turn: " << tmpTurn + 1 << endl;
		log_ofs << endl << "[" << asctime(local) << "start to fix reads in phase1, turn: " << tmpTurn + 1 << endl;			
		
		omp_set_num_threads(threads_num);
		#pragma omp parallel for schedule(dynamic)
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{
			int threadNO = omp_get_thread_num();
			QuerySeq_Info tmpQuerySeqInfo_1;
			tmpQuerySeqInfo_1.initiate(readName1Vec[tmpOpenMP], readSeq1Vec[tmpOpenMP]);
			
			vector<valueT> tmpKmerValueVec;
			if((readSeq1Vec[tmpOpenMP]).length() >= kmer_length)
				tmpQuerySeqInfo_1.querySeq_returnKmerValueVec(tmpKmerValueVec, moth, kmer_length, jump_length); 
			int tmpQuerySeqKmerNum_1 = tmpKmerValueVec.size();
			if(!SE_or_PE_bool)
			{
				QuerySeq_Info tmpQuerySeqInfo_2;
				tmpQuerySeqInfo_2.initiate(readName2Vec[tmpOpenMP], readSeq2Vec[tmpOpenMP]);
				if((readSeq2Vec[tmpOpenMP]).length() >= kmer_length)
					tmpQuerySeqInfo_2.querySeq_returnKmerValueVec(tmpKmerValueVec, moth, kmer_length, jump_length);
			}

			vector<valueT> tmpKmerValueVec_filtered;
			reorder_removeIsolatedKmerValue(tmpKmerValueVec_filtered, tmpKmerValueVec, tmpQuerySeqKmerNum_1);
			#ifdef PRINT_DEBUG
			cout << "tmpKmerValueVec info: " << endl;
			for(int tmp = 0; tmp < tmpKmerValueVec.size(); tmp ++)
				cout << tmpKmerValueVec[tmp] << ",";
			cout << endl; 
			cout << "tmpKmerValueVec_filtered info: " << endl;
			for(int tmp = 0; tmp < tmpKmerValueVec_filtered.size(); tmp ++)
				cout << tmpKmerValueVec_filtered[tmp] << ",";
			cout << endl;				 
			#endif
			#ifdef FUSION_GENE_READ_DETECTION
			int fusionAssignment_geneReissuedId_1 = -1;
			int fusionAssignment_geneReissuedId_2 = -1;
			int fusionAssignment_KmerCount_1 = 0;
			int fusionAssignment_KmerCount_2 = 0;
			int fusionAssignment_fusionSite = -1;
			int fusionAssignment_kmerClusterDistance = 999;
			int fusionAssignment_KmerCount_1_compatible = 0;
			int fusionAssignment_KmerCount_2_compatible = 0;
			int fusionAssignment_fusionSite_compatible = -1;
			int fusionAssignment_kmerClusterDistance_compatible = 999;

			fusionAssignmentInfo.getFusionAssignment(fusionAssignment_geneReissuedId_1, fusionAssignment_geneReissuedId_2,
				
				fusionAssignment_KmerCount_1, fusionAssignment_KmerCount_2, 
				fusionAssignment_fusionSite, fusionAssignment_kmerClusterDistance,
				
				fusionAssignment_KmerCount_1_compatible, fusionAssignment_KmerCount_2_compatible, 
				fusionAssignment_fusionSite_compatible, fusionAssignment_kmerClusterDistance_compatible,
				
				tmpKmerValueVec_filtered, tmpQuerySeqKmerNum_1, kmer_length);

			assignment_fusionReissuedIdPairVec[tmpOpenMP].first = fusionAssignment_geneReissuedId_1;
			assignment_fusionReissuedIdPairVec[tmpOpenMP].second = fusionAssignment_geneReissuedId_2;
			
			kmerCountPairVec[tmpOpenMP].first = fusionAssignment_KmerCount_1;
			kmerCountPairVec[tmpOpenMP].second = fusionAssignment_KmerCount_2;
			fusionSiteVec[tmpOpenMP] = fusionAssignment_fusionSite;
			kmerClusterDistanceVec[tmpOpenMP] = fusionAssignment_kmerClusterDistance;

			kmerCountPairVec_compatible[tmpOpenMP].first = fusionAssignment_KmerCount_1_compatible;
			kmerCountPairVec_compatible[tmpOpenMP].second = fusionAssignment_KmerCount_2_compatible;
			fusionSiteVec_compatible[tmpOpenMP] = fusionAssignment_fusionSite_compatible;
			kmerClusterDistanceVec_compatible[tmpOpenMP] = fusionAssignment_kmerClusterDistance_compatible;			
			#ifdef PRINT_DEBUG
			cout << "fusionAssignment_geneReissuedId_1: " << fusionAssignment_geneReissuedId_1 << endl;
			cout << "fusionAssignment_geneReissuedId_2: " << fusionAssignment_geneReissuedId_2 << endl;
			cout << "fusionAssignment_KmerCount_1: " << fusionAssignment_KmerCount_1 << endl;
			cout << "fusionAssignment_KmerCount_2: " << fusionAssignment_KmerCount_2 << endl;
			cout << "fusionAssignment_fusionSite: " << fusionAssignment_fusionSite << endl;
			cout << "fusionAssignment_KmerCount_1_compatible: " << fusionAssignment_KmerCount_1_compatible << endl;
			cout << "fusionAssignment_KmerCount_2_compatible: " << fusionAssignment_KmerCount_2_compatible << endl;
			cout << "fusionAssignment_fusionSite_compatible: " << fusionAssignment_fusionSite_compatible << endl;
			cout << "fusionAssignment_kmerClusterDistance_compatible: " << fusionAssignment_kmerClusterDistance_compatible << endl;
			#endif

			#else
			int tmpGeneAssignementId = geneAssignmentInfo.getGeneAssignment(tmpKmerValueVec_filtered);
			assignment_geneReissuedIdVec[tmpOpenMP] = tmpGeneAssignementId;
			#endif
		}
		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		cout << endl << "[" << asctime(local) << "start to output ... turn: " << tmpTurn+1 << endl;
		log_ofs << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		log_ofs << endl << "[" << asctime(local) << "start to output ... turn: " << tmpTurn+1 << endl;
		//cout << "realRecordNum: " << realRecordNum << endl;
		for(int tmp = 0; tmp < realRecordNum; tmp++)
		{
			#ifdef FUSION_GENE_READ_DETECTION
			int tmpFusionGenePairAssignmentId_1 = assignment_fusionReissuedIdPairVec[tmp].first;
			int tmpFusionGenePairAssignmentId_2 = assignment_fusionReissuedIdPairVec[tmp].second;
			
			int tmpKmerCount_1 = kmerCountPairVec[tmp].first;
			int tmpKmerCount_2 = kmerCountPairVec[tmp].second;
			int tmpFusionSite = fusionSiteVec[tmp];
			int tmpKmerClusterDistance = kmerClusterDistanceVec[tmp];

			int tmpKmerCount_1_compatible = kmerCountPairVec_compatible[tmp].first;
			int tmpKmerCount_2_compatible = kmerCountPairVec_compatible[tmp].second;
			int tmpFusionSite_compatible = fusionSiteVec_compatible[tmp];
			int tmpKmerClusterDistance_compatible = kmerClusterDistanceVec_compatible[tmp];

			fprintf(fusion_assignment_pOut, "%s\t%d\t%d\t%s\t%s\n", readName1Vec[tmp].c_str(), 
				tmpFusionGenePairAssignmentId_1, tmpFusionGenePairAssignmentId_2, 
				readSeq1Vec[tmp].c_str(), readSeq2Vec[tmp].c_str());
			if(((tmpFusionGenePairAssignmentId_1 < 1)||(tmpFusionGenePairAssignmentId_1 > geneTotalNum))
				&&((tmpFusionGenePairAssignmentId_2 < 1)||(tmpFusionGenePairAssignmentId_2 > geneTotalNum)))
				fprintf(fusion_assignment_pOut_invalid, "%s\t%s\t%s\n", readName1Vec[tmp].c_str(), 
					readSeq1Vec[tmp].c_str(), readSeq2Vec[tmp].c_str());
			else
			{
				fusionCountInfo.addFusionCount(tmpFusionGenePairAssignmentId_1, tmpFusionGenePairAssignmentId_2);
				fprintf(fusion_assignment_pOut_valid, "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n", readName1Vec[tmp].c_str(), 
					(fusionCountInfo.return_geneIdStr(tmpFusionGenePairAssignmentId_1)).c_str(), 
					(fusionCountInfo.return_geneIdStr(tmpFusionGenePairAssignmentId_2)).c_str(),
					tmpKmerCount_1, tmpKmerCount_2, tmpFusionSite, tmpKmerClusterDistance,
					tmpKmerCount_1_compatible, tmpKmerCount_2_compatible, tmpFusionSite_compatible, tmpKmerClusterDistance_compatible,
					readSeq1Vec[tmp].c_str(), readSeq2Vec[tmp].c_str());
			}
			#else
			int tmpGeneAssignementId = assignment_geneReissuedIdVec[tmp];
			fprintf(gene_assignment_pOut, "%s\t%d\t%s\t%s\n", readName1Vec[tmp].c_str(), tmpGeneAssignementId,
				readSeq1Vec[tmp].c_str(), readSeq2Vec[tmp].c_str());
			if((tmpGeneAssignementId < 1)||(tmpGeneAssignementId > geneTotalNum))
			{
				geneCountInfo.addUnassignedReadCount();
				fprintf(gene_assignment_pOut_invalid, "%s\t%s\t%s\n", readName1Vec[tmp].c_str(),
					readSeq1Vec[tmp].c_str(), readSeq2Vec[tmp].c_str());
			}
			else
			{
				geneCountInfo.addGeneCount(tmpGeneAssignementId);
				fprintf(gene_assignment_pOut_valid, "%s\t%s\t%s\t%s\n", readName1Vec[tmp].c_str(), 
					(geneCountInfo.return_geneIdStr(tmpGeneAssignementId)).c_str(),
					readSeq1Vec[tmp].c_str(), readSeq2Vec[tmp].c_str());	
			}
			#endif
		}
	}
	#ifdef FUSION_GENE_READ_DETECTION
	string fusionCount_file = outputDirStr + "fusionCount.txt";
	fusionCountInfo.printFusionCountInfo(fusionCount_file);
	fclose(fusion_assignment_pOut);
	fclose(fusion_assignment_pOut_valid);
	fclose(fusion_assignment_pOut_invalid);
	#else
	string geneCount_file = outputDirStr + "geneCount.txt";
	string assignStats_file = outputDirStr + "assignStats.txt";
	geneCountInfo.printGeneCountInfo(geneCount_file);
	geneCountInfo.printAssignStatsInfo(assignStats_file);
	fclose(gene_assignment_pOut);
	fclose(gene_assignment_pOut_valid);
	fclose(gene_assignment_pOut_invalid);
	#endif
	inputRead_ifs.close();
	inputRead_PE_ifs.close();
    nowtime = time(NULL);
    local = localtime(&nowtime);    
    cout << endl << "[" << asctime(local) << "end of doing query" << endl;
    log_ofs << endl << "[" << asctime(local) << "end of doing query" << endl;
    log_ofs.close();
    delete helper;
    delete moth;
    return 0;
}