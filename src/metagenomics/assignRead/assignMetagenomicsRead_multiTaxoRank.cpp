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
//#include "../../genomeRegionAssignment/general/regionAssignment_info.h"
//#include "../general/reissuedGenomeID2TaxoID_info.h"
#include "../general/NCBIfullTaxoID2Name_info.h"
#include "../general/bacterialTaxo_info.h"
#include "../general/reissuedGenomeID2TaxoID_info.h"
#include "../general/MSKW_info.h"
#include "../general/taxoClassAssignment_info.h"

using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if((argc != 13)&&(argc != 14))
	{
		cout << "#1 inputLothelloNodeIndex" << endl;
		cout << "#2 outputFolder" << endl; 
		cout << "#3 Kmer_length" << endl;
		cout << "#4 threads_num " << endl;
		cout << "#5 fa_or_fq" << endl; 
		cout << "#6 SE_or_PE" << endl;
		cout << "#7 inputSpeciesTaxoIdFile" << endl;
		cout << "#8 NCBIfullTaxoId2NameFile" << endl;
		cout << "#9 Rank_Level(Phylum/Genus/Species)" << endl;
		cout << "#10 assignment_KmerCount_min" << endl;
		cout << "#11 assignment_KmerCountSquareSum_min" << endl;
		cout << "#12 inputReadFile_1" << endl;
		cout << "(#13 inputReadFile_2)" << endl;
		exit(1);
	}
	int jump_length = 1;
	int normalRecordNum = 500000;//0;
	int taxo_rank;

	string outputDirStr = argv[2];
	string kmer_length_str = argv[3];
    string threads_num_str = argv[4];
    string fa_or_fq_str = argv[5];
    string SE_or_PE_str = argv[6];
	string inputSpeciesTaxoIdFile = argv[7];
	string NCBIfullTaxoId2NameFile = argv[8];
	string taxo_level_name = argv[9];
	string assignment_KmerCount_min_str = argv[10];
	string assignment_KmerCountSquareSum_min_str = argv[11];

	double assignment_KmerCount_min = atof(assignment_KmerCount_min_str.c_str()); //0.05;
	double assignment_KmerCountSquareSum_min = atof(assignment_KmerCountSquareSum_min_str.c_str()); //0.5;
	//double assignment_score_min = 0.6;

	if(taxo_level_name == "Phylum")
		taxo_rank = 3;
	else if(taxo_level_name == "Class")
		taxo_rank = 4;
	else if(taxo_level_name == "Order")
		taxo_rank = 5;
	else if(taxo_level_name == "Family")
		taxo_rank = 6;
	else if(taxo_level_name == "Genus")
		taxo_rank = 7;
	else if(taxo_level_name == "Species")
		taxo_rank = 8;
	else
	{
		cout << "invalid taxo_level_name: " << taxo_level_name << endl;
		exit(1);
	}

	outputDirStr += "/";
  	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());
   	string log_file = outputDirStr + "log.txt";
   	ofstream log_ofs(log_file.c_str());
   	string stats_file = outputDirStr + "stats.txt";
   	ofstream stats_ofs(stats_file.c_str());   	

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

	BacterialTaxo_Info bacterialTaxoInfo;
	bacterialTaxoInfo.initiate_bacterialTaxoFile_NCBIfullTaxoId2NameFile(
		inputSpeciesTaxoIdFile, NCBIfullTaxoId2NameFile);
	bacterialTaxoInfo.reissueTaxoIdName_all();
	bacterialTaxoInfo.initiate_lowerRankTaxoReissuedId_to_higherRankTaxoReissuedId_array();
	string outputDir_taxoInfo = outputDirStr + "taxoInfo/";
	bacterialTaxoInfo.print(outputDir_taxoInfo);

	int taxo_num_specific = bacterialTaxoInfo.return_taxo_num(taxo_rank);
    int taxo_num_total = bacterialTaxoInfo.return_taxo_num_total();
    int repetitive_class_Id = taxo_num_total + 1;
    int alien_class_Id = 0;

    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "end of initiating bacterialTaxoInfo" << endl;

	vector< vector<unsigned long long> > taxoCountVecVec_thread;
	vector<unsigned long long> repetitiveCountVec_thread;
	for(int tmp = 0; tmp < threads_num; tmp++)
	{
		vector<unsigned long long> tmp_taxoCountVec_thread;
		for(int tmp2 = 0; tmp2 < taxo_num_specific; tmp2++)
			tmp_taxoCountVec_thread.push_back(0);
		taxoCountVecVec_thread.push_back(tmp_taxoCountVec_thread);
		repetitiveCountVec_thread.push_back(0);
	}

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

    if((SE_or_PE_bool && (argc != 13)) || ((!SE_or_PE_bool) && (argc != 14)))
    {
    	cout << "inconsistent SE_or_PE_bool and argc" << endl;
    	cout << "SE_or_PE_bool: " << SE_or_PE_bool << endl;
    	cout << "argc: " << argc << endl;
    	exit(1);
    }
    
    cout << endl << "[" << asctime(local) << "start to do query" << endl;
    string InputReadFile = argv[12];
	ifstream inputRead_ifs(InputReadFile.c_str());
	string InputReadFile_PE;
	if(!SE_or_PE_bool)
		InputReadFile_PE = argv[13];
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

	//vector<int> readAssignedTaxoClassVec(normalRecordNum);
	//vector<string> readAssignedTaxoClassNameVec(normalRecordNum);
	vector<int> readAssignedTaxoClassIdVec(normalRecordNum);
	vector<int> readAssignedTaxoClassIdVec_squareSum(normalRecordNum);
	vector<int> readAssignedTaxoClassIdVec_KmerCount(normalRecordNum);
	vector<int> readAssignedTaxoClassIdVec_KmerCountSquareSum(normalRecordNum);

	#ifdef ASSIGN_INFO
	if(threads_num != 1)
	{
		cout << "if ASSIGN_INFO defined, threads_num should be 1" << endl;
		cout << "threads_num: " << threads_num << endl;
		exit(1);
	}
	string output_taxo_assignment_file_assignInfo = outputDirStr + "taxo_assignment.txt.assignInfo";
	ofstream taxo_assignment_ofs_assignInfo(output_taxo_assignment_file_assignInfo.c_str());
	string output_taxo_assignment_file_assignInfo_detail = outputDirStr + "taxo_assignment.txt.assignInfo.detail";
	ofstream taxo_assignment_ofs_assignInfo_detail(output_taxo_assignment_file_assignInfo_detail.c_str());
	#endif
	string output_taxo_assignment_file = outputDirStr + "taxo_assignment.txt";
	ofstream taxo_assignment_ofs(output_taxo_assignment_file.c_str());

	TaxoClassAssignment_Info taxoClassAssignmentInfo;
	taxoClassAssignmentInfo.initiate_multiTaxoRank(bacterialTaxoInfo);
	taxoClassAssignmentInfo.initiate_groupIdPair_inLineOrNot_array(bacterialTaxoInfo);
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

			if(!SE_or_PE_bool)
			{
				QuerySeq_Info tmpQuerySeqInfo_2;
				tmpQuerySeqInfo_2.initiate(readName2Vec[tmpOpenMP], readSeq2Vec[tmpOpenMP]);
				if((readSeq2Vec[tmpOpenMP]).length() >= kmer_length)
					tmpQuerySeqInfo_2.querySeq_returnKmerValueVec(tmpKmerValueVec, moth, kmer_length, jump_length);
			}

			int tmpInvalidKmerCount = 0;
			int tmpRepetitiveKmerCount = 0;
			int tmpDiscriminativeKmerCount = 0;
			MaximalSpecificKmerWindow_Info tmpMSKWinfo;
			taxoClassAssignmentInfo.get_taxoId2windowVecMap_specific(tmpMSKWinfo, tmpKmerValueVec, kmer_length, 
				taxo_rank, bacterialTaxoInfo, tmpInvalidKmerCount, tmpRepetitiveKmerCount, tmpDiscriminativeKmerCount);

			int tmp_taxo_reissuedId_best, tmpMaximalWindowKmerCount_best, 
				tmp_taxo_reissuedId_secondBest, tmpMaximalWindowKmerCount_secondBest,
				tmp_taxo_reissuedId_best_squareSum, tmpMaximalWindowKmerCount_best_squareSum,
				tmp_taxo_reissuedId_secondBest_squareSum, tmpMaximalWindowKmerCount_secondBest_squareSum;
			tmpMSKWinfo.get_taxoReissuedId_maximalWindow_best_secondBest(
				tmp_taxo_reissuedId_best, tmpMaximalWindowKmerCount_best, 
				tmp_taxo_reissuedId_secondBest, tmpMaximalWindowKmerCount_secondBest, 
				assignment_KmerCount_min);
			tmpMSKWinfo.get_taxoReissuedId_maximalWindowSquareSum_best_secondBest(
				tmp_taxo_reissuedId_best_squareSum, tmpMaximalWindowKmerCount_best_squareSum,
				tmp_taxo_reissuedId_secondBest_squareSum, tmpMaximalWindowKmerCount_secondBest_squareSum,
				assignment_KmerCountSquareSum_min);
			//cout << "tmp_taxo_reissuedId_best: " << tmp_taxo_reissuedId_best << endl;
			//cout << "tmp_taxo_reissuedId_secondBest: " << tmp_taxo_reissuedId_secondBest << endl;
			if((tmpMaximalWindowKmerCount_best == tmpMaximalWindowKmerCount_secondBest)
				&&(tmp_taxo_reissuedId_best != tmp_taxo_reissuedId_secondBest))
				tmp_taxo_reissuedId_best = -2;
			if(tmpMaximalWindowKmerCount_best_squareSum == tmpMaximalWindowKmerCount_secondBest_squareSum)
				tmp_taxo_reissuedId_best_squareSum = -2;
			//cout << "tmp_taxo_reissuedId_best: " << tmp_taxo_reissuedId_best << endl;
			//cout << "tmp_taxo_reissuedId_secondBest: " << tmp_taxo_reissuedId_secondBest << endl;
			int tmpTaxoClassId = bacterialTaxoInfo.return_taxo_Id(taxo_rank, tmp_taxo_reissuedId_best);
			int tmpTaxoClassId_squareSum = bacterialTaxoInfo.return_taxo_Id(taxo_rank, tmp_taxo_reissuedId_best_squareSum);

			readAssignedTaxoClassIdVec[tmpOpenMP] = tmpTaxoClassId;
			readAssignedTaxoClassIdVec_squareSum[tmpOpenMP] = tmpTaxoClassId_squareSum;

			readAssignedTaxoClassIdVec_KmerCount[tmpOpenMP] = tmpMaximalWindowKmerCount_best;
			readAssignedTaxoClassIdVec_KmerCountSquareSum[tmpOpenMP] = tmpMaximalWindowKmerCount_best_squareSum;			
		
			#ifdef ASSIGN_INFO
			taxo_assignment_ofs_assignInfo_detail << endl << readName1Vec[tmpOpenMP] << endl << readSeq1Vec[tmpOpenMP] << endl << readQualSeq1Vec[tmpOpenMP] << endl;
			taxo_assignment_ofs_assignInfo_detail << "GroupId info:" << endl;
			if(!SE_or_PE_bool)
				taxo_assignment_ofs_assignInfo_detail << readName2Vec[tmpOpenMP] << endl << readSeq2Vec[tmpOpenMP] << endl << readQualSeq2Vec[tmpOpenMP] << endl;
			for(int tmp = 0; tmp < tmpKmerValueVec.size(); tmp++)
				taxo_assignment_ofs_assignInfo_detail << tmpKmerValueVec[tmp] << ",";
			taxo_assignment_ofs_assignInfo_detail << endl;
			vector<int> tmpSpecificRankTaxoReissuedIdVec;
			taxoClassAssignmentInfo.convertGroupIdVec2reissuedIdVec(tmpKmerValueVec, tmpSpecificRankTaxoReissuedIdVec, taxo_rank, bacterialTaxoInfo);
			taxo_assignment_ofs_assignInfo_detail << endl << "Taxo Reissued id: " << endl;
			for(int tmp = 0; tmp < tmpSpecificRankTaxoReissuedIdVec.size(); tmp++)
				taxo_assignment_ofs_assignInfo_detail << tmpSpecificRankTaxoReissuedIdVec[tmp] << ",";
			taxo_assignment_ofs_assignInfo_detail << endl;
			vector<int> tmpRawTaxoId;
			bacterialTaxoInfo.get_taxoIdVec_from_taxoReissuedIdVec(tmpRawTaxoId, tmpSpecificRankTaxoReissuedIdVec, taxo_rank);
			taxo_assignment_ofs_assignInfo_detail << "Taxo raw id: " << endl;
			for(int tmp = 0; tmp < tmpRawTaxoId.size(); tmp++)
				taxo_assignment_ofs_assignInfo_detail << tmpRawTaxoId[tmp] << ",";
			taxo_assignment_ofs_assignInfo_detail << endl;
			string tmpStr_taxoId2windowMap = tmpMSKWinfo.print_taxoId2windowMap_2_string();
			taxo_assignment_ofs_assignInfo_detail << "TaxoId2windowMap: " << endl << tmpStr_taxoId2windowMap << endl;
			taxo_assignment_ofs_assignInfo_detail << "Specific:\t" << tmpDiscriminativeKmerCount 
				<< "\tRepetitive:\t" << tmpRepetitiveKmerCount << "\tInvalid:\t" << tmpInvalidKmerCount << endl;;
			taxo_assignment_ofs_assignInfo_detail << "Best_maxKmerCount        : " << tmpTaxoClassId << "-" << tmpMaximalWindowKmerCount_best << endl;
			taxo_assignment_ofs_assignInfo_detail << "Best_maxKmerCountSqureSum: " << tmpTaxoClassId_squareSum << "-" << tmpMaximalWindowKmerCount_best_squareSum << endl;
			#endif
		}
		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		cout << endl << "[" << asctime(local) << "start to output ... turn: " << tmpTurn+1 << endl;
		log_ofs << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		log_ofs << endl << "[" << asctime(local) << "start to output ... turn: " << tmpTurn+1 << endl;			
		for(int tmp = 0; tmp < realRecordNum; tmp++)
			taxo_assignment_ofs << readName1Vec[tmp].c_str() << "\t" << readAssignedTaxoClassIdVec[tmp] 
				<< "\t" << readAssignedTaxoClassIdVec_squareSum[tmp] << endl;
		
		#ifdef ASSIGN_INFO
		for(int tmp = 0; tmp < realRecordNum; tmp++)
			taxo_assignment_ofs_assignInfo << readName1Vec[tmp].c_str() << "\t" 
				<< readAssignedTaxoClassIdVec[tmp] << "\t" << readAssignedTaxoClassIdVec_KmerCount[tmp] << "\t" 
				<< readAssignedTaxoClassIdVec_squareSum[tmp] << "\t" << readAssignedTaxoClassIdVec_KmerCountSquareSum[tmp] << endl;
		#endif
	}

	#ifdef ASSIGN_INFO
	taxo_assignment_ofs_assignInfo_detail.close();
	taxo_assignment_ofs_assignInfo.close();
	#endif
	taxo_assignment_ofs.close();

	// vector<unsigned long long> taxoCountVec;
	// for(int tmp = 0; tmp < taxo_num_specific; tmp++)
	// 	taxoCountVec.push_back(0);
	// for(int tmp = 0; tmp < threads_num; tmp++)
	// {
	// 	for(int tmp2 = 0; tmp2 < taxo_num_specific; tmp2++)
	// 		taxoCountVec[tmp2] += (taxoCountVecVec_thread[tmp])[tmp2];
	// }
	// for(int tmp = 0; tmp < taxo_num_specific; tmp++)
	// 	stats_ofs << tmp << "\t" << bacterialTaxoInfo.return_taxo_name(taxo_rank, tmp) << "\t" << taxoCountVec[tmp] << endl;

	// unsigned long long repetitiveCount = 0;
	// for(int tmp = 0; tmp < threads_num; tmp++)
	// 	repetitiveCount += repetitiveCountVec_thread[tmp];
	// stats_ofs << endl << "Repetitive_Count:\t" << repetitiveCount << endl;

	inputRead_ifs.close();
	inputRead_PE_ifs.close();
    nowtime = time(NULL);
    local = localtime(&nowtime);    
    cout << endl << "[" << asctime(local) << "end of doing query" << endl;
    log_ofs << endl << "[" << asctime(local) << "end of doing query" << endl;
    log_ofs.close();
    stats_ofs.close();
    delete helper;
    delete moth;
    return 0;
}