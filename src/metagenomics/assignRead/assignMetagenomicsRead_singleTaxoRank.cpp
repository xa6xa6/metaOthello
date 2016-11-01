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
#include "../general/taxoClassAssignment_info.h"

using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if((argc != 11)&&(argc != 12))
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
		cout << "#10 inputReadFile_1" << endl;
		cout << "(#11 inputReadFile_2)" << endl;
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
	string outputDir_taxoInfo = outputDirStr + "taxoInfo/";
	bacterialTaxoInfo.print(outputDir_taxoInfo);

    int total_taxo_num = bacterialTaxoInfo.return_taxo_num(taxo_rank);
    int taxo_num = total_taxo_num;
    int repetitiveClassId = taxo_num + 1;
    int repetitive_class_id = total_taxo_num + 1;
    int alienClassId = 0;

    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "end of initiating bacterialTaxoInfo" << endl;

	vector< vector<unsigned long long> > taxoCountVecVec_thread;
	vector<unsigned long long> repetitiveCountVec_thread;
	for(int tmp = 0; tmp < threads_num; tmp++)
	{
		vector<unsigned long long> tmp_taxoCountVec_thread;
		for(int tmp2 = 0; tmp2 < total_taxo_num; tmp2++)
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

    if((SE_or_PE_bool && (argc != 11)) || ((!SE_or_PE_bool) && (argc != 12)))
    {
    	cout << "inconsistent SE_or_PE_bool and argc" << endl;
    	cout << "SE_or_PE_bool: " << SE_or_PE_bool << endl;
    	cout << "argc: " << argc << endl;
    	exit(1);
    }
    
    cout << endl << "[" << asctime(local) << "start to do query" << endl;
    string InputReadFile = argv[10];
	ifstream inputRead_ifs(InputReadFile.c_str());
	string InputReadFile_PE;
	if(!SE_or_PE_bool)
		InputReadFile_PE = argv[11];
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
	#ifdef ASSIGN_INFO
	if(threads_num != 1)
	{
		cout << "if ASSIGN_INFO defined, threads_num should be 1" << endl;
		cout << "threads_num: " << threads_num << endl;
		exit(1);
	}
	vector<int> invalidKmerCountVec(normalRecordNum);
	vector<int> repetitiveKmerCountVec(normalRecordNum);
	vector<int> discriminativeKmerCountVec(normalRecordNum);
	vector<int> KmerCountVec_best(normalRecordNum);
	vector<int> KmerCountVec_secondBest(normalRecordNum);
	vector<int> KmerIndexVec_best(normalRecordNum);
	vector<int> KmerIndexVec_secondBest(normalRecordNum);
	string output_taxo_assignment_file_assignInfo = outputDirStr + "taxo_assignment.txt.assignInfo";
	ofstream taxo_assignment_ofs_assignInfo(output_taxo_assignment_file_assignInfo.c_str());
	string output_taxo_assignment_file_assignInfo_detail = outputDirStr + "taxo_assignment.txt.assignInfo.detail";
	ofstream taxo_assignment_ofs_assignInfo_detail(output_taxo_assignment_file_assignInfo_detail.c_str());
	#endif
	string output_taxo_assignment_file = outputDirStr + "taxo_assignment.txt";
	ofstream taxo_assignment_ofs(output_taxo_assignment_file.c_str());

	TaxoClassAssignment_Info taxoClassAssignmentInfo;
	taxoClassAssignmentInfo.initiate_singleTaxoRank(bacterialTaxoInfo, taxo_rank);
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
			vector<valueT> tmpKmerValueVec_1;
			tmpQuerySeqInfo_1.querySeq_returnKmerValueVec(
				tmpKmerValueVec_1, moth, kmer_length, jump_length); 

			int tmpInvalidKmerCount = 0; 
			int tmpRepetitiveKmerCount = 0; 
			int tmpDiscriminativeKmerCount = 0;
			int tmpKmerCount_best = 0;
			int tmpKmerCount_secondBest = 0;
			int tmpKmerIndex_best = -1;
			int tmpKmerIndex_secondBest = -1;

			vector< pair<int,int> > tmpTaxoIdCountPairVec;

			if(SE_or_PE_bool)
			{	
				taxoClassAssignmentInfo.get_rawTaxoIdCountPairVec_from_rawKmerValueVec_SE(tmpKmerValueVec_1, 
					tmpInvalidKmerCount, tmpRepetitiveKmerCount, tmpDiscriminativeKmerCount, tmpTaxoIdCountPairVec);
			}
			else
			{
				QuerySeq_Info tmpQuerySeqInfo_2;
				tmpQuerySeqInfo_2.initiate(readName2Vec[tmpOpenMP], readSeq2Vec[tmpOpenMP]);
				vector<valueT> tmpKmerValueVec_2;
				tmpQuerySeqInfo_2.querySeq_returnKmerValueVec(tmpKmerValueVec_2, moth, kmer_length, jump_length);
				taxoClassAssignmentInfo.get_rawTaxoIdCountPairVec_from_rawKmerValueVec_PE(tmpKmerValueVec_1, tmpKmerValueVec_2, 
					tmpInvalidKmerCount, tmpRepetitiveKmerCount, tmpDiscriminativeKmerCount, tmpTaxoIdCountPairVec);

				#ifdef ASSIGN_INFO
				taxo_assignment_ofs_assignInfo_detail << endl << readName1Vec[tmpOpenMP] << endl << readSeq1Vec[tmpOpenMP] 
					<< endl << readQualSeq1Vec[tmpOpenMP] << endl << readName2Vec[tmpOpenMP] << endl << readSeq2Vec[tmpOpenMP] 
					<< endl << readQualSeq2Vec[tmpOpenMP] << endl;
				for(int tmp = 0; tmp < tmpKmerValueVec_1.size(); tmp++)
					taxo_assignment_ofs_assignInfo_detail << tmpKmerValueVec_1[tmp] << ",";
				for(int tmp = 0; tmp < tmpKmerValueVec_2.size(); tmp++)
					taxo_assignment_ofs_assignInfo_detail << tmpKmerValueVec_2[tmp] << ",";
				taxo_assignment_ofs_assignInfo_detail << endl;
				#endif
			}
			taxoClassAssignmentInfo.get_taxoIdCount_best_secondBest(tmpTaxoIdCountPairVec,
				tmpInvalidKmerCount, tmpRepetitiveKmerCount, tmpDiscriminativeKmerCount, 
				tmpKmerIndex_best, tmpKmerCount_best, tmpKmerIndex_secondBest, tmpKmerCount_secondBest);
			int tmp_taxo_reissuedId_best = tmpKmerIndex_best - 1;
			if(tmp_taxo_reissuedId_best < 0)
				repetitiveCountVec_thread[threadNO] ++;
			else
				(taxoCountVecVec_thread[threadNO])[tmp_taxo_reissuedId_best] ++;
			#ifdef ASSIGN_INFO
			invalidKmerCountVec[tmpOpenMP] = tmpInvalidKmerCount;
			repetitiveKmerCountVec[tmpOpenMP] = tmpRepetitiveKmerCount;
			discriminativeKmerCountVec[tmpOpenMP] = tmpDiscriminativeKmerCount;
			KmerCountVec_best[tmpOpenMP] = tmpKmerCount_best;
			KmerCountVec_secondBest[tmpOpenMP] = tmpKmerCount_secondBest;
			KmerIndexVec_best[tmpOpenMP] = tmpKmerIndex_best;
			KmerIndexVec_secondBest[tmpOpenMP] = tmpKmerIndex_secondBest;

			taxo_assignment_ofs_assignInfo_detail << endl;
			taxo_assignment_ofs_assignInfo_detail << KmerIndexVec_best[tmpOpenMP] << "-" << KmerCountVec_best[tmpOpenMP] 
				<< "\t" << KmerIndexVec_secondBest[tmpOpenMP] << "-" << KmerCountVec_secondBest[tmpOpenMP] << "\t" 
				<< discriminativeKmerCountVec[tmpOpenMP] << "," << repetitiveKmerCountVec[tmpOpenMP] << "," << invalidKmerCountVec[tmpOpenMP] << endl;
			#endif
			//string tmpTaxoClassName = bacterialTaxoInfo.return_taxo_name(taxo_rank, tmp_taxo_reissuedId_best);
			//readAssignedTaxoClassNameVec[tmpOpenMP] = tmpTaxoClassName;
			int tmpTaxoClassId = bacterialTaxoInfo.return_taxo_Id(taxo_rank, tmp_taxo_reissuedId_best);
			readAssignedTaxoClassIdVec[tmpOpenMP] = tmpTaxoClassId;			
		}
		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		cout << endl << "[" << asctime(local) << "start to output ... turn: " << tmpTurn+1 << endl;
		log_ofs << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		log_ofs << endl << "[" << asctime(local) << "start to output ... turn: " << tmpTurn+1 << endl;			
		for(int tmp = 0; tmp < realRecordNum; tmp++)
			taxo_assignment_ofs << readName1Vec[tmp].c_str() << "\t" << readAssignedTaxoClassIdVec[tmp] << endl;
		
		#ifdef ASSIGN_INFO
		for(int tmp = 0; tmp < realRecordNum; tmp++)
			taxo_assignment_ofs_assignInfo << readName1Vec[tmp].c_str() << "\t" << readAssignedTaxoClassIdVec[tmp]
				<< "\t" << KmerIndexVec_best[tmp] << "-" << KmerCountVec_best[tmp] << "\t"
				<< KmerIndexVec_secondBest[tmp] << "-" << KmerCountVec_secondBest[tmp] << "\t"
				<< discriminativeKmerCountVec[tmp] << "," << repetitiveKmerCountVec[tmp] << "," << invalidKmerCountVec[tmp] << endl;		
		#endif
	}

	#ifdef ASSIGN_INFO
	taxo_assignment_ofs_assignInfo_detail.close();
	taxo_assignment_ofs_assignInfo.close();
	#endif
	taxo_assignment_ofs.close();

	vector<unsigned long long> taxoCountVec;
	for(int tmp = 0; tmp < taxo_num; tmp++)
		taxoCountVec.push_back(0);
	for(int tmp = 0; tmp < threads_num; tmp++)
	{
		for(int tmp2 = 0; tmp2 < taxo_num; tmp2++)
			taxoCountVec[tmp2] += (taxoCountVecVec_thread[tmp])[tmp2];
	}
	for(int tmp = 0; tmp < taxo_num; tmp++)
		stats_ofs << tmp << "\t" << bacterialTaxoInfo.return_taxo_name(taxo_rank, tmp) << "\t" << taxoCountVec[tmp] << endl;

	unsigned long long repetitiveCount = 0;
	for(int tmp = 0; tmp < threads_num; tmp++)
		repetitiveCount += repetitiveCountVec_thread[tmp];
	stats_ofs << endl << "Repetitive_Count:\t" << repetitiveCount << endl;

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