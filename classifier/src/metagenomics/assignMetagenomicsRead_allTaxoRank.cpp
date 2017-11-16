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
#include "../othello/othello.h"
#include "../othello/muloth.h"
#include <chrono>
#include <inttypes.h>
#include "../othello/io_helper.h"
#include "../mps3Lib/read_block_test.h"
#include "../mps3Lib/otherFunc.h"
#include "../mps3Lib/index_info.h"
#include "../query/general/queryConstantDef.h"
#include "../query/general/querySeq_info.h"
#include "general/NCBIfullTaxoID2Name_info.h"
#include "general/bacterialTaxo_info.h"
#include "general/reissuedGenomeID2TaxoID_info.h"
#include "general/MSKW_info.h"
#include "general/taxoClassAssignment_info.h"

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

void print_groupId_rawTaxoIdVecVec(ofstream& taxo_assignment_ofs_assignInfo_detail, 
	BacterialTaxo_Info& bacterialTaxoInfo,
	TaxoClassAssignment_Info& taxoClassAssignmentInfo, bool SE_or_PE_bool,
	string& tmp_readName_1, string& tmp_readSeq_1, string& tmp_readQual_1, 
	string& tmp_readName_2, string& tmp_readSeq_2, string& tmp_readQual_2, 
	vector<valueT>& tmpKmerValueVec, vector<valueT>& tmpKmerValueVec_filtered, 
	vector<int>& tmpDiscriminativeKmerCountVec, vector<int>& tmpRepetitiveKmerCountVec, vector<int>& tmpInvalidKmerCountVec)
{
	taxo_assignment_ofs_assignInfo_detail << endl << tmp_readName_1 << endl << tmp_readSeq_1 << endl << tmp_readQual_1 << endl;
	if(!SE_or_PE_bool)
		taxo_assignment_ofs_assignInfo_detail << tmp_readName_2 << endl << tmp_readSeq_2 << endl << tmp_readQual_2 << endl;
	taxo_assignment_ofs_assignInfo_detail << "GroupId info:" << endl;
	for(int tmp = 0; tmp < tmpKmerValueVec.size(); tmp++)
		taxo_assignment_ofs_assignInfo_detail << tmpKmerValueVec[tmp] << ",";
	taxo_assignment_ofs_assignInfo_detail << endl;

	taxo_assignment_ofs_assignInfo_detail << "Filtred GroupId info:" << endl;
	for(int tmp = 0; tmp < tmpKmerValueVec_filtered.size(); tmp++)
		taxo_assignment_ofs_assignInfo_detail << tmpKmerValueVec_filtered[tmp] << ",";
	taxo_assignment_ofs_assignInfo_detail << endl;

	vector<int> tmp_taxoReissuedIdVec_species, tmp_taxoReissuedIdVec_genus, tmp_taxoReissuedIdVec_family,
		tmp_taxoReissuedIdVec_order, tmp_taxoReissuedIdVec_class, tmp_taxoReissuedIdVec_phylum;
	taxoClassAssignmentInfo.convertGroupIdVec2reissuedIdVec(tmpKmerValueVec_filtered, tmp_taxoReissuedIdVec_species, 8, bacterialTaxoInfo);
	taxoClassAssignmentInfo.convertGroupIdVec2reissuedIdVec(tmpKmerValueVec_filtered, tmp_taxoReissuedIdVec_genus, 7, bacterialTaxoInfo);	
	taxoClassAssignmentInfo.convertGroupIdVec2reissuedIdVec(tmpKmerValueVec_filtered, tmp_taxoReissuedIdVec_family, 6, bacterialTaxoInfo);
	taxoClassAssignmentInfo.convertGroupIdVec2reissuedIdVec(tmpKmerValueVec_filtered, tmp_taxoReissuedIdVec_order, 5, bacterialTaxoInfo);
	taxoClassAssignmentInfo.convertGroupIdVec2reissuedIdVec(tmpKmerValueVec_filtered, tmp_taxoReissuedIdVec_class, 4, bacterialTaxoInfo);
	taxoClassAssignmentInfo.convertGroupIdVec2reissuedIdVec(tmpKmerValueVec_filtered, tmp_taxoReissuedIdVec_phylum, 3, bacterialTaxoInfo);		

	vector<int> tmpRawTaxoId_species, tmpRawTaxoId_genus, tmpRawTaxoId_family, 
		tmpRawTaxoId_order, tmpRawTaxoId_class, tmpRawTaxoId_phylum;
	bacterialTaxoInfo.get_taxoIdVec_from_taxoReissuedIdVec(tmpRawTaxoId_species, tmp_taxoReissuedIdVec_species, 8);
	bacterialTaxoInfo.get_taxoIdVec_from_taxoReissuedIdVec(tmpRawTaxoId_genus, tmp_taxoReissuedIdVec_genus, 7);
	bacterialTaxoInfo.get_taxoIdVec_from_taxoReissuedIdVec(tmpRawTaxoId_family, tmp_taxoReissuedIdVec_family, 6);
	bacterialTaxoInfo.get_taxoIdVec_from_taxoReissuedIdVec(tmpRawTaxoId_order, tmp_taxoReissuedIdVec_order, 5);
	bacterialTaxoInfo.get_taxoIdVec_from_taxoReissuedIdVec(tmpRawTaxoId_class, tmp_taxoReissuedIdVec_class, 4);
	bacterialTaxoInfo.get_taxoIdVec_from_taxoReissuedIdVec(tmpRawTaxoId_phylum, tmp_taxoReissuedIdVec_phylum, 3);

	taxo_assignment_ofs_assignInfo_detail << endl << "SPECIES-level taxo reissued id: " << "Discriminative: " << tmpDiscriminativeKmerCountVec[0] 
		<< " Repetitive: " << tmpRepetitiveKmerCountVec[0] << " Invalid: " << tmpInvalidKmerCountVec[0] << endl;
	for(int tmp = 0; tmp < tmpRawTaxoId_species.size(); tmp++)
		taxo_assignment_ofs_assignInfo_detail << tmpRawTaxoId_species[tmp] << ",";
	taxo_assignment_ofs_assignInfo_detail << endl << endl << "GENUS-level taxo reissued id: " << "Discriminative: " << tmpDiscriminativeKmerCountVec[1] 
		<< " Repetitive: " << tmpRepetitiveKmerCountVec[1] << " Invalid: " << tmpInvalidKmerCountVec[1] << endl;
	for(int tmp = 0; tmp < tmpRawTaxoId_genus.size(); tmp++)
		taxo_assignment_ofs_assignInfo_detail << tmpRawTaxoId_genus[tmp] << ",";
	taxo_assignment_ofs_assignInfo_detail << endl << endl << "FAMILY-level taxo reissued id: " << "Discriminative: " << tmpDiscriminativeKmerCountVec[2] 
	<< " Repetitive: " << tmpRepetitiveKmerCountVec[2] << " Invalid: " << tmpInvalidKmerCountVec[2] << endl;
	for(int tmp = 0; tmp < tmpRawTaxoId_family.size(); tmp++)
		taxo_assignment_ofs_assignInfo_detail << tmpRawTaxoId_family[tmp] << ",";
	taxo_assignment_ofs_assignInfo_detail << endl << endl << "ORDER-level taxo reissued id: " << "Discriminative: " << tmpDiscriminativeKmerCountVec[3] 
	<< " Repetitive: " << tmpRepetitiveKmerCountVec[3] << " Invalid: " << tmpInvalidKmerCountVec[3] << endl;
	for(int tmp = 0; tmp < tmpRawTaxoId_order.size(); tmp++)
		taxo_assignment_ofs_assignInfo_detail << tmpRawTaxoId_order[tmp] << ",";
	taxo_assignment_ofs_assignInfo_detail << endl << endl << "CLASS-level taxo reissued id: " << "Discriminative: " << tmpDiscriminativeKmerCountVec[4] 
	<< " Repetitive: " << tmpRepetitiveKmerCountVec[4] << " Invalid: " << tmpInvalidKmerCountVec[4] << endl;
	for(int tmp = 0; tmp < tmpRawTaxoId_class.size(); tmp++)
		taxo_assignment_ofs_assignInfo_detail << tmpRawTaxoId_class[tmp] << ",";
	taxo_assignment_ofs_assignInfo_detail << endl << endl << "PHYLUM-level taxo reissued id: " << "Discriminative: " << tmpDiscriminativeKmerCountVec[5] 
	<< " Repetitive: " << tmpRepetitiveKmerCountVec[5] << " Invalid: " << tmpInvalidKmerCountVec[5] << endl;
	for(int tmp = 0; tmp < tmpRawTaxoId_phylum.size(); tmp++)
		taxo_assignment_ofs_assignInfo_detail << tmpRawTaxoId_phylum[tmp] << ",";
	taxo_assignment_ofs_assignInfo_detail << endl;
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
		cout << "#7 inputSpeciesTaxoIdFile" << endl;
		cout << "#8 NCBIfullTaxoId2NameFile" << endl;
		cout << "#9 inputReadFile_1" << endl;
		cout << "(#10 inputReadFile_2)" << endl;
		exit(1);
	}
	int jump_length = 1;
	int normalRecordNum = 1000000;// 1;

	string inputLothelloNodeIndexStr = argv[1];
	string outputDirStr = argv[2];
	string kmer_length_str = argv[3];
    string threads_num_str = argv[4];
    string fa_or_fq_str = argv[5];
    string SE_or_PE_str = argv[6];
	string inputSpeciesTaxoIdFile = argv[7];
	string NCBIfullTaxoId2NameFile = argv[8];
	string inputReadFile_1 = argv[9];
	string inputReadFile_2;
	if(argc == 11)
		inputReadFile_2 = argv[10];

	//double assignment_windowSize_min = 0;
	double assignment_windowSizeSquareSum_min = 0;

	outputDirStr += "/";
  	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());
   	string log_file = outputDirStr + "log.txt";
   	ofstream log_ofs(log_file.c_str());
   	string stats_file = outputDirStr + "stats.txt";
   	ofstream stats_ofs(stats_file.c_str());   	

   	string command_line_str = "Command_Line:\n#1 inputLothelloNodeIndex: " + inputLothelloNodeIndexStr
   		+ "\n#2 outputFolder: " + outputDirStr + "\n#3 Kmer_length: " + kmer_length_str + "\n#4 threads_num: "
   		+ threads_num_str + "\n#5 fa_or_fq: " + fa_or_fq_str + "\n#6 SE_or_PE: " + SE_or_PE_str 
   		+ "\n#7 inputSpeciesTaxoIdFile: " + inputSpeciesTaxoIdFile + "\n#8 NCBIfullTaxoId2NameFile: "
   		+ NCBIfullTaxoId2NameFile + "\n#9 inputReadFile_1: " + inputReadFile_1;
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

	BacterialTaxo_Info bacterialTaxoInfo;
	bacterialTaxoInfo.initiate_bacterialTaxoFile_NCBIfullTaxoId2NameFile(inputSpeciesTaxoIdFile, NCBIfullTaxoId2NameFile);
	bacterialTaxoInfo.reissueTaxoIdName_all();
	bacterialTaxoInfo.initiate_lowerRankTaxoReissuedId_to_higherRankTaxoReissuedId_array();
	string outputDir_taxoInfo = outputDirStr + "taxoInfo/";
	bacterialTaxoInfo.print(outputDir_taxoInfo);

    int taxo_num_total = bacterialTaxoInfo.return_taxo_num_total();
    int repetitive_class_Id = taxo_num_total + 1;
    int alien_class_Id = 0;

    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "end of initiating bacterialTaxoInfo" << endl;

    bool fa_or_fq_bool, SE_or_PE_bool;
    get_Fa_or_Fq_SE_or_PE_bool(fa_or_fq_bool, SE_or_PE_bool, fa_or_fq_str, SE_or_PE_str, argc);
	bool InputAsFastq = (!fa_or_fq_bool);

    cout << endl << "[" << asctime(local) << "start to do query" << endl;
    string InputReadFile = argv[9];
	ifstream inputRead_ifs(InputReadFile.c_str());
	string InputReadFile_PE;
	if(!SE_or_PE_bool)
		InputReadFile_PE = argv[10];
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

	vector<int> assignment_taxoIdVec_species(normalRecordNum);
	vector<int> assignment_taxoIdVec_genus(normalRecordNum);
	vector<int> assignment_taxoIdVec_family(normalRecordNum);
	vector<int> assignment_taxoIdVec_order(normalRecordNum);
	vector<int> assignment_taxoIdVec_class(normalRecordNum);
	vector<int> assignment_taxoIdVec_phylum(normalRecordNum);					

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
	//ofstream taxo_assignment_ofs(output_taxo_assignment_file.c_str());
	FILE * taxo_assignment_pOut = fopen(output_taxo_assignment_file.c_str(), "w");

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
			
			#ifndef NO_FILTERING
			vector<valueT> tmpKmerValueVec_filtered;
			taxoClassAssignmentInfo.filter_KmerValueVec(tmpKmerValueVec_filtered, tmpKmerValueVec);
			#endif

			vector<MaximalSpecificKmerWindow_Info> tmpMSKWinfoVec;
			vector<int> tmpInvalidKmerCountVec;
			vector<int> tmpRepetitiveKmerCountVec;
			vector<int> tmpDiscriminativeKmerCountVec;
			
			#ifndef NO_FILTERING
			taxoClassAssignmentInfo.get_taxoId2windowVecMapVec_specific_allTaxoRank(
			 	tmpMSKWinfoVec, tmpKmerValueVec_filtered, kmer_length, bacterialTaxoInfo, 
			 	tmpInvalidKmerCountVec, tmpRepetitiveKmerCountVec, tmpDiscriminativeKmerCountVec);
			#else
			taxoClassAssignmentInfo.get_taxoId2windowVecMapVec_specific_allTaxoRank(
			 	tmpMSKWinfoVec, tmpKmerValueVec, kmer_length, bacterialTaxoInfo, 
			 	tmpInvalidKmerCountVec, tmpRepetitiveKmerCountVec, tmpDiscriminativeKmerCountVec);
			#endif

			#ifdef ASSIGN_INFO
			#ifndef NO_FILTERING
			print_groupId_rawTaxoIdVecVec(taxo_assignment_ofs_assignInfo_detail, bacterialTaxoInfo,
				taxoClassAssignmentInfo, SE_or_PE_bool,
				readName1Vec[tmpOpenMP], readSeq1Vec[tmpOpenMP], readQualSeq1Vec[tmpOpenMP],
				readName2Vec[tmpOpenMP], readSeq2Vec[tmpOpenMP], readQualSeq2Vec[tmpOpenMP], 
				tmpKmerValueVec, tmpKmerValueVec_filtered,
				tmpDiscriminativeKmerCountVec, tmpRepetitiveKmerCountVec, tmpInvalidKmerCountVec);
			#else
			print_groupId_rawTaxoIdVecVec(taxo_assignment_ofs_assignInfo_detail, bacterialTaxoInfo,
				taxoClassAssignmentInfo, SE_or_PE_bool,
				readName1Vec[tmpOpenMP], readSeq1Vec[tmpOpenMP], readQualSeq1Vec[tmpOpenMP],
				readName2Vec[tmpOpenMP], readSeq2Vec[tmpOpenMP], readQualSeq2Vec[tmpOpenMP], 
				tmpKmerValueVec, tmpKmerValueVec,
				tmpDiscriminativeKmerCountVec, tmpRepetitiveKmerCountVec, tmpInvalidKmerCountVec);			
			#endif
			#endif

			int tmp_assignment_taxoId_species, tmp_assignment_taxoId_genus, tmp_assignment_taxoId_family,
				tmp_assignment_taxoId_order, tmp_assignment_taxoId_class, tmp_assignment_taxoId_phylum; 

			taxoClassAssignmentInfo.get_assignment_taxoId_allTaxoRank(
				#ifdef ASSIGN_INFO
				taxo_assignment_ofs_assignInfo_detail,
				#endif
				bacterialTaxoInfo, tmpMSKWinfoVec, 
				//assignment_windowSize_min, 
				assignment_windowSizeSquareSum_min,
				tmpInvalidKmerCountVec, tmpRepetitiveKmerCountVec, tmpDiscriminativeKmerCountVec,
				tmp_assignment_taxoId_species, tmp_assignment_taxoId_genus, tmp_assignment_taxoId_family,
				tmp_assignment_taxoId_order, tmp_assignment_taxoId_class, tmp_assignment_taxoId_phylum);

			assignment_taxoIdVec_species[tmpOpenMP] = tmp_assignment_taxoId_species;
			assignment_taxoIdVec_genus[tmpOpenMP] = tmp_assignment_taxoId_genus;
			assignment_taxoIdVec_family[tmpOpenMP] = tmp_assignment_taxoId_family;
			assignment_taxoIdVec_order[tmpOpenMP] = tmp_assignment_taxoId_order;
			assignment_taxoIdVec_class[tmpOpenMP] = tmp_assignment_taxoId_class;
			assignment_taxoIdVec_phylum[tmpOpenMP] = tmp_assignment_taxoId_phylum;
		}
		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		cout << endl << "[" << asctime(local) << "start to output ... turn: " << tmpTurn+1 << endl;
		log_ofs << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		log_ofs << endl << "[" << asctime(local) << "start to output ... turn: " << tmpTurn+1 << endl;			
		// for(int tmp = 0; tmp < realRecordNum; tmp++)
		// 	taxo_assignment_ofs << readName1Vec[tmp].c_str() << "\t" 
		// 		<< assignment_taxoIdVec_species[tmp] << "\t" << assignment_taxoIdVec_genus[tmp] << "\t"
		// 		<< assignment_taxoIdVec_family[tmp] << "\t" << assignment_taxoIdVec_order[tmp] << "\t"
		// 		<< assignment_taxoIdVec_class[tmp] << "\t" << assignment_taxoIdVec_phylum[tmp] << endl;
		for(int tmp = 0; tmp < realRecordNum; tmp++)
			fprintf(taxo_assignment_pOut, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", readName1Vec[tmp].c_str(),
				assignment_taxoIdVec_species[tmp], assignment_taxoIdVec_genus[tmp], assignment_taxoIdVec_family[tmp],
				assignment_taxoIdVec_order[tmp], assignment_taxoIdVec_class[tmp], assignment_taxoIdVec_phylum[tmp]);
		#ifdef ASSIGN_INFO
		for(int tmp = 0; tmp < realRecordNum; tmp++)
			taxo_assignment_ofs_assignInfo << readName1Vec[tmp].c_str() << "\t" 
				<< assignment_taxoIdVec_species[tmp] << "\t" << assignment_taxoIdVec_genus[tmp] << "\t"
				<< assignment_taxoIdVec_family[tmp] << "\t" << assignment_taxoIdVec_order[tmp] << "\t"
				<< assignment_taxoIdVec_class[tmp] << "\t" << assignment_taxoIdVec_phylum[tmp] << endl;
		#endif
	}

	#ifdef ASSIGN_INFO
	taxo_assignment_ofs_assignInfo_detail.close();
	taxo_assignment_ofs_assignInfo.close();
	#endif
	//taxo_assignment_ofs.close();
	fclose(taxo_assignment_pOut);

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