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
	if((argc != 14)&&(argc != 15))
	{
		cout << "#1 inputLothelloNodeIndex" << endl;
		cout << "#2 outputFolder" << endl; 
		cout << "#3 Kmer_length" << endl;
		cout << "#4 threads_num " << endl;
		cout << "#5 genome_num" << endl;
		cout << "#6 fa_or_fq" << endl; 
		cout << "#7 SE_or_PE" << endl;
		cout << "#8 DiscrimitiveKmerProportionMin_Best_SecondBest_Diff" << endl;
		cout << "#9 reissuedId2oriGenomeIdFile" << endl;
		cout << "#10 oriGenomeid2taxoIdFile" << endl;
		cout << "#11 NCBIfullTaxoId2NameFile" << endl;
		cout << "#12 Rank_Level(Phylum/Genus/Species)" << endl;
		cout << "#13 inputReadFile_1" << endl;
		cout << "(#14 inputReadFile_2)" << endl;
		exit(1);
	}
	int jump_length = 1;
	int normalRecordNum = 1000000;//0;
	int rank;
	string rank_str = argv[12];
	if(rank_str == "Phylum")
		rank = 3;
	else if(rank_str == "Genus")
		rank = 7;
	else if(rank_str == "Species")
		rank = 8;
	else
	{
		cout << "invalid rank: " << rank_str << endl;
		exit(1);
	}

	string outputDirStr = argv[2];
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
    string kmer_length_str = argv[3];
    int kmer_length = atoi(kmer_length_str.c_str());
    string threads_num_str = argv[4];
    int threads_num = atoi(threads_num_str.c_str());
    string genome_num_str = argv[5];
    int genome_num = atoi(genome_num_str.c_str());
    int repetitiveClassId = genome_num + 1;
    int alienClassId = 0;

    int splitbit = 6;
    IOHelper<keyT,valueT> *helper;
    helper = new ConstantLengthKmerHelper<keyT,valueT>(kmer_length,splitbit);
    MulOth<keyT> * moth;
    moth = new MulOth<keyT>(argv[1], helper);    

    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "end of loading index" << endl;

	string reissuedId2oriGenomeIdFile = argv[9];
	string oriGenomeid2taxoIdFile = argv[10];
	string NCBIfullTaxoId2NameFile = argv[11];

	ReissuedGenomeID2TaxoID_Info reissuedGenomeId2TaxoIdInfo;
	reissuedGenomeId2TaxoIdInfo.initiate_reissuedId2oriGenomeIdFile_oriGenomeId2taxoIdFile_NCBIfullTaxoId2NameFile(
		reissuedId2oriGenomeIdFile, oriGenomeid2taxoIdFile, NCBIfullTaxoId2NameFile);
	reissuedGenomeId2TaxoIdInfo.generate_reissuedTaxoIdVec_taxoIdNamePairVec();
	string output_reissuedGenomeID2TaxoIDName = outputDirStr + "reissuedId_to_taxoIdName.txt";
	string output_species_info = outputDirStr + "species_info.txt";
	string output_genus_info = outputDirStr + "genus_info.txt";
	string output_phylum_info = outputDirStr + "phylum_info.txt";
	reissuedGenomeId2TaxoIdInfo.print_reissuedGenomeId2taxoIdName(output_reissuedGenomeID2TaxoIDName,
		output_species_info, output_genus_info, output_phylum_info);
	int reissuedGenomeIDnum = reissuedGenomeId2TaxoIdInfo.return_genome_num();
	int species_num = reissuedGenomeId2TaxoIdInfo.return_species_num();
	int genus_num = reissuedGenomeId2TaxoIdInfo.return_genus_num();
	int phylum_num = reissuedGenomeId2TaxoIdInfo.return_phylum_num();
	log_ofs << "Genome #:\t" << reissuedGenomeIDnum << endl;
	log_ofs << "Species #:\t" << species_num << endl;
	log_ofs << "Genus #:\t" << genus_num << endl;
	log_ofs << "Phylum #:\t" << phylum_num << endl;

    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "end of initiating reissuedGenomeId2TaxoIdInfo" << endl;

	vector< vector<unsigned long long> > speciesCountVecVec_thread;
	vector< vector<unsigned long long> > genusCountVecVec_thread;
	vector< vector<unsigned long long> > phylumCountVecVec_thread;
	vector<unsigned long long> repetitiveCountVec_thread;
	for(int tmp = 0; tmp < threads_num; tmp++)
	{
		vector<unsigned long long> tmp_speciesCountVec_thread;
		for(int tmp2 = 0; tmp2 < species_num; tmp2++)
			tmp_speciesCountVec_thread.push_back(0);
		speciesCountVecVec_thread.push_back(tmp_speciesCountVec_thread);
		vector<unsigned long long> tmp_genusCountVec_thread;
		for(int tmp2 = 0; tmp2 < genus_num; tmp2++)
			tmp_genusCountVec_thread.push_back(0);
		genusCountVecVec_thread.push_back(tmp_genusCountVec_thread);
		vector<unsigned long long> tmp_phylumCountVec_thread;
		for(int tmp2 = 0; tmp2 < phylum_num; tmp2++)
			tmp_phylumCountVec_thread.push_back(0);
		phylumCountVecVec_thread.push_back(tmp_phylumCountVec_thread);
		repetitiveCountVec_thread.push_back(0);
	}

	if(genome_num != reissuedGenomeIDnum)
	{
		cout << "genome_num set in command line is not consistent with reissuedGenomeIDnum" << endl;
		cout << "genome_num: " << genome_num << endl;
		cout << "reissuedGenomeIDnum: " << reissuedGenomeIDnum << endl;
		exit(1);
	}	
	string DiscrimitiveKmerProportionMinStr_best_secondBest_diff = argv[8];
	double DiscrimitiveKmerProportionMin_best_secondBest_diff = atof(DiscrimitiveKmerProportionMinStr_best_secondBest_diff.c_str());
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

    if((SE_or_PE_bool && (argc != 14)) || ((!SE_or_PE_bool) && (argc != 15)))
    {
    	cout << "inconsistent SE_or_PE_bool and argc" << endl;
    	cout << "SE_or_PE_bool: " << SE_or_PE_bool << endl;
    	cout << "argc: " << argc << endl;
    	exit(1);
    }
    
    cout << endl << "[" << asctime(local) << "start to do query" << endl;
    //string output_query_file = outputDirStr + "query.txt";
    //ofstream query_ofs(output_query_file.c_str());
    string InputReadFile = argv[13];
	ifstream inputRead_ifs(InputReadFile.c_str());
	string InputReadFile_PE;
	if(!SE_or_PE_bool)
		InputReadFile_PE = argv[14];
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
	vector<string> readAssignedTaxoClassNameVec(normalRecordNum);
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
	string output_genome_assignment_file_assignInfo = outputDirStr + "taxo_assignment.txt.assignInfo";
	ofstream genome_assignment_ofs_assignInfo(output_genome_assignment_file_assignInfo.c_str());
	string output_genome_assignment_file_assignInfo_detail = outputDirStr + "taxo_assignment.txt.assignInfo.detail";
	ofstream genome_assignment_ofs_assignInfo_detail(output_genome_assignment_file_assignInfo_detail.c_str());
	#endif
	string output_genome_assignment_file = outputDirStr + "taxo_assignment.txt";
	ofstream genome_assignment_ofs(output_genome_assignment_file.c_str());

	TaxoClassAssignment_Info taxoClassAssignmentInfo;
	taxoClassAssignmentInfo.initiate(reissuedGenomeId2TaxoIdInfo, rank);
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

			#ifdef SUMMARIZE_TAXOCLASS_INSTEADOF_GENOME_COUNT
			vector< pair<int,int> > tmpTaxoClassIndexCountPairVec;
			#else
			vector< pair<int,int> > tmpGenomeIdCountPairVec;
			#endif
			if(SE_or_PE_bool)
			{
				#ifdef SUMMARIZE_TAXOCLASS_INSTEADOF_GENOME_COUNT
				taxoClassAssignmentInfo.get_taxoClassIndexCountPairVec_from_rawKmerValueVec_SE(tmpKmerValueVec_1,
					tmpInvalidKmerCount, tmpRepetitiveKmerCount, tmpDiscriminativeKmerCount, tmpTaxoClassIndexCountPairVec);
				#else
				taxoClassAssignmentInfo.get_genomeIdCountPairVec_from_rawKmerValueVec_SE(tmpKmerValueVec_1,
					tmpInvalidKmerCount, tmpRepetitiveKmerCount, tmpDiscriminativeKmerCount, tmpGenomeIdCountPairVec);				
				#endif
			}
			else
			{
				QuerySeq_Info tmpQuerySeqInfo_2;
				tmpQuerySeqInfo_2.initiate(readName2Vec[tmpOpenMP], readSeq2Vec[tmpOpenMP]);
				vector<valueT> tmpKmerValueVec_2;
				tmpQuerySeqInfo_2.querySeq_returnKmerValueVec(tmpKmerValueVec_2, moth, kmer_length, jump_length);
				#ifdef SUMMARIZE_TAXOCLASS_INSTEADOF_GENOME_COUNT
				taxoClassAssignmentInfo.get_taxoClassIndexCountPairVec_from_rawKmerValueVec_PE(tmpKmerValueVec_1, tmpKmerValueVec_2,
					tmpInvalidKmerCount, tmpRepetitiveKmerCount, tmpDiscriminativeKmerCount, tmpTaxoClassIndexCountPairVec);
				#else
				taxoClassAssignmentInfo.get_genomeIdCountPairVec_from_rawKmerValueVec_PE(tmpKmerValueVec_1, tmpKmerValueVec_2,
					tmpInvalidKmerCount, tmpRepetitiveKmerCount, tmpDiscriminativeKmerCount, tmpGenomeIdCountPairVec);
				#endif
				#ifdef ASSIGN_INFO
				genome_assignment_ofs_assignInfo_detail << endl << readName1Vec[tmpOpenMP] << endl << readSeq1Vec[tmpOpenMP] 
					<< endl << readQualSeq1Vec[tmpOpenMP] << endl << readName2Vec[tmpOpenMP] << endl << readSeq2Vec[tmpOpenMP] 
					<< endl << readQualSeq2Vec[tmpOpenMP] << endl;
				for(int tmp = 0; tmp < tmpKmerValueVec_1.size(); tmp++)
					genome_assignment_ofs_assignInfo_detail << tmpKmerValueVec_1[tmp] << ",";
				for(int tmp = 0; tmp < tmpKmerValueVec_2.size(); tmp++)
					genome_assignment_ofs_assignInfo_detail << tmpKmerValueVec_2[tmp] << ",";
				genome_assignment_ofs_assignInfo_detail << endl;			
				for(int tmp = 0; tmp < tmpKmerValueVec_1.size(); tmp++)
					genome_assignment_ofs_assignInfo_detail << reissuedGenomeId2TaxoIdInfo.return_reissuedTaxoId_from_genomeId(tmpKmerValueVec_1[tmp], rank) << ",";
				for(int tmp = 0; tmp < tmpKmerValueVec_2.size(); tmp++)
					genome_assignment_ofs_assignInfo_detail << reissuedGenomeId2TaxoIdInfo.return_reissuedTaxoId_from_genomeId(tmpKmerValueVec_2[tmp], rank) << ",";
				#endif
				// debug ....
			}
			#ifdef SUMMARIZE_TAXOCLASS_INSTEADOF_GENOME_COUNT
			int tmpReissuedId = taxoClassAssignmentInfo.determineAssignedTaxoClassIndex_from_taxoClassIndexCountPairVec(
				DiscrimitiveKmerProportionMin_best_secondBest_diff, 2, tmpInvalidKmerCount, tmpRepetitiveKmerCount, 
				tmpDiscriminativeKmerCount, tmpTaxoClassIndexCountPairVec, 
				tmpKmerIndex_best, tmpKmerCount_best, tmpKmerIndex_secondBest, tmpKmerCount_secondBest);
			#else
			int tmpReissuedId = taxoClassAssignmentInfo.determineAssignedTaxoClassIndex_from_genomeIdCountPairVec(
				DiscrimitiveKmerProportionMin_best_secondBest_diff, 2, tmpInvalidKmerCount, tmpRepetitiveKmerCount,
				tmpDiscriminativeKmerCount, tmpGenomeIdCountPairVec, reissuedGenomeId2TaxoIdInfo,
				tmpKmerIndex_best, tmpKmerCount_best, tmpKmerIndex_secondBest, tmpKmerCount_secondBest);
			#endif

			//readAssignedTaxoClassVec[tmpOpenMP] = tmpReissuedId;
			#ifdef ASSIGN_INFO
			invalidKmerCountVec[tmpOpenMP] = tmpInvalidKmerCount;
			repetitiveKmerCountVec[tmpOpenMP] = tmpRepetitiveKmerCount;
			discriminativeKmerCountVec[tmpOpenMP] = tmpDiscriminativeKmerCount;
			KmerCountVec_best[tmpOpenMP] = tmpKmerCount_best;
			KmerCountVec_secondBest[tmpOpenMP] = tmpKmerCount_secondBest;
			KmerIndexVec_best[tmpOpenMP] = tmpKmerIndex_best;
			KmerIndexVec_secondBest[tmpOpenMP] = tmpKmerIndex_secondBest;

			genome_assignment_ofs_assignInfo_detail << endl;
			genome_assignment_ofs_assignInfo_detail << KmerIndexVec_best[tmpOpenMP] << "-" 
				<< KmerCountVec_best[tmpOpenMP] << "\t" << KmerIndexVec_secondBest[tmpOpenMP] << "-" 
				<< KmerCountVec_secondBest[tmpOpenMP] << "\t" << discriminativeKmerCountVec[tmpOpenMP] 
				<< "," << repetitiveKmerCountVec[tmpOpenMP] << "," << invalidKmerCountVec[tmpOpenMP] << endl;
			#endif
			if(tmpReissuedId < 0)
				repetitiveCountVec_thread[threadNO] ++;
			int tmpTaxoClassId;
			string tmpTaxoClassName;
			if(rank == 3){
				if(tmpReissuedId >= 0)
					(phylumCountVecVec_thread[threadNO])[tmpReissuedId] ++;				
				reissuedGenomeId2TaxoIdInfo.return_phylum_Id_name_from_reissuedId(tmpReissuedId, tmpTaxoClassId, tmpTaxoClassName);
			}
			else if(rank == 7){
				if(tmpReissuedId >= 0)
					(genusCountVecVec_thread[threadNO])[tmpReissuedId] ++;					
				reissuedGenomeId2TaxoIdInfo.return_genus_Id_name_from_reissuedId(tmpReissuedId, tmpTaxoClassId, tmpTaxoClassName);
			}
			else if(rank == 8){
				if(tmpReissuedId >= 0)
					(speciesCountVecVec_thread[threadNO])[tmpReissuedId] ++;	
				reissuedGenomeId2TaxoIdInfo.return_species_Id_name_from_reissuedId(tmpReissuedId, tmpTaxoClassId, tmpTaxoClassName);			
			}
			else{
				cout << "Not application for rank: " << rank << endl;
				exit(1);
			}
			readAssignedTaxoClassNameVec[tmpOpenMP] = tmpTaxoClassName;
		}
		nowtime = time(NULL);
		local = localtime(&nowtime);
		cout << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		cout << endl << "[" << asctime(local) << "start to output ... turn: " << tmpTurn+1 << endl;
		log_ofs << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		log_ofs << endl << "[" << asctime(local) << "start to output ... turn: " << tmpTurn+1 << endl;			
		//for(int tmp = 0; tmp < realRecordNum; tmp++)
		//	genome_assignment_ofs << readName1Vec[tmp].c_str() << "\t" << readAssignedTaxoClassNameVec[tmp] << endl;
		#ifdef ASSIGN_INFO
		genome_assignment_ofs_assignInfo_detail.close();
		for(int tmp = 0; tmp < realRecordNum; tmp++)
			genome_assignment_ofs_assignInfo << readName1Vec[tmp].c_str() << "\t" << readAssignedTaxoClassNameVec[tmp]
				<< "\t" << KmerIndexVec_best[tmp] << "-" << KmerCountVec_best[tmp] << "\t"
				<< KmerIndexVec_secondBest[tmp] << "-" << KmerCountVec_secondBest[tmp] << "\t"
				<< discriminativeKmerCountVec[tmp] << "," << repetitiveKmerCountVec[tmp] << "," << invalidKmerCountVec[tmp] << endl;		
		#endif
	}

	#ifdef ASSIGN_INFO
	genome_assignment_ofs_assignInfo.close();
	#endif
	genome_assignment_ofs.close();
	
	int tmpTaxoId; 
	string tmpTaxoName;
	if(rank == 3)
	{
		vector<unsigned long long> phylumCountVec;
		for(int tmp = 0; tmp < phylum_num; tmp++)
			phylumCountVec.push_back(0);
		for(int tmp = 0; tmp < threads_num; tmp++){
			for(int tmp2 = 0; tmp2 < phylum_num; tmp2++)
				phylumCountVec[tmp2] += (phylumCountVecVec_thread[tmp])[tmp2];
		}
		for(int tmpReissuedId = 0; tmpReissuedId < phylum_num; tmpReissuedId++){
			reissuedGenomeId2TaxoIdInfo.return_phylum_Id_name_from_reissuedId(tmpReissuedId, tmpTaxoId, tmpTaxoName);
			stats_ofs << tmpReissuedId << "\t" << tmpTaxoId << "\t" << tmpTaxoName << "\t" << phylumCountVec[tmpReissuedId] << endl;		
		}
	}
	else if(rank == 7)
	{	
		vector<unsigned long long> genusCountVec;
		for(int tmp = 0; tmp < genus_num; tmp++)
			genusCountVec.push_back(0);
		for(int tmp = 0; tmp < threads_num; tmp++){
			for(int tmp2 = 0; tmp2 < genus_num; tmp2++)
				genusCountVec[tmp2] += (genusCountVecVec_thread[tmp])[tmp2];
		}
		for(int tmpReissuedId = 0; tmpReissuedId < genus_num; tmpReissuedId++){
			reissuedGenomeId2TaxoIdInfo.return_genus_Id_name_from_reissuedId(tmpReissuedId, tmpTaxoId, tmpTaxoName);
			stats_ofs << tmpReissuedId << "\t" << tmpTaxoId << "\t" << tmpTaxoName << "\t" << genusCountVec[tmpReissuedId] << endl;
		}
	}
	else if(rank == 8)
	{		
		vector<unsigned long long> speciesCountVec;
		for(int tmp = 0; tmp < species_num; tmp++)
			speciesCountVec.push_back(0);
		for(int tmp = 0; tmp < threads_num; tmp++){
			for(int tmp2 = 0; tmp2 < species_num; tmp2++)
				speciesCountVec[tmp2] += (speciesCountVecVec_thread[tmp])[tmp2];
		}
		for(int tmpReissuedId = 0; tmpReissuedId < species_num; tmpReissuedId++){
			reissuedGenomeId2TaxoIdInfo.return_species_Id_name_from_reissuedId(tmpReissuedId, tmpTaxoId, tmpTaxoName);
			stats_ofs << tmpReissuedId << "\t" << tmpTaxoId << "\t" << tmpTaxoName << "\t" << speciesCountVec[tmpReissuedId] << endl;	
		}
	}
	else
	{
		cout << "Not applicable for rank: " << rank << endl;
		exit(1);
	}
	
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