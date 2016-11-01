#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <bitset>
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
#include "../../mps3Lib/read_block_test.h"
#include "../../mps3Lib/otherFunc.h"
#include "../../mps3Lib/index_info.h"
#include "../general/NCBIfullTaxoID2Name_info.h"
#include "../general/bacterialTaxo_info.h"
using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 7)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputArtBin" << endl;
		// /scratch/lcph222/Xinan/lothelloClassifier/metagenomics/total/total_taxo_info/speciesFaFile2taxoInfo.txt
		cout << "#2 speciesFaFile2taxoInfo" << endl;
		// /scratch/lcph222/Xinan/lothelloClassifier/metagenomics/total/total_taxo_info/names.dmp.scientific
		cout << "#3 NCBIfullTaxoId2NameFile" << endl;
		cout << "#4 inputSpeciesCountFile" << endl;
		cout << "#5 outputDir" << endl;
		cout << "#6 SE_MiSeq_or_PE_HiSeq" << endl;
		exit(1);
	}
	string inputArtBin = argv[1];
	string speciesFaFile2taxoInfo = argv[2];
	string NCBIfullTaxoId2NameFile = argv[3];
	string inputSpeciesCountFile = argv[4];
	string outputDir = argv[5];
	inputArtBin += "/";
	outputDir += "/";
	string mkdir_output = "mkdir " + outputDir;
	system(mkdir_output.c_str());
	string log_file = outputDir + "log.txt";
	ofstream log_ofs(log_file.c_str());
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    log_ofs << endl << "[" << asctime(local) << "start to simulate reads ..." << endl;
	
	log_ofs << "Command Line:" << endl;
	for(int tmp = 0; tmp < argc; tmp++)
		log_ofs << "#" << tmp << "\t" << argv[tmp] << endl; 

	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    log_ofs << endl << "[" << asctime(local) << "start to initiate options for art" << endl;
    string SE_MiSeq_or_PE_HiSeq_str = argv[6];
    bool SE_MiSeq_or_PE_HiSeq_bool;
    if(SE_MiSeq_or_PE_HiSeq_str == "SE_MiSeq")
    	SE_MiSeq_or_PE_HiSeq_bool = true;
    else if(SE_MiSeq_or_PE_HiSeq_str == "PE_HiSeq")
    	SE_MiSeq_or_PE_HiSeq_bool = false;
    else
    {
    	cout << "Invalid parameter for SE_MiSeq_or_PE_HiSeq_str: " << SE_MiSeq_or_PE_HiSeq_str << endl;
    	exit(1);
    }
	string opt_sequencing_platform, opt_length;
	if(SE_MiSeq_or_PE_HiSeq_bool)
	{	
		opt_sequencing_platform = "MSv3";
		opt_length = "250";
	}
	else
	{
		opt_sequencing_platform = "HS25";
		opt_length = "100 -m 200 -s 10 -p";
	}

	////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////
	vector<string> speciesFaPathVec;
	vector<int> speciesIdVec;
	vector<int> genusIdVec;
	vector<int> familyIdVec;
	vector<int> orderIdVec;
	vector<int> classIdVec;
	vector<int> phylumIdVec;
	vector<int> speciesCountVec;

	string inputSpeciesTaxoIdFile = outputDir + "speciesId2taxoInfo.txt";
	string inputSpeciesFaListFile = outputDir + "speciesFaList.txt";
	string cut_speciesFaFile2taxoInfoFile_2_speciesId2taxoInfoFile_cmd
		= "cut -f 2,3,4,5,6,7 " + speciesFaFile2taxoInfo + " > " + inputSpeciesTaxoIdFile;
	system(cut_speciesFaFile2taxoInfoFile_2_speciesId2taxoInfoFile_cmd.c_str());
	string cut_speciesFaFile2taxoInfoFile_2_speciesFaListFile_cmd
		= "cut -f 1 " + speciesFaFile2taxoInfo + " > " + inputSpeciesFaListFile;
	system(cut_speciesFaFile2taxoInfoFile_2_speciesFaListFile_cmd.c_str());

	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    log_ofs << endl << "[" << asctime(local) << "start to initiate speciesFaPathVec" << endl;
    cout << endl << "[" << asctime(local) << "start to initiate speciesFaPathVec" << endl;
	ifstream speciesFaList_ifs(inputSpeciesFaListFile.c_str());
	while(!speciesFaList_ifs.eof())
	{
		string tmpStr;
		getline(speciesFaList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		speciesFaPathVec.push_back(tmpStr);
	}
	speciesFaList_ifs.close();

	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    log_ofs << endl << "[" << asctime(local) << "start to initiate bacterialTaxoInfo" << endl;
    cout << endl << "[" << asctime(local) << "start to initiate bacterialTaxoInfo" << endl;
	BacterialTaxo_Info bacterialTaxoInfo;
	bacterialTaxoInfo.initiate_bacterialTaxoFile_NCBIfullTaxoId2NameFile(
		inputSpeciesTaxoIdFile, NCBIfullTaxoId2NameFile);
	bacterialTaxoInfo.reissueTaxoIdName_all();
	string taxo_dir = outputDir + "taxo_info";
	bacterialTaxoInfo.print(taxo_dir);
	bacterialTaxoInfo.print2Vec_species_genus_family_order_class_phylum(
		speciesIdVec, genusIdVec, familyIdVec, orderIdVec, classIdVec, phylumIdVec);
	// int species_num = bacterialTaxoInfo.return_taxo_num(8);
	// int genus_num = bacterialTaxoInfo.return_taxo_num(7);
	// int family_num = bacterialTaxoInfo.return_taxo_num(6);
	// int order_num = bacterialTaxoInfo.return_taxo_num(5);
	// int class_num = bacterialTaxoInfo.return_taxo_num(4);
	// int phylum_num = bacterialTaxoInfo.return_taxo_num(3);
	// log_ofs << endl << "species_num: " << species_num << endl << "genus_num: " << genus_num << endl
	// 	<< "family_num: " << family_num << endl << "order_num: " << order_num << endl
	// 	<< "class_num: " << class_num << endl << "phylum_num: " << phylum_num << endl;
	// cout << endl << "species_num: " << species_num << endl << "genus_num: " << genus_num << endl
	// 	<< "family_num: " << family_num << endl << "order_num: " << order_num << endl
	// 	<< "class_num: " << class_num << endl << "phylum_num: " << phylum_num << endl;		

	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    log_ofs << endl << "[" << asctime(local) << "start to generate speciesCountVec" << endl;
    cout << endl << "[" << asctime(local) << "start to generate speciesCountVec" << endl;
    for(int tmp = 0; tmp < speciesIdVec.size(); tmp++)
    	speciesCountVec.push_back(0);
    ifstream speciesCount_ifs(inputSpeciesCountFile.c_str());
    while(!speciesCount_ifs.eof())
    {
    	string tmpStr;
    	getline(speciesCount_ifs, tmpStr);
    	if(tmpStr == "")
    		break;
    	int tabLoc = tmpStr.find("\t");
    	string tmpSpeciesIdStr = tmpStr.substr(0, tabLoc);
    	string tmpCountStr = tmpStr.substr(tabLoc + 1);
    	int tmpSpeciesId = atoi(tmpSpeciesIdStr.c_str());
    	int tmpCount = atoi(tmpCountStr.c_str());
    	cout << "tmpSpeciesId: " << tmpSpeciesId << endl;
    	cout << "tmpCount: " << tmpCount << endl;
    	for(int tmp = 0; tmp < speciesIdVec.size(); tmp++)
    	{
    		if(speciesIdVec[tmp] == tmpSpeciesId)
    		{
    			cout << "speciesIdVec[tmp] == tmpSpeciesId" << endl;
    			speciesCountVec[tmp] = tmpCount;
    		}
    	}
    }
    speciesCount_ifs.close();

	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	string outputDir_rawSimulatedFa = outputDir + "raw_simulated_fa/";
	string mkdir_rawSimulatedFa = "mkdir " + outputDir_rawSimulatedFa;
	system(mkdir_rawSimulatedFa.c_str());
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    log_ofs << endl << "[" << asctime(local) << "start to simulate reads for each fa" << endl;
    cout << endl << "[" << asctime(local) << "start to simulate reads for each fa" << endl;
    for(int tmp = 0; tmp < speciesIdVec.size(); tmp++)
    {
		string tmpSpeciesFa = speciesFaPathVec[tmp];
		int tmpSpeciesId = speciesIdVec[tmp];
		int tmpGenusId = genusIdVec[tmp];
		int tmpFamilyId = familyIdVec[tmp];
		int tmpOrderId = orderIdVec[tmp];
		int tmpClassId = classIdVec[tmp];
		int tmpPhylumId = phylumIdVec[tmp];
		int tmpSpeciesCount = speciesCountVec[tmp];
		if(tmpSpeciesCount == 0)
		{
			log_ofs << "tmpSpeciesCount == 0" << endl;
			log_ofs << "tmpSpeciesFa: " << tmpSpeciesFa << endl;
			log_ofs << "tmpTaxoInfo: " << tmpSpeciesId << "_" << tmpGenusId << "_" << tmpFamilyId
				<< "_" << tmpOrderId << "_" << tmpClassId << "_" << tmpPhylumId << endl;
			continue;
		}
		else
		{
			log_ofs << "tmpSpeciesCount: " << tmpSpeciesCount << endl; 
			log_ofs << "tmpSpeciesFa: " << tmpSpeciesFa << endl;
			log_ofs << "tmpTaxoInfo: " << tmpSpeciesId << "_" << tmpGenusId << "_" << tmpFamilyId
				<< "_" << tmpOrderId << "_" << tmpClassId << "_" << tmpPhylumId << endl;			
		}
		string tmpTargetFa = outputDir_rawSimulatedFa + int_to_str(tmpSpeciesId);
		string tmp_cmd_simulate_read = inputArtBin + "art_illumina -ss " + opt_sequencing_platform
			+ " -na -l " + opt_length + " -i " + tmpSpeciesFa + " -c " + int_to_str(tmpSpeciesCount)
			+ " -o " + tmpTargetFa;
		if(!SE_MiSeq_or_PE_HiSeq_bool)
			tmp_cmd_simulate_read += ".";
		cout << "tmp_cmd_simulate_read: " << tmp_cmd_simulate_read << endl;
		system(tmp_cmd_simulate_read.c_str());
    	
    }

	log_ofs.close();
	return 0;
}