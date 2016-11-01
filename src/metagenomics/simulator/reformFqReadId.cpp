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
		// /scratch/lcph222/Xinan/lothelloClassifier/metagenomics/total/total_taxo_info/speciesFaFile2taxoInfo.txt
		cout << "#1 speciesFaFile2taxoInfo" << endl;
		// /scratch/lcph222/Xinan/lothelloClassifier/metagenomics/total/total_taxo_info/names.dmp.scientific
		cout << "#2 NCBIfullTaxoId2NameFile" << endl;
		cout << "#3 inputSpeciesCountFile" << endl;
		cout << "#4 inputFqDir (+=.fq for SE reads; += .1/2.fq for PE reads)" << endl;		
		cout << "#5 outputDir" << endl;
		cout << "#6 SE_or_PE" << endl;
		exit(1);
	}
	string speciesFaFile2taxoInfo = argv[1];
	string NCBIfullTaxoId2NameFile = argv[2];
	string inputSpeciesCountFile = argv[3];
	string inputFqDir = argv[4];
	string outputDir = argv[5];
    string SE_or_PE_str = argv[6];

	outputDir += "/";
	string mkdir_output = "mkdir " + outputDir;
	system(mkdir_output.c_str());
	string log_file = outputDir + "log.txt";
	ofstream log_ofs(log_file.c_str());
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    log_ofs << endl << "[" << asctime(local) << "start to reform fq read id ..." << endl;
	log_ofs << "Command Line:" << endl;
	for(int tmp = 0; tmp < argc; tmp++)
		log_ofs << "#" << tmp << "\t" << argv[tmp] << endl;

    bool SE_or_PE_bool;
    if(SE_or_PE_str == "SE")
    	SE_or_PE_bool = true;
    else if(SE_or_PE_str == "PE")
    	SE_or_PE_bool = false;
    else
    {
    	cout << "Invalid parameter for SE_or_PE_str: " << SE_or_PE_str << endl;
    	exit(1);
    }
	////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////
	//vector<string> speciesFaPathVec;
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
    inputFqDir += "/";
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    log_ofs << endl << "[" << asctime(local) << "start to reform raw read id and merge 2 single fq file" << endl;
    cout << endl << "[" << asctime(local) << "start to reform raw read id and merge 2 single fq file" << endl;
    vector<int> speciesCountVec_inFq;
   	for(int tmp = 0; tmp < speciesIdVec.size(); tmp++)
   		speciesCountVec_inFq.push_back(0);
	string outputSpeciesCountFileInFqFile = outputDir + "speciesCount_inFqFile.txt";
	ofstream speciesCountInFq_ofs(outputSpeciesCountFileInFqFile.c_str());
	speciesCountInFq_ofs << "Species_Id\tGenus_Id\tFamily_Id\tOrder_Id\tClass_Id\tPhylum_Id\tCount" << endl;
    if(SE_or_PE_bool)
    {
    	string output_fq_SE = outputDir + "merged.SE.fq";
    	ofstream SE_ofs(output_fq_SE.c_str());
	    for(int tmp = 0; tmp < speciesIdVec.size(); tmp++)
	    {
			int tmpSpeciesId = speciesIdVec[tmp];
			int tmpGenusId = genusIdVec[tmp];
			int tmpFamilyId = familyIdVec[tmp];
			int tmpOrderId = orderIdVec[tmp];
			int tmpClassId = classIdVec[tmp];
			int tmpPhylumId = phylumIdVec[tmp];
			int tmpSpeciesCount_ori = speciesCountVec[tmp];
			speciesCountInFq_ofs << speciesIdVec[tmp] << "\t" << genusIdVec[tmp] << "\t"
				<< familyIdVec[tmp] << "\t" << orderIdVec[tmp] << "\t" << classIdVec[tmp]
				<< "\t" << phylumIdVec[tmp] << "\t";
			string tmpTaxoInfo_id = int_to_str(tmpSpeciesId) + "_" + int_to_str(tmpGenusId) 
				+ "_" + int_to_str(tmpFamilyId) + "_" + int_to_str(tmpOrderId) + "_"
				+ int_to_str(tmpClassId) + "_" + int_to_str(tmpPhylumId);
			if(tmpSpeciesCount_ori == 0)
			{
				speciesCountInFq_ofs << "0" << endl;
				continue;
			}
			string tmpInputRawSeFqFilePath = inputFqDir + int_to_str(tmpSpeciesId) + ".fq";
			ifstream tmpRaw_ifs(tmpInputRawSeFqFilePath.c_str());
			int tmpFqReadNum = 0;
			while(!tmpRaw_ifs.eof())
			{
				string tmpStr_1, tmpStr_2, tmpStr_3, tmpStr_4;
				getline(tmpRaw_ifs, tmpStr_1);
				if(tmpStr_1 == "")
					break;
				getline(tmpRaw_ifs, tmpStr_2);
				getline(tmpRaw_ifs, tmpStr_3);
				getline(tmpRaw_ifs, tmpStr_4);
				tmpFqReadNum ++;
				SE_ofs << "@" << tmpTaxoInfo_id << "_" << int_to_str(tmpFqReadNum) << endl
					<< tmpStr_2 << endl << tmpStr_3 << endl << tmpStr_4 << endl;
			}
			speciesCountInFq_ofs << tmpFqReadNum << endl;
			speciesCountVec_inFq[tmp] = tmpFqReadNum;
			tmpRaw_ifs.close();
		}    
    	SE_ofs.close();
    }
    else
    {	
	    string output_fq_PE_1 = outputDir + "merged.PE.1.fq";
	    string output_fq_PE_2 = outputDir + "merged.PE.2.fq";
	    ofstream PE_1_ofs(output_fq_PE_1.c_str());
	    ofstream PE_2_ofs(output_fq_PE_2.c_str());
	    for(int tmp = 0; tmp < speciesIdVec.size(); tmp++)
	    {
			int tmpSpeciesId = speciesIdVec[tmp];
			int tmpGenusId = genusIdVec[tmp];
			int tmpFamilyId = familyIdVec[tmp];
			int tmpOrderId = orderIdVec[tmp];
			int tmpClassId = classIdVec[tmp];
			int tmpPhylumId = phylumIdVec[tmp];
			int tmpSpeciesCount_ori = speciesCountVec[tmp];
			speciesCountInFq_ofs << speciesIdVec[tmp] << "\t" << genusIdVec[tmp] << "\t"
				<< familyIdVec[tmp] << "\t" << orderIdVec[tmp] << "\t" << classIdVec[tmp]
				<< "\t" << phylumIdVec[tmp] << "\t";
			string tmpTaxoInfo_id = int_to_str(tmpSpeciesId) + "_" + int_to_str(tmpGenusId) 
				+ "_" + int_to_str(tmpFamilyId) + "_" + int_to_str(tmpOrderId) + "_"
				+ int_to_str(tmpClassId) + "_" + int_to_str(tmpPhylumId);
			if(tmpSpeciesCount_ori == 0)
			{
				speciesCountInFq_ofs << "0" << endl;
				continue;
			}
			string tmpInputRawPeFqFilePath_1 = inputFqDir + int_to_str(tmpSpeciesId) + ".1.fq";
			string tmpInputRawPeFqFilePath_2 = inputFqDir + int_to_str(tmpSpeciesId) + ".2.fq";
			ifstream tmpRaw_ifs_1(tmpInputRawPeFqFilePath_1.c_str());
			ifstream tmpRaw_ifs_2(tmpInputRawPeFqFilePath_2.c_str());
			int tmpFqReadNum = 0;
			while((!tmpRaw_ifs_1.eof())&&(!tmpRaw_ifs_2.eof()))
			{
				string tmpStr_1, tmpStr_2, tmpStr_3, tmpStr_4;
				string tmpStr_5, tmpStr_6, tmpStr_7, tmpStr_8;
				getline(tmpRaw_ifs_1, tmpStr_1);
				getline(tmpRaw_ifs_2, tmpStr_5);
				if((tmpStr_1 == "")||(tmpStr_5 == ""))
					break;
				getline(tmpRaw_ifs_1, tmpStr_2);
				getline(tmpRaw_ifs_1, tmpStr_3);
				getline(tmpRaw_ifs_1, tmpStr_4);
				getline(tmpRaw_ifs_2, tmpStr_6);
				getline(tmpRaw_ifs_2, tmpStr_7);
				getline(tmpRaw_ifs_2, tmpStr_8);				
				tmpFqReadNum ++;
				PE_1_ofs << "@" << tmpTaxoInfo_id << "_" << int_to_str(tmpFqReadNum) << "/1" 
					<< endl << tmpStr_2 << endl << tmpStr_3 << endl << tmpStr_4 << endl;
				PE_2_ofs << "@" << tmpTaxoInfo_id << "_" << int_to_str(tmpFqReadNum) << "/2" 
					<< endl << tmpStr_6 << endl << tmpStr_7 << endl << tmpStr_8 << endl; 
			}
			speciesCountInFq_ofs << tmpFqReadNum << endl;
			speciesCountVec_inFq[tmp] = tmpFqReadNum;
			tmpRaw_ifs_1.close();
			tmpRaw_ifs_2.close();
		}
		PE_1_ofs.close();
		PE_2_ofs.close();
	}
	speciesCountInFq_ofs.close();
    log_ofs.close();
	return 0;
}