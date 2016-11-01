#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>
#include "../mps3Lib/read_block_test.h"
#include "../mps3Lib/otherFunc.h"
#include "../mps3Lib/index_info.h"
#include "general/chromosomeSeq_info_vec.h"
#include "general/species_info.h"
#include "general/genus_info.h"
#include "general/phylum_info.h"
time_t nowtime;
struct tm *local;

using namespace std;

void getSpecificTaxoLevelJfFile2sourceChromosomeSeqJfVecPairVec(
	vector< pair<string, vector<string> > >& specificTaxoLevelJfFile2sourceChromosomeSeqJfVecPairVec, string& outputFolder_jf, 
	vector< pair<int, vector<int> > >& taxoId2sourceChromosomeSeqIdVecPairVec, string& inputJfFolder)
{
	outputFolder_jf += "/";
	inputJfFolder += "/"
	for(int tmp = 0; tmp < taxoId2sourceChromosomeSeqIdVecPairVec.size(); tmp++)
	{
		int tmpOutputJf_id = taxoId2sourceChromosomeSeqIdVecPairVec[tmp].first;
		string tmpOutputJfFile = outputFolder_jf + int_to_str(tmpOutputJf_id) + ".jf";
		vector<string> tmpInputJfFileVec;
		for(int tmp2 = 0; tmp2 < (taxoId2sourceChromosomeSeqIdVecPairVec[tmp].second).size(); tmp2++)
		{
			int tmpInputJf_id = (taxoId2sourceChromosomeSeqIdVecPairVec[tmp].second)[tmp2];
			string tmpInputSourceJf = inputJfFolder + int_to_str(tmpInputJf_id) + ".jf";
			tmpInputJfFileVec.push_back(tmpInputSourceJf);
		}
		specificTaxoLevelJfFile2sourceChromosomeSeqJfVecPairVec.push_back(pair<string, vector<string> >(
			tmpOutputJfFile, tmpInputJfFileVec));
	}

}

int main(int argc, char** argv)
{
	if(argc != 5)
	{	
		cout << "Executable JellyFishBin inputChromosomeSeq2taxoInfoFile inputJfFolder outputFolder taxo_level(species/genus/phylum)" << endl;
		exit(1);
	}
	string outputFolder = argv[4];
	outputFolder += "/";
	string cmd_mkdir = "mkdir " + outputFolder;
	system(cmd_mkdir.c_str());

	string log_file = outputFolder + "log";
	ofstream log_ofs(log_file.c_str());
	log_ofs << "Command Line:" << endl;
	for(int tmp = 0; tmp < argc; tmp++)
		log_ofs << argv[1] << endl;
	log_ofs << endl;

	string JellyFishBin = argv[1];
	JellyFishBin += "/";
	string inputChromosomeSeq2taxoInfoFile = argv[2];
	string inputJfFolder = argv[3];
	inputJfFolder += "/";

	string taxo_level_str = argv[5];
	int taxo_level_int = -1; // 1-phylum, 2-class, 3-order, 4-family, 5-genus, 6-species
	if((taxo_level_str == "phylum")||(taxo_level_str == "Phylum")||(taxo_level_str == "PHYLUM"))
		taxo_level_int = 1;
	else if((taxo_level_str == "genus")||(taxo_level_str == "Genus")||(taxo_level_str == "GENUS"))
		taxo_level_int = 5;
	else if((taxo_level_str == "species")||(taxo_level_str == "Species")||(taxo_level_str == "SPECIES"))
		taxo_level_int = 6;
	else
	{
		cout << "invalid taxo_level_str: " << taxo_level_str << endl;
		cout << "Please set as: Phylum or Genus or Species" << endl;
		exit(1); 
	}


	cout << "start to initiate ChromosomeSeq_Info_Vec" << endl;
	ChromosomeSeq_Info_Vec tmpChrSeqInfoVec;
	tmpChrSeqInfoVec.initiate(inputChromosomeSeq2taxoInfoFile);	

	cout << "start to initiate Spcies_Info" << endl;
	Species_Info tmpSpeciesInfo;
	tmpSpeciesInfo.initiate(tmpChrSeqInfoVec);

	cout << "start to initiate Genus_Info" << endl;
	Genus_Info tmpGenusInfo;
	tmpGenusInfo.initiate(tmpChrSeqInfoVec);

	cout << "start to initiate Species_Info" << endl;
	Phylum_Info tmpPhylumInfo;
	tmpPhylumInfo.initiate(tmpChrSeqInfoVec);

	log_ofs << endl << "start to generateTaxoId2sourceChromosomeSeqIdVecPairVec" << endl << endl;
	vector< pair<int, vector<int> > > taxoId2sourceChromosomeSeqIdVecPairVec;
	if(taxo_level_int == 1) // phylum
		tmpPhylumInfo.generateTaxoId2sourceChromosomeSeqIdVecPairVec(taxoId2sourceChromosomeSeqIdVecPairVec);
	else if(taxo_level_int == 5) // genus
		tmpGenusInfo.generateTaxoId2sourceChromosomeSeqIdVecPairVec(taxoId2sourceChromosomeSeqIdVecPairVec);
	else // species, taxo_level_int == 6
		tmpSpeciesInfo.generateTaxoId2sourceChromosomeSeqIdVecPairVec(taxoId2sourceChromosomeSeqIdVecPairVec);

	string outputFolder_jf = outputFolder + "_jf/";
	string outputFolder_Kmer = outputFolder + "_Kmer/";
	string outputFolder_sortedKmer = outputFolder + "_sortedKmer/";
	string cmd_mkdir_jf = "mdkir " + outputFolder_jf;
	string cmd_mkdir_Kmer = "mdkir " + outputFolder_Kmer;
	string cmd_mkdir_sortedKmer = "mdkir " + outputFolder_sortedKmer;
	system(cmd_mkdir_jf.c_str());
	system(cmd_mkdir_Kmer.c_str());
	system(cmd_mkdir_sortedKmer.c_str());	

	log_ofs << endl << "start to getSpecificTaxoLevelJfFile2sourceChromosomeSeqJfVecPairVec" << endl << endl;
	vector< pair<string, vector<string> > > specificTaxoLevelJfFile2sourceChromosomeSeqJfVecPairVec;
	getSpecificTaxoLevelJfFile2sourceChromosomeSeqJfVecPairVec(specificTaxoLevelJfFile2sourceChromosomeSeqJfVecPairVec,
		outputFolder_jf, taxoId2sourceChromosomeSeqIdVecPairVec, inputJfFolder);
	
	// get Jf files
	log_ofs << endl << "start to get jf files" << endl << endl;
	for(int tmp = 0; tmp < specificTaxoLevelJfFile2sourceChromosomeSeqJfVecPairVec.size(); tmp++)
	{
		string tmpOutputJfFile = specificTaxoLevelJfFile2sourceChromosomeSeqJfVecPairVec[tmp].first;
		int tmpInputJfNum = (specificTaxoLevelJfFile2sourceChromosomeSeqJfVecPairVec[tmp].second).size();
		if(tmpInputJfNum == 1)
		{
			string tmpSingleInputJfFile = (specificTaxoLevelJfFile2sourceChromosomeSeqJfVecPairVec[tmp].second)[0];
			string cmd_cp_singleJf = "cp " + tmpSingleInputJfFile + " " + tmpOutputJfFile;
			log_ofs << "tmpCmd: " << cmd_cp_singleJf << endl;
			system(cmd_cp_singleJf.c_str());
		}
		else if(tmpInputJfNum == 2)
		{
			string tmpInputJfFile_1 = (specificTaxoLevelJfFile2sourceChromosomeSeqJfVecPairVec[tmp].second)[0];
			string tmpInputJfFile_2 = (specificTaxoLevelJfFile2sourceChromosomeSeqJfVecPairVec[tmp].second)[1];
			string cmd_mergeJf = JellyFishBin + "jellyfish merge -o " + tmpOutputJfFile + " " 
				+ tmpInputJfFile_1 + " " + tmpInputJfFile_2;
			log_ofs << "tmpCmd: " << cmd_mergeJf << endl;
			system(cmd_mergeJf.c_str());
		}
		else if(tmpInputJfNum == 0)
		{
			cout << "tmpInputJfNum == 0, error !" << endl;
			exit(1);
		}
		else // tmpInputJfNum > 2
		{
			string cmd_mergeJf = JellyFishBin + "jellyfish merge -o " + tmpOutputJfFile;
			for(int tmp2 = 0; tmp2 < tmpInputJfNum; tmp2++)
			{
				cmd_mergeJf += " ";
				cmd_mergeJf += (specificTaxoLevelJfFile2sourceChromosomeSeqJfVecPairVec[tmp].second)[tmp2];
			}
			log_ofs << "tmpCmd: " << cmd_mergeJf << endl;
			system(cmd_mergeJf.c_str());
		}
	}

	// get Kmer files
	log_ofs << endl << "start to get Kmer files" << endl << endl;
	for(int tmp = 0; tmp < specificTaxoLevelJfFile2sourceChromosomeSeqJfVecPairVec.size(); tmp++)
	{
		string tmpInputJfFile = specificTaxoLevelJfFile2sourceChromosomeSeqJfVecPairVec[tmp].first;
		string tmpOutputKmerFile = outputFolder_Kmer + int_to_str(taxoId2sourceChromosomeSeqIdVecPairVec[tmp].first) + ".Kmer";
		string cmd_jf2Kmer = JellyFishBin + "jellyfish dump -t -c -o " + tmpOutputKmerFile + " " + tmpInputJfFile;
		log_ofs << "tmpCmd: " << cmd_jf2Kmer << endl;
		system(cmd_jf2Kmer.c_str());
	}

	// get sortedKmer files
	log_ofs << endl << "start to get sortedKmer files" << endl << endl;
	string outputFolder_tmpSortDir = outputFolder + "tmpSortDir/";
	string cmd_mkdir_tmpSortDir = "mkdir " + outputFolder_tmpSortDir;
	system(cmd_mkdir_tmpSortDir.c_str());
	for(int tmp = 0; tmp < taxoId2sourceChromosomeSeqIdVecPairVec.size(); tmp++)
	{
		string tmpInputKmerFile = outputFolder_Kmer + int_to_str(taxoId2sourceChromosomeSeqIdVecPairVec[tmp].first) + ".Kmer";
		string tmpOutputSortedKmerFile = outputFolder_sortedKmer + int_to_str(taxoId2sourceChromosomeSeqIdVecPairVec[tmp].first) + ".sortedKmer";
		string cmd_sortKmer = "sort -k1 -T " + outputFolder_tmpSortDir + " " + tmpInputKmerFile + " > " + tmpOutputSortedKmerFile;
		log_ofs << "tmpCmd: " << cmd_sortKmer << endl;
		system(cmd_sortKmer.c_str());
	}
	log_ofs << "All jobs done!" << endl;
	log_ofs.close();
	return 0;
}