// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#ifndef CHROMOSOMESEQ_INFO_VEC_H
#define CHROMOSOMESEQ_IFNO_VEC_H
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
#include "../../mps3Lib/read_block_test.h"
#include "../../mps3Lib/otherFunc.h"
#include "chromosomeSeq_info.h"

using namespace std;

class ChromosomeSeq_Info_Vec
{
private:
	vector<ChromosomeSeq_Info> chromosomeSeqInfoVec;
	vector<string> invalidChromosomeSeqPathVec;

	string dir_fa;
	vector< pair<int,string> > dirVec_jf; // <KmerLength, path>
	vector< pair<int,string> > dirVec_Kmer; // <KmerLength, path>
	vector< pair<int,string> > dirVec_sortedKmer; // <KmerLength, path>
public:
	ChromosomeSeq_Info_Vec()
	{}

	int return_chromosomeSeqInfoVecSize()
	{
		return chromosomeSeqInfoVec.size();
	}

	int return_chromosomeSeqInfo_accession_id(int index)
	{	
		return chromosomeSeqInfoVec[index].return_accession_id();
	}

	int return_chromosomeSeqInfo_species_id(int index)
	{	
		return chromosomeSeqInfoVec[index].return_species_id();
	}

	int return_chromosomeSeqInfo_genus_id(int index)
	{	
		return chromosomeSeqInfoVec[index].return_genus_id();
	}	

	int return_chromosomeSeqInfo_phylum_id(int index)
	{	
		return chromosomeSeqInfoVec[index].return_phylum_id();
	}

	string return_chromosomeSeqInfo_species_name(int index)
	{	
		return chromosomeSeqInfoVec[index].return_species_name();
	}

	string return_chromosomeSeqInfo_genus_name(int index)
	{	
		return chromosomeSeqInfoVec[index].return_genus_name();
	}	

	string return_chromosomeSeqInfo_phylum_name(int index)
	{	
		return chromosomeSeqInfoVec[index].return_phylum_name();
	}

	void initiate_fa_fileVec_creatDir(string& outputDir_fa)
	{
		outputDir_fa += "/";
		dir_fa = outputDir_fa;
		string cmd_mkdir = "mkdir " + outputDir_fa;
		system(cmd_mkdir.c_str());
		for(int tmp = 0; tmp < chromosomeSeqInfoVec.size(); tmp++)
		{
			string tmpMergedFaFile = outputDir_fa + int_to_str(chromosomeSeqInfoVec[tmp].return_accession_id()) + ".fa" ;
			chromosomeSeqInfoVec[tmp].assign_mergedGenome_fa_path(tmpMergedFaFile);
		}
	}

	void initiate_fa_fileVec_creatDir_mergedSpeciesFa(string& outputDir_fa)
	{
		this->initiate_fa_fileVec_creatDir(outputDir_fa);
		// outputDir_fa += "/";
		// dir_fa = outputDir_fa;
		// string cmd_mkdir = "mkdir " + outputDir_fa;
		// system(cmd_mkdir.c_str());
		// for(int tmp = 0; tmp < chromosomeSeqInfoVec.size(); tmp++)
		// {
		// 	string tmpMergedFaFile = outputDir_fa + int_to_str(chromosomeSeqInfoVec[tmp].return_accession_id()) + ".fa" ;
		// 	chromosomeSeqInfoVec[tmp].assign_mergedSpecies_fa_path(tmpMergedFaFile);
		// }
	}

	void initiate_fa_fileVec_keepDir(string& outputDir_fa)
	{
		outputDir_fa += "/";
		dir_fa = outputDir_fa;
		//string cmd_mkdir = "mkdir " + outputDir_fa;
		//system(cmd_mkdir.c_str());
		for(int tmp = 0; tmp < chromosomeSeqInfoVec.size(); tmp++)
		{
			string tmpMergedFaFile = outputDir_fa + int_to_str(chromosomeSeqInfoVec[tmp].return_accession_id()) + ".fa" ;
			chromosomeSeqInfoVec[tmp].assign_mergedGenome_fa_path(tmpMergedFaFile);
		}
	}

	void generateMergedGenome_fa()
	{
		//string outputDir_fa = dir_fa;
		for(int tmp = 0; tmp < chromosomeSeqInfoVec.size(); tmp++)
		{
			string tmpMergeFaFile_cmd = "cat";
			int tmpSourceFaFileNum = chromosomeSeqInfoVec[tmp].return_chromosomeSeqFileNum();
			for(int tmp2 = 0; tmp2 < tmpSourceFaFileNum; tmp2++)
			{
				string tmpSourceFaFilePath = chromosomeSeqInfoVec[tmp].return_chromosomeSeqFilePath(tmp2);
				tmpMergeFaFile_cmd += " ";
				tmpMergeFaFile_cmd += tmpSourceFaFilePath;
			}
			tmpMergeFaFile_cmd += " > ";
			tmpMergeFaFile_cmd += chromosomeSeqInfoVec[tmp].return_mergedGenome_fa_path();
			system(tmpMergeFaFile_cmd.c_str());
		}
	}

	void generateMergedSpecies_fa()
	{
		this->generateMergedGenome_fa();
	}

	void initiate_jf_fileVec_creatDir(string& outputDir_jf, int Kmer_length)
	{
		outputDir_jf += "/";
		string dir_jf = outputDir_jf;
		dirVec_jf.push_back(pair<int,string>(Kmer_length, dir_jf));
		string cmd_mkdir = "mkdir " + outputDir_jf;
		system(cmd_mkdir.c_str());
		for(int tmp = 0; tmp < chromosomeSeqInfoVec.size(); tmp++)
		{		
			string tmpJfFile = outputDir_jf + int_to_str(chromosomeSeqInfoVec[tmp].return_accession_id()) 
				+ "_" + int_to_str(Kmer_length) + "mer.jf";
			chromosomeSeqInfoVec[tmp].assign_mergedGenome_jf_path(tmpJfFile, Kmer_length);
		}
	}

	void generateMergedGenome_jf(string& jellyFishBin, int Kmer_length, int threads_num, string& bitMapSizeStr)
	{
		for(int tmp = 0; tmp < chromosomeSeqInfoVec.size(); tmp++)
		{
			string tmpSourceFaFilePath = chromosomeSeqInfoVec[tmp].return_mergedGenome_fa_path();
			string tmpJfFilePath = chromosomeSeqInfoVec[tmp].return_mergedGenome_jf_path_with_KmerLength(Kmer_length);
			string tmpGetJfFile_cmd = jellyFishBin + "/jellyfish count -o " + tmpJfFilePath + " -m " + int_to_str(Kmer_length)
				+ " -t " + int_to_str(threads_num) + " -s " + bitMapSizeStr + " -C " + tmpSourceFaFilePath;
			cout << "tmpGetJfFile_cmd: " << endl << tmpGetJfFile_cmd << endl;
			system(tmpGetJfFile_cmd.c_str());
		}
	}

	void initiate_Kmer_fileVec_creatDir(string& outputDir_Kmer, int Kmer_length)
	{
		outputDir_Kmer += "/";
		string dir_Kmer = outputDir_Kmer;
		dirVec_Kmer.push_back(pair<int,string>(Kmer_length, dir_Kmer));
		string cmd_mkdir = "mkdir " + outputDir_Kmer;
		system(cmd_mkdir.c_str());
		for(int tmp = 0; tmp < chromosomeSeqInfoVec.size(); tmp++)
		{		
			string tmpKmerFile = outputDir_Kmer + int_to_str(chromosomeSeqInfoVec[tmp].return_accession_id()) 
				+ "_" + int_to_str(Kmer_length) + "mer.Kmer";
			chromosomeSeqInfoVec[tmp].assign_mergedGenome_Kmer_path(tmpKmerFile, Kmer_length);
		}		
	}

	void generateMergedGenome_Kmer(string& jellyFishBin, int Kmer_length)
	{
		for(int tmp = 0; tmp < chromosomeSeqInfoVec.size(); tmp++)
		{
			string tmpSourceJfFilePath = chromosomeSeqInfoVec[tmp].return_mergedGenome_jf_path_with_KmerLength(Kmer_length);
			string tmpKmerFilePath = chromosomeSeqInfoVec[tmp].return_mergedGenome_Kmer_path_with_KmerLength(Kmer_length);
			string tmpDump_cmd = jellyFishBin + "/jellyfish dump -t -c -o " + tmpKmerFilePath
				+ " " + tmpSourceJfFilePath;
			cout << "tmpDump_cmd: " << tmpDump_cmd << endl;
			system(tmpDump_cmd.c_str());
		}		
	}

	void initiate_sortedKmer_fileVec_creatDir(string& outputDir_sortedKmer, int Kmer_length)
	{
		outputDir_sortedKmer += "/";
		string dir_sortedKmer = outputDir_sortedKmer;
		dirVec_sortedKmer.push_back(pair<int,string>(Kmer_length, dir_sortedKmer));
		string cmd_mkdir = "mkdir " + outputDir_sortedKmer;
		system(cmd_mkdir.c_str());
		for(int tmp = 0; tmp < chromosomeSeqInfoVec.size(); tmp++)
		{		
			string tmpSortedKmerFile = outputDir_sortedKmer + int_to_str(chromosomeSeqInfoVec[tmp].return_accession_id()) 
				+ "_" + int_to_str(Kmer_length) + "mer.sortedKmer";
			chromosomeSeqInfoVec[tmp].assign_mergedGenome_sortedKmer_path(tmpSortedKmerFile, Kmer_length);
		}			
	}

	void generateMergedGenome_sortedKmer(string& tmpSortDir, int Kmer_length)
	{
		string cmd_mkdir_tmpSortDir = "mkdir " + tmpSortDir;
		system(cmd_mkdir_tmpSortDir.c_str());
		for(int tmp = 0; tmp < chromosomeSeqInfoVec.size(); tmp++)
		{
			string tmpSourceKmerFilePath = chromosomeSeqInfoVec[tmp].return_mergedGenome_Kmer_path_with_KmerLength(Kmer_length);	
			string tmpSortedKmerFilePath = chromosomeSeqInfoVec[tmp].return_mergedGenome_sortedKmer_path_with_KmerLength(Kmer_length);
			string tmpSortKmer_cmd = "sort -k1 -T " + tmpSortDir + " " + tmpSourceKmerFilePath + " > " + tmpSortedKmerFilePath;
			cout << "tmpSortKmer_cmd: " << tmpSortKmer_cmd << endl;
			system(tmpSortKmer_cmd.c_str());
		}
	}

	bool accessionIdAlreadyExists_bool(int tmpAccessionId, int& existingIndexInVec)
	{
		for(int tmp = 0; tmp < chromosomeSeqInfoVec.size(); tmp++)
		{
			int tmpSpeciesAccessionIdInVec = chromosomeSeqInfoVec[tmp].return_accession_id();
			if(tmpAccessionId == tmpSpeciesAccessionIdInVec)
			{
				existingIndexInVec = tmp;
				return true;
			}
		}
		return false;
	}

	void initiate(string& bacteriaChromosomeSeqPath2TaxIdFile, Taxonomy_Info& taxonomyInfo)
	{
		ifstream chromosomeSeqPath2TaxId_ifs(bacteriaChromosomeSeqPath2TaxIdFile.c_str());
		while(!chromosomeSeqPath2TaxId_ifs.eof())
		{
			string tmpStr;
			getline(chromosomeSeqPath2TaxId_ifs, tmpStr);
			if(tmpStr == "")
				break;
			int tabLoc_1 = tmpStr.find("\t");
			int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
			string tmpChromosomeSeqPath = tmpStr.substr(0, tabLoc_1);
			string tmpAccessionIdStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
			string tmpTaxonomyInfoStr = tmpStr.substr(tabLoc_2 + 1);
			if((tmpAccessionIdStr == "-1")||(tmpAccessionIdStr == "UNKNOWN"))
				invalidChromosomeSeqPathVec.push_back(tmpChromosomeSeqPath);
			int tmpAccessionId = atoi(tmpAccessionIdStr.c_str());
			int tmpExistingIndexInVec;
			bool tmpAccessionId_alreadyExistsBool = this->accessionIdAlreadyExists_bool(
				tmpAccessionId, tmpExistingIndexInVec);
			if(tmpAccessionId_alreadyExistsBool)
				chromosomeSeqInfoVec[tmpExistingIndexInVec].update_newChromosomeSeq(tmpChromosomeSeqPath);
			else
			{
				ChromosomeSeq_Info tmpChromosomeSeqInfo;
				tmpChromosomeSeqInfo.initiate_1stChromosomeSeq(tmpChromosomeSeqPath, tmpAccessionId, 
					tmpTaxonomyInfoStr, taxonomyInfo);
				chromosomeSeqInfoVec.push_back(tmpChromosomeSeqInfo);
			}
		}
		chromosomeSeqPath2TaxId_ifs.close();
	}
	
	void initiate(string& bacteriaChromosomeSeqPath2TaxIdFile)
	{
		ifstream chromosomeSeqPath2TaxId_ifs(bacteriaChromosomeSeqPath2TaxIdFile.c_str());
		while(!chromosomeSeqPath2TaxId_ifs.eof())
		{
			string tmpStr;
			getline(chromosomeSeqPath2TaxId_ifs, tmpStr);
			if(tmpStr == "")
				break;
			int tabLoc_1 = tmpStr.find("\t");
			int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
			string tmpChromosomeSeqPath = tmpStr.substr(0, tabLoc_1);
			string tmpAccessionIdStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
			string tmpTaxonomyInfoStr = tmpStr.substr(tabLoc_2 + 1);
			if((tmpAccessionIdStr == "-1")||(tmpAccessionIdStr == "UNKNOWN"))
			{
				invalidChromosomeSeqPathVec.push_back(tmpChromosomeSeqPath);
				continue;
			}
			int tmpAccessionId = atoi(tmpAccessionIdStr.c_str());
			int tmpExistingIndexInVec;
			bool tmpAccessionId_alreadyExistsBool = this->accessionIdAlreadyExists_bool(
				tmpAccessionId, tmpExistingIndexInVec);
			if(tmpAccessionId_alreadyExistsBool)
				chromosomeSeqInfoVec[tmpExistingIndexInVec].update_newChromosomeSeq(tmpChromosomeSeqPath);
			else
			{
				ChromosomeSeq_Info tmpChromosomeSeqInfo;
				tmpChromosomeSeqInfo.initiate_1stChromosomeSeq(tmpChromosomeSeqPath, tmpAccessionId, tmpTaxonomyInfoStr);
				chromosomeSeqInfoVec.push_back(tmpChromosomeSeqInfo);
			}
		}
		chromosomeSeqPath2TaxId_ifs.close();
	}

	void initiate_mergedSpeciesFa(string& bacteriaChromosomeSeqPath2TaxIdFile)
	{
		ifstream chromosomeSeqPath2TaxId_ifs(bacteriaChromosomeSeqPath2TaxIdFile.c_str());
		while(!chromosomeSeqPath2TaxId_ifs.eof())
		{
			string tmpStr;
			getline(chromosomeSeqPath2TaxId_ifs, tmpStr);
			if(tmpStr == "")
				break;
			int tabLoc_1 = tmpStr.find("\t");
			int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
			int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);

			string tmpChromosomeSeqPath = tmpStr.substr(0, tabLoc_1);
			string tmpAccessionIdStr = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
			string tmpTaxonomyInfoStr = tmpStr.substr(tabLoc_2 + 1);
			if((tmpAccessionIdStr == "-1")||(tmpAccessionIdStr == "UNKNOWN"))
			{
				invalidChromosomeSeqPathVec.push_back(tmpChromosomeSeqPath);
				continue;
			}
			int tmpAccessionId = atoi(tmpAccessionIdStr.c_str());
			int tmpExistingIndexInVec;
			bool tmpAccessionId_alreadyExistsBool = this->accessionIdAlreadyExists_bool(
				tmpAccessionId, tmpExistingIndexInVec);
			if(tmpAccessionId_alreadyExistsBool)
				chromosomeSeqInfoVec[tmpExistingIndexInVec].update_newChromosomeSeq(tmpChromosomeSeqPath);
			else
			{
				ChromosomeSeq_Info tmpChromosomeSeqInfo;
				tmpChromosomeSeqInfo.initiate_1stChromosomeSeq(tmpChromosomeSeqPath, tmpAccessionId, tmpTaxonomyInfoStr);
				chromosomeSeqInfoVec.push_back(tmpChromosomeSeqInfo);
			}
		}
		chromosomeSeqPath2TaxId_ifs.close();
	}	

	void printMergedGenome_fa_info(string& mergedGenome_fa_info_file_path)
	{
		ofstream mergedGenome_fa_info_ofs(mergedGenome_fa_info_file_path.c_str());
		for(int tmp = 0; tmp < chromosomeSeqInfoVec.size(); tmp++)
		{
			string tmpMergedGenome_fa_path = chromosomeSeqInfoVec[tmp].return_mergedGenome_fa_path();
			int tmp_accession_id = chromosomeSeqInfoVec[tmp].return_accession_id();
			int tmp_species_id = chromosomeSeqInfoVec[tmp].return_species_id();
			int tmp_genus_id = chromosomeSeqInfoVec[tmp].return_genus_id();
			int tmp_family_id = chromosomeSeqInfoVec[tmp].return_family_id();
			int tmp_order_id = chromosomeSeqInfoVec[tmp].return_order_id();
			int tmp_class_id = chromosomeSeqInfoVec[tmp].return_class_id();
			int tmp_phylum_id = chromosomeSeqInfoVec[tmp].return_phylum_id();
			mergedGenome_fa_info_ofs << tmpMergedGenome_fa_path << "\t" << tmp_accession_id << "\t"
				<< tmp_species_id << "\t" << tmp_genus_id << "\t" << tmp_family_id << "\t"
				<< tmp_order_id << "\t" << tmp_class_id << "\t" << tmp_phylum_id << endl;
		}
		mergedGenome_fa_info_ofs.close();
	}

	void printMergedSpecies_fa_info(string& mergedSpecies_fa_info_file_path)
	{
		this->printMergedGenome_fa_info(mergedSpecies_fa_info_file_path);
	}

	void printMergedGenome_jf_info(string& mergedGenome_jf_info_file_path, int Kmer_length)
	{
		ofstream mergedGenome_jf_info_ofs(mergedGenome_jf_info_file_path.c_str());
		for(int tmp = 0; tmp < chromosomeSeqInfoVec.size(); tmp++)
		{
			string tmpMergedGenome_jf_path = chromosomeSeqInfoVec[tmp].return_mergedGenome_jf_path_with_KmerLength(Kmer_length);
			int tmp_accession_id = chromosomeSeqInfoVec[tmp].return_accession_id();
			int tmp_species_id = chromosomeSeqInfoVec[tmp].return_species_id();
			int tmp_genus_id = chromosomeSeqInfoVec[tmp].return_genus_id();
			int tmp_family_id = chromosomeSeqInfoVec[tmp].return_family_id();
			int tmp_order_id = chromosomeSeqInfoVec[tmp].return_order_id();
			int tmp_class_id = chromosomeSeqInfoVec[tmp].return_class_id();
			int tmp_phylum_id = chromosomeSeqInfoVec[tmp].return_phylum_id();
			mergedGenome_jf_info_ofs << tmpMergedGenome_jf_path << "\t" << Kmer_length << "\t" 
				<< tmp_accession_id << "\t" << tmp_species_id << "\t" << tmp_genus_id << "\t" << tmp_family_id 
				<< "\t" << tmp_order_id << "\t" << tmp_class_id << "\t" << tmp_phylum_id << endl;
		}
		mergedGenome_jf_info_ofs.close();
	}

	void printMergedGenome_Kmer_info(string& mergedGenome_Kmer_info_file_path, int Kmer_length)
	{
		ofstream mergedGenome_Kmer_info_ofs(mergedGenome_Kmer_info_file_path.c_str());
		for(int tmp = 0; tmp < chromosomeSeqInfoVec.size(); tmp++)
		{
			string tmpMergedGenome_Kmer_path = chromosomeSeqInfoVec[tmp].return_mergedGenome_Kmer_path_with_KmerLength(Kmer_length);
			int tmp_accession_id = chromosomeSeqInfoVec[tmp].return_accession_id();
			int tmp_species_id = chromosomeSeqInfoVec[tmp].return_species_id();
			int tmp_genus_id = chromosomeSeqInfoVec[tmp].return_genus_id();
			int tmp_family_id = chromosomeSeqInfoVec[tmp].return_family_id();
			int tmp_order_id = chromosomeSeqInfoVec[tmp].return_order_id();
			int tmp_class_id = chromosomeSeqInfoVec[tmp].return_class_id();
			int tmp_phylum_id = chromosomeSeqInfoVec[tmp].return_phylum_id();
			mergedGenome_Kmer_info_ofs << tmpMergedGenome_Kmer_path << "\t" << Kmer_length << "\t" 
				<< tmp_accession_id << "\t" << tmp_species_id << "\t" << tmp_genus_id << "\t" << tmp_family_id 
				<< "\t" << tmp_order_id << "\t" << tmp_class_id << "\t" << tmp_phylum_id << endl;
		}
		mergedGenome_Kmer_info_ofs.close();
	}

	void printMergedGenome_sortedKmer_info(string& mergedGenome_sortedKmer_info_file_path, int Kmer_length)
	{
		ofstream mergedGenome_sortedKmer_info_ofs(mergedGenome_sortedKmer_info_file_path.c_str());
		for(int tmp = 0; tmp < chromosomeSeqInfoVec.size(); tmp++)
		{
			string tmpMergedGenome_sortedKmer_path = chromosomeSeqInfoVec[tmp].return_mergedGenome_sortedKmer_path_with_KmerLength(Kmer_length);
			int tmp_accession_id = chromosomeSeqInfoVec[tmp].return_accession_id();
			int tmp_species_id = chromosomeSeqInfoVec[tmp].return_species_id();
			int tmp_genus_id = chromosomeSeqInfoVec[tmp].return_genus_id();
			int tmp_family_id = chromosomeSeqInfoVec[tmp].return_family_id();
			int tmp_order_id = chromosomeSeqInfoVec[tmp].return_order_id();
			int tmp_class_id = chromosomeSeqInfoVec[tmp].return_class_id();
			int tmp_phylum_id = chromosomeSeqInfoVec[tmp].return_phylum_id();
			mergedGenome_sortedKmer_info_ofs << tmpMergedGenome_sortedKmer_path << "\t" << Kmer_length << "\t" 
				<< tmp_accession_id << "\t" << tmp_species_id << "\t" << tmp_genus_id << "\t" << tmp_family_id 
				<< "\t" << tmp_order_id << "\t" << tmp_class_id << "\t" << tmp_phylum_id << endl;
		}
		mergedGenome_sortedKmer_info_ofs.close();
	}
};
#endif