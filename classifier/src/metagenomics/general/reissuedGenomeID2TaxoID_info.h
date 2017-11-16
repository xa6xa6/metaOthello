// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#ifndef REISSUEDGENOMEID2TAXOID_INFO_H
#define REISSUEDGENOMEID2TAXOID_IFNO_H
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
//#include "NCBIfullTaxoID2Name_info.h"
using namespace std;

class ReissuedGenomeID2TaxoID_Info
{
private:
	vector<int> oriGenomeIdVec;

	vector<int> speciesIdVec;
	vector<int> genusIdVec;
	vector<int> phylumIdVec;

	vector<string> speciesNameVec;
	vector<string> genusNameVec;
	vector<string> phylumNameVec;

	vector<int> reissuedSpeciesIdVec;
	vector<int> reissuedGenusIdVec;
	vector<int> reissuedPhylumIdVec;

	// reissued id 2 taxoId and name
	vector< pair<int, string> > speciesIdNamePairVec;
	vector< pair<int, string> > genusIdNamePairVec;
	vector< pair<int, string> > phylumIdNamePairVec;

public:
	ReissuedGenomeID2TaxoID_Info()
	{}

	void get_KmerSetSpecificClass(int tmpReissuedId_genome,
		int& tmpKmerSetSpecificClassId_genome, int& tmpKmerSetSpecificClassId_species,
		int& tmpKmerSetSpecificClassId_genus, int& tmpKmerSetSpecificClassId_phylum)
	{
		int tmpReissuedId_species = reissuedSpeciesIdVec[tmpReissuedId_genome];
		int tmpReissuedId_genus = reissuedGenusIdVec[tmpReissuedId_genome];
		int tmpReissuedId_phylum = reissuedPhylumIdVec[tmpReissuedId_genome];

		int total_genome_num = this->return_genome_num();
		int total_species_num = this->return_species_num();
		int total_genus_num = this->return_genus_num();
		int total_phylum_num = this->return_phylum_num();

		tmpKmerSetSpecificClassId_genome = tmpReissuedId_genome + 1;
		tmpKmerSetSpecificClassId_species = total_genome_num + 1 + tmpReissuedId_species;
		tmpKmerSetSpecificClassId_genus = total_genome_num + total_species_num + 1 + tmpReissuedId_genus;
		tmpKmerSetSpecificClassId_phylum = total_genome_num + total_species_num + total_genus_num + 1 + tmpReissuedId_phylum;
	}

	void generate_reissuedTaxoIdVec_taxoIdNamePairVec()
	{
		for(int tmp = 0; tmp < oriGenomeIdVec.size(); tmp++)
		{
			//cout << endl << "tmpGenomeIndex: " << tmp << endl;
			//cout << "tmpGenomeId: " << tmp + 1 << endl;
			int tmpSpeciesId = speciesIdVec[tmp];
			int tmpGenusId = genusIdVec[tmp];
			int tmpPhylumId = phylumIdVec[tmp];
			string tmpSpeciesName = speciesNameVec[tmp];
			string tmpGenusName = genusNameVec[tmp];
			string tmpPhylumName = phylumNameVec[tmp];
			//cout << "tmpPhylumId: " << tmpPhylumId << endl;
			//cout << "tmpPhylumName: " << tmpPhylumName << endl;

			// species
			int currentSpeciesVecSize = speciesIdNamePairVec.size();
			bool tmpSpeciesAlreadyExists_bool = false;
			for(int tmp2 = 0; tmp2 < currentSpeciesVecSize; tmp2++)
			{
				int tmpCurrentSpeciesId = speciesIdNamePairVec[tmp2].first;
				if(tmpSpeciesId == tmpCurrentSpeciesId)
				{	
					reissuedSpeciesIdVec.push_back(tmp2);
					tmpSpeciesAlreadyExists_bool = true;
					break;
				}
			}
			if(!tmpSpeciesAlreadyExists_bool)
			{
				reissuedSpeciesIdVec.push_back(currentSpeciesVecSize);
				speciesIdNamePairVec.push_back(pair<int,string>(tmpSpeciesId, tmpSpeciesName));
			}
			// genus 
			int currentGenusVecSize = genusIdNamePairVec.size();
			bool tmpGenusAlreadyExists_bool = false;
			for(int tmp2 = 0; tmp2 < currentGenusVecSize; tmp2++)
			{
				int tmpCurrentGenusId = genusIdNamePairVec[tmp2].first;
				if(tmpGenusId == tmpCurrentGenusId)
				{
					reissuedGenusIdVec.push_back(tmp2);
					tmpGenusAlreadyExists_bool = true;
					break;
				}
			}
			if(!tmpGenusAlreadyExists_bool)
			{
				reissuedGenusIdVec.push_back(currentGenusVecSize);
				genusIdNamePairVec.push_back(pair<int,string>(tmpGenusId, tmpGenusName));
			}
			// phylum
			int currentPhylumVecSize = phylumIdNamePairVec.size();
			//cout << "currentPhylumVecSize: " << currentPhylumVecSize << endl;
			bool tmpPhylumAlreadyExists_bool = false;
			for(int tmp2 = 0; tmp2 < currentPhylumVecSize; tmp2++)
			{
				int tmpCurrentPhylumId = phylumIdNamePairVec[tmp2].first;
				if(tmpPhylumId == tmpCurrentPhylumId)
				{
					//cout << "reissued phylumId: " << tmp2 << endl;
					reissuedPhylumIdVec.push_back(tmp2);
					tmpPhylumAlreadyExists_bool = true;
					break;
				}
			}
			if(!tmpPhylumAlreadyExists_bool)
			{
				reissuedPhylumIdVec.push_back(currentPhylumVecSize);
				//cout << "reissued phylumId: " << currentPhylumVecSize << endl;
				phylumIdNamePairVec.push_back(pair<int,string>(tmpPhylumId, tmpPhylumName));
			}
		}
	}

	void return_phylum_Id_name_from_reissuedId(int tmpReissuedId, int& tmpPhylumId, string& tmpPhylumName)
	{
		if(tmpReissuedId < 0)
		{
			tmpPhylumId = -999;
			tmpPhylumName = "NA";
			return;
		}
		tmpPhylumId = phylumIdNamePairVec[tmpReissuedId].first;
		tmpPhylumName = phylumIdNamePairVec[tmpReissuedId].second;
		return;
	}

	void return_genus_Id_name_from_reissuedId(int tmpReissuedId, int& tmpGenusId, string& tmpGenusName)
	{
		if(tmpReissuedId < 0)
		{
			tmpGenusId = -999;
			tmpGenusName = "NA";
			return;
		}
		tmpGenusId = genusIdNamePairVec[tmpReissuedId].first;
		tmpGenusName = genusIdNamePairVec[tmpReissuedId].second;
		return;
	}

	void return_species_Id_name_from_reissuedId(int tmpReissuedId, int& tmpSpeciesId, string& tmpSpeciesName)
	{
		if(tmpReissuedId < 0)
		{
			tmpSpeciesId = -999;
			tmpSpeciesName = "NA";
			return;
		}
		tmpSpeciesId = speciesIdNamePairVec[tmpReissuedId].first;
		tmpSpeciesName = speciesIdNamePairVec[tmpReissuedId].second;
		return;
	}	

	int return_reissuedTaxoId_from_genomeId(int tmpGenomeId, int rank)
	{
		if((tmpGenomeId <= 0)||(tmpGenomeId > oriGenomeIdVec.size()))
			return -1;

		if(rank == 3)
			return this->return_reissuedPhylumId_from_indexInOriGenomeIdVec(tmpGenomeId - 1);
		else if(rank == 7)
			return this->return_reissuedGenusId_from_indexInOriGenomeIdVec(tmpGenomeId - 1);
		else if(rank == 8)
			return this->return_reissuedSpeciesId_from_indexInOriGenomeIdVec(tmpGenomeId - 1);
		else
		{
			cout << "invalid rank: " << rank << endl;
			exit(1);
		}
	}

	int return_reissuedSpeciesId_from_indexInOriGenomeIdVec(int tmpIndex)
	{
		return reissuedSpeciesIdVec[tmpIndex];
	}

	int return_reissuedGenusId_from_indexInOriGenomeIdVec(int tmpIndex)
	{
		return reissuedGenusIdVec[tmpIndex];
	}	

	int return_reissuedPhylumId_from_indexInOriGenomeIdVec(int tmpIndex)
	{
		return reissuedPhylumIdVec[tmpIndex];
	}		

	int return_species_num()
	{
		return speciesIdNamePairVec.size();
	}

	int return_genus_num()
	{
		return genusIdNamePairVec.size();
	}	

	int return_phylum_num()
	{
		return phylumIdNamePairVec.size();
	}

	int return_genome_num()
	{
		return oriGenomeIdVec.size();
	}

	void print_taxoInfoDir_Jflist_reissuedId_oriId_genome_taxo(string& outputDir_taxoInfo)
	{
		this->print_taxoInfoDir_KmerSet_reissuedId_oriId_genome_taxo(outputDir_taxoInfo);
	}

	void print_taxoInfoDir_KmerSet_reissuedId_oriId_genome_taxo(string& outputDir_taxoInfo)
	{
		string mkdir_outputDir_taxoInfo = "mkdir " + outputDir_taxoInfo;
		system(mkdir_outputDir_taxoInfo.c_str());
		//string outputFile_KmerSet_reissuedId_oriId_genome_taxo = outputDir_taxoInfo + "KmerSet_reissuedId_oriId_genome_taxo.txt";
		//this->print_KmerSet_reissuedId_oriId_genome_taxo(inputKmer2taxoInfoFile, outputFile_KmerSet_reissuedId_oriId_genome_taxo);
		string outputFile_reissuedId_to_oriId_genome = outputDir_taxoInfo + "reissuedId_to_oriId.genome.txt";
		string outputFile_reissuedId_to_oriId_species = outputDir_taxoInfo + "reissuedId_to_oriId.species.txt";
		string outputFile_reissuedId_to_oriId_genus = outputDir_taxoInfo + "reissuedId_to_oriId.genus.txt";
		string outputFile_reissuedId_to_oriId_phylum = outputDir_taxoInfo + "reissuedId_to_oriId.phylum.txt";
		this->print_reissuedId_to_oriId_genome(outputFile_reissuedId_to_oriId_genome);
		this->print_reissuedId_to_oriId_species(outputFile_reissuedId_to_oriId_species);
		this->print_reissuedId_to_oriId_genus(outputFile_reissuedId_to_oriId_genus);
		this->print_reissuedId_to_oriId_phylum(outputFile_reissuedId_to_oriId_phylum);
	}

	// void print_KmerSet_reissuedId_oriId_genome_taxo(string& inputKmer2taxoInfoFile, string& outputFile_KmerSet_reissuedId_oriId_genome_taxo)
	// {
	// 	ofstream taxoInfo_ofs(outputFile_KmerSet_reissuedId_oriId_genome_taxo.c_str());
	// 	ifstream Kmer2taxoInfo_ifs(Kmer2taxoInfoFile.c_str());
	// 	taxoInfo_ofs << "KmerSetFile\tGenome_ID_reissued\tGenome_ID_ori\t";
	// 	taxoInfo_ofs << "Species_ID_reissued\tSpecies_ID_ori\tSpecies_Name\t";
	// 	taxoInfo_ofs << "Genus_ID_reissued\tGenus_ID_ori\tGenus_Name\t";
	// 	taxoInfo_ofs << "Phylum_ID_reissued\tPhylum_ID_ori\tPhylum_Name" << endl;
	// 	int tmpReissuedId = 0;
	// 	while(!Kmer2taxoInfo_ifs.eof())
	// 	{
	// 		string tmpStr;
	// 		getline(Kmer2taxoInfo_ifs, tmpStr);
	// 		if(tmpStr == "")
	// 			break;
	// 		vector<string> tmpFieldVec;
	// 		parseStr2fieldVec(tmpFieldVec, tmpStr);
	// 		string tmpKmerSetFile = tmpFieldVec[0];
	// 		taxoInfo_ofs << tmpKmerSetFile << "\t" << tmpReissuedId << "\t" << oriGenomeIdVec[tmpReissuedId] << "\t";
	// 		taxoInfo_ofs << 


	// 		int tmpOriGenomeId = atoi(tmpFieldVec[1].c_str());
	// 		int tmpSpeciesId = atoi(tmpFieldVec[2].c_str());
	// 		int tmpGenusId = atoi(tmpFieldVec[3].c_str());
	// 		int tmpPhylumId = atoi(tmpFieldVec[7].c_str());						
	// 		oriGenomeIdVec.push_back(tmpOriGenomeId);
	// 		speciesIdVec.push_back(tmpSpeciesId);
	// 		genusIdVec.push_back(tmpGenusId);
	// 		phylumIdVec.push_back(tmpPhylumId);		
	// 		tmpReissuedId ++;
	// 	}
	// 	Kmer2taxoInfo_ifs.close();
	// 	taxoInfo_ofs.close();
	// }

	void print_reissuedId_to_oriId_genome(string& tmp_reissuedId_to_oriId_genome_file)
	{
		ofstream genome_ofs(tmp_reissuedId_to_oriId_genome_file.c_str());
		genome_ofs << "Reissued_index\tGenome_ID" << endl;
		for(int tmp = 0; tmp < oriGenomeIdVec.size(); tmp++)
			genome_ofs << tmp << "\t" << oriGenomeIdVec[tmp] << endl;
		genome_ofs.close();
	}

	void print_reissuedId_to_oriId_species(string& tmp_reissuedId_to_oriId_species_file)
	{
		ofstream species_ofs(tmp_reissuedId_to_oriId_species_file.c_str());
		species_ofs << "Reissued_index\tSpecies_ID\tSpecies_name" << endl;
		for(int tmp = 0; tmp < speciesIdNamePairVec.size(); tmp++)
			species_ofs << tmp << "\t" << speciesIdNamePairVec[tmp].first << "\t" << speciesIdNamePairVec[tmp].second << endl;
		species_ofs.close();
	}	

	void print_reissuedId_to_oriId_genus(string& tmp_reissuedId_to_oriId_genus_file)
	{
		ofstream genus_ofs(tmp_reissuedId_to_oriId_genus_file.c_str());
		genus_ofs << "Reissued_index\tGenus_ID\tGenus_name" << endl;
		for(int tmp = 0; tmp < genusIdNamePairVec.size(); tmp++)
			genus_ofs << tmp << "\t" << genusIdNamePairVec[tmp].first << "\t" << genusIdNamePairVec[tmp].second << endl;
		genus_ofs.close();	
	}	

	void print_reissuedId_to_oriId_phylum(string& tmp_reissuedId_to_oriId_phylum_file)
	{
		ofstream phylum_ofs(tmp_reissuedId_to_oriId_phylum_file.c_str());
		phylum_ofs << "Reissued_index\tPhylum_ID\tPhylum_name" << endl;
		for(int tmp = 0; tmp < phylumIdNamePairVec.size(); tmp++)
			phylum_ofs << tmp << "\t" << phylumIdNamePairVec[tmp].first << "\t" << phylumIdNamePairVec[tmp].second << endl;
		phylum_ofs.close();
	}	

	void print_reissuedGenomeId2taxoIdName(string& tmpFile, string& tmpSpeciesFile, string& tmpGenusFile, string& tmpPhylumFile)
	{
		ofstream out_ofs(tmpFile.c_str());
		out_ofs << "Total genome #: " << oriGenomeIdVec.size() << endl;
		out_ofs << "ReissuedGenomeIndex\tOriGenomeId\tSpeciesID\tSpeciesName\tGenusID\tGenusName\tPhylumID\tPhylumNam" << endl;
		for(int tmp = 0; tmp < oriGenomeIdVec.size(); tmp++)
			out_ofs << tmp << "\t" << oriGenomeIdVec[tmp] << "\t" 
				<< speciesIdVec[tmp] << "\t" << speciesNameVec[tmp] << "\t"
				<< genusIdVec[tmp] << "\t" << genusNameVec[tmp] << "\t"
				<< phylumIdVec[tmp] << "\t" << phylumNameVec[tmp] << endl;
		out_ofs.close();

		ofstream species_ofs(tmpSpeciesFile.c_str());
		species_ofs << "Reissued_index\tSpecies_ID\tSpecies_name" << endl;
		for(int tmp = 0; tmp < speciesIdNamePairVec.size(); tmp++)
			species_ofs << tmp << "\t" << speciesIdNamePairVec[tmp].first << "\t" << speciesIdNamePairVec[tmp].second << endl;
		species_ofs.close();

		ofstream genus_ofs(tmpGenusFile.c_str());
		genus_ofs << "Reissued_index\tGenus_ID\tGenus_name" << endl;
		for(int tmp = 0; tmp < genusIdNamePairVec.size(); tmp++)
			genus_ofs << tmp << "\t" << genusIdNamePairVec[tmp].first << "\t" << genusIdNamePairVec[tmp].second << endl;
		genus_ofs.close();		

		ofstream phylum_ofs(tmpPhylumFile.c_str());
		phylum_ofs << "Reissued_index\tPhylum_ID\tPhylum_name" << endl;
		for(int tmp = 0; tmp < phylumIdNamePairVec.size(); tmp++)
			phylum_ofs << tmp << "\t" << phylumIdNamePairVec[tmp].first << "\t" << phylumIdNamePairVec[tmp].second << endl;
		phylum_ofs.close();		
	}

	int returnTaxoId_fromReissuedId_species(int tmpReissuedId)
	{
		return speciesIdVec[tmpReissuedId - 1];
	}

	int returnTaxoId_fromReissuedId_genus(int tmpReissuedId)
	{
		return genusIdVec[tmpReissuedId - 1];
	}

	int returnTaxoId_fromReissuedId_phylum(int tmpReissuedId)
	{
		return phylumIdVec[tmpReissuedId - 1];
	}

	string returnTaxoName_fromReissuedId_species(int tmpReissuedId)
	{
		return speciesNameVec[tmpReissuedId - 1];
	}

	string returnTaxoName_fromReissuedId_genus(int tmpReissuedId)
	{
		return genusNameVec[tmpReissuedId - 1];
	}

	string returnTaxoName_fromReissuedId_phylum(int tmpReissuedId)
	{
		return phylumNameVec[tmpReissuedId - 1];
	}

	int returnTaxoId_fromReissuedId(int tmpReissuedId, int tmpRank)
	//1-domain, 2-kindom, 3-phylum, 4-class, 
	//5-order, 6-family, 7-genus, 8-species
	{
		if(tmpRank == 3)
			return this->returnTaxoId_fromReissuedId_phylum(tmpReissuedId);
		else if(tmpRank = 7)
			return this->returnTaxoId_fromReissuedId_genus(tmpReissuedId);
		else if(tmpRank = 8)
			return this->returnTaxoId_fromReissuedId_species(tmpReissuedId);		
		else
		{
			cout << "invalid rank: " << tmpRank << endl;
			exit(1);
		}
	}

	string returnTaxoName_fromReissuedId(int tmpReissuedId, int tmpRank)
	//1-domain, 2-kindom, 3-phylum, 4-class, 
	//5-order, 6-family, 7-genus, 8-species
	{
		if(tmpRank == 3)
			return this->returnTaxoName_fromReissuedId_phylum(tmpReissuedId);
		else if(tmpRank = 7)
			return this->returnTaxoName_fromReissuedId_genus(tmpReissuedId);
		else if(tmpRank = 8)
			return this->returnTaxoName_fromReissuedId_species(tmpReissuedId);		
		else
		{
			cout << "invalid rank: " << tmpRank << endl;
			exit(1);
		}
	}

	void initiate_Jf2taxoInfoFile_NCBIfullTaxoId2NameFile(
		string& Jf2taxoInfoFile, string& NCBIfullTaxoId2NameFile)
	{
		//cout << "start to initiate_Jf2taxoInfoFile_NCBIfullTaxoId2NameFile(" << endl;
		this->initiate_Kmer2taxoInfoFile_NCBIfullTaxoId2NameFile(
			Jf2taxoInfoFile, NCBIfullTaxoId2NameFile);
	}

	void initiate_Kmer2taxoInfoFile_NCBIfullTaxoId2NameFile(
		string& Kmer2taxoInfoFile, string& NCBIfullTaxoId2NameFile)
	{
		cout << "start to initiate_Kmer2taxoInfoFile_NCBIfullTaxoId2NameFile" << endl;
		ifstream Kmer2taxoInfo_ifs(Kmer2taxoInfoFile.c_str());
		int tmpReissuedId = 0;
		while(!Kmer2taxoInfo_ifs.eof())
		{
			string tmpStr;
			getline(Kmer2taxoInfo_ifs, tmpStr);
			if(tmpStr == "")
				break;
			vector<string> tmpFieldVec;
			parseStr2fieldVec(tmpFieldVec, tmpStr);
			int tmpOriGenomeId = atoi(tmpFieldVec[1].c_str());
			int tmpSpeciesId = atoi(tmpFieldVec[2].c_str());
			int tmpGenusId = atoi(tmpFieldVec[3].c_str());
			int tmpPhylumId = atoi(tmpFieldVec[7].c_str());
			cout << "tmpOriGenomeId: " << tmpOriGenomeId << endl;
			cout << "tmpSpeciesId: " << tmpSpeciesId << endl;
			cout << "tmpGenusId: " << tmpGenusId << endl;
			cout << "tmpPhylumId: " << tmpPhylumId << endl;					
			oriGenomeIdVec.push_back(tmpOriGenomeId);
			speciesIdVec.push_back(tmpSpeciesId);
			genusIdVec.push_back(tmpGenusId);
			phylumIdVec.push_back(tmpPhylumId);		
			tmpReissuedId ++;
		}
		Kmer2taxoInfo_ifs.close();
		cout << "start to initiate NCBIfullTaxoId2Name" << endl;
		this->initiateNameVec_idVec_NCBIfullTaxoId2NameFile(NCBIfullTaxoId2NameFile);
	}

	void initiate_reissuedId2oriGenomeIdFile_oriGenomeId2taxoIdFile_NCBIfullTaxoId2NameFile(
		string& reissuedId2oriGenomeIdFile, 
		string& oriGenomeid2taxoIdFile,
		string& NCBIfullTaxoId2NameFile)
	{
		cout << "start to initiate oriGenomeid2taxoId" << endl;
		ifstream oriGenomeid2taxoId_ifs(oriGenomeid2taxoIdFile.c_str());
		vector<int> interOriGenomeIdVec;
		vector<int> interSpeciesIdVec;
		vector<int> interGenusIdVec;
		vector<int> interPhylumIdVec;
		while(!oriGenomeid2taxoId_ifs.eof())
		{
			string tmpStr;
			getline(oriGenomeid2taxoId_ifs, tmpStr);
			if(tmpStr == "")
				break;
			vector<string> tmpFieldVec;
			parseStr2fieldVec(tmpFieldVec, tmpStr);
			int tmpOriGenomeId = atoi(tmpFieldVec[0].c_str());
			int tmpSpeciesId = atoi(tmpFieldVec[1].c_str());
			int tmpGenusId = atoi(tmpFieldVec[2].c_str());
			int tmpPhylumId = atoi(tmpFieldVec[6].c_str());
			interOriGenomeIdVec.push_back(tmpOriGenomeId);
			interSpeciesIdVec.push_back(tmpSpeciesId);
			interGenusIdVec.push_back(tmpGenusId);
			interPhylumIdVec.push_back(tmpPhylumId);
		}
		oriGenomeid2taxoId_ifs.close();
		cout << "start to initiate reissuedId2oriGenomeIdFile" << endl;
		int interOriGenomeIdVecSize = interOriGenomeIdVec.size();
		ifstream reissuedId2oriGenomeId_ifs(reissuedId2oriGenomeIdFile.c_str());
		while(!reissuedId2oriGenomeId_ifs.eof())
		{
			string tmpStr;
			getline(reissuedId2oriGenomeId_ifs, tmpStr);
			if(tmpStr == "")
				break;
			int tabLoc = tmpStr.find("\t");
			string tmpReissuedIdStr = tmpStr.substr(0, tabLoc);
			string tmpOriGenomeIdStr = tmpStr.substr(tabLoc + 1);
			int tmpReissuedId = atoi(tmpReissuedIdStr.c_str());
			int tmpOriGenomeId = atoi(tmpOriGenomeIdStr.c_str());
			oriGenomeIdVec.push_back(tmpOriGenomeId);

			bool tmpOriGenomeId_exist_bool = false;
			for(int tmpInterVecIndex = 0; tmpInterVecIndex < interOriGenomeIdVecSize; tmpInterVecIndex ++)
			{
				if(interOriGenomeIdVec[tmpInterVecIndex] == tmpOriGenomeId)
				{
					speciesIdVec.push_back(interSpeciesIdVec[tmpInterVecIndex]);
					genusIdVec.push_back(interGenusIdVec[tmpInterVecIndex]);
					phylumIdVec.push_back(interPhylumIdVec[tmpInterVecIndex]);
					tmpOriGenomeId_exist_bool = true;
					break;
				}
			}
			if(!tmpOriGenomeId_exist_bool)
			{
				cout << "tmpOriGenomeId not found: " << tmpOriGenomeId << endl;
				exit(1);
			}
		}
		reissuedId2oriGenomeId_ifs.close();
		cout << "start to initiate NCBIfullTaxoId2Name" << endl;
		this->initiateNameVec_idVec_NCBIfullTaxoId2NameFile(NCBIfullTaxoId2NameFile);
	}

	void initiateNameVec_idVec_NCBIfullTaxoId2NameFile(string& NCBIfullTaxoId2NameFile)
	{
		NCBIfullTaxoID2Name_Info fullTaxoId2NameInfo;
		cout << "start to initiate fullTaxoId2NameInfo" << endl;
		fullTaxoId2NameInfo.initiate_taxoID2NameFile(NCBIfullTaxoId2NameFile);
		int NCBIfullTaxoIdMax = fullTaxoId2NameInfo.return_NCBIfullTaxoIdMax();
		int genomeIdVecSize = oriGenomeIdVec.size();
		cout << "genomeIdVecSize: " << genomeIdVecSize << endl;
		for(int tmp = 0; tmp < genomeIdVecSize; tmp++)
		{
			//cout << "tmp: " << tmp << endl;
			int tmpSpeciesId = speciesIdVec[tmp];
			int tmpGenusId = genusIdVec[tmp];
			int tmpPhylumId = phylumIdVec[tmp];
			//cout << "tmpSpeciesId: " << tmpSpeciesId << endl;
			//cout << "tmpGenusId: " << tmpGenusId << endl;
			//cout << "tmpPhylumId: " << tmpPhylumId << endl;
			if((tmpSpeciesId < 0)||(tmpSpeciesId > NCBIfullTaxoIdMax))
				speciesNameVec.push_back("NULL");
			else
				speciesNameVec.push_back(fullTaxoId2NameInfo.return_taxoName(tmpSpeciesId));

			if((tmpGenusId < 0)||(tmpGenusId > NCBIfullTaxoIdMax))
				genusNameVec.push_back("NULL");
			else
				genusNameVec.push_back(fullTaxoId2NameInfo.return_taxoName(tmpGenusId));

			if((tmpPhylumId < 0)||(tmpPhylumId > NCBIfullTaxoIdMax))
				phylumNameVec.push_back("NULL");
			else
				phylumNameVec.push_back(fullTaxoId2NameInfo.return_taxoName(tmpPhylumId));			
		}
	}

	void parseStr2fieldVec(vector<string>& tmpFieldVec, string& tmpStr)
	{
		int startLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tabLoc = tmpStr.find("\t", startLoc);
			if(tabLoc == string::npos)
				break;
			string tmpField = tmpStr.substr(startLoc, tabLoc-startLoc);
			tmpFieldVec.push_back(tmpField);
			startLoc = tabLoc + 1;
		}
		tmpFieldVec.push_back(tmpStr.substr(startLoc));
	}
};
#endif