#ifndef CHROMOSOMESEQ_INFO_H
#define CHROMOSOMESEQ_IFNO_H
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
#include "taxonomy_info.h"
using namespace std;

class ChromosomeSeq_Info
{
private:
	vector<string> chromosomeSeqFilePathVec;
	int accession_id;
	
	bool species_assignedOrNot_bool;
	int species_id;
	string species_name;

	bool genus_assignedOrNot_bool;
	int genus_id;
	string genus_name;	

	bool family_assignedOrNot_bool;
	int family_id;
	string family_name;

	bool order_assignedOrNot_bool;
	int order_id;
	string order_name;

	bool class_assignedOrNot_bool;
	int class_id;
	string class_name;

	bool phylum_assignedOrNot_bool;
	int phylum_id;
	string phylum_name;

	////assigned fa, jf, Kmer, sortedKmer files
	string mergedGenome_fa_path;
	vector< pair< int, string> > mergedGenome_jf_pathVec;
	vector< pair< int, string> > mergedGenome_Kmer_pathVec;
	vector< pair< int, string> > mergedGenome_sortedKmer_pathVec;
public:
	ChromosomeSeq_Info()
	{
		accession_id = -1;
	
		species_assignedOrNot_bool = false;
		species_id = -1;
		species_name = "UNKNOWN";

		genus_assignedOrNot_bool = false;
		genus_id = -1;
		genus_name = "UNKNOWN";

		family_assignedOrNot_bool = false;
		family_id = -1;
		family_name = "UNKNOWN";

		order_assignedOrNot_bool = false;
		order_id = -1;
		order_name = "UNKNOWN";

		class_assignedOrNot_bool = false;
		class_id = -1;
		class_name = "UNKNOWN";

		phylum_assignedOrNot_bool = false;
		phylum_id = -1;
		phylum_name = "UNKNOWN";				
	}

	string return_mergedGenome_jf_path_with_KmerLength(int tmp_Kmer_length)
	{
		for(int tmp = 0; tmp < mergedGenome_jf_pathVec.size(); tmp++)
		{
			if(mergedGenome_jf_pathVec[tmp].first == tmp_Kmer_length)
				return mergedGenome_jf_pathVec[tmp].second;
		}
		cout << "no jf file was found for Kmer length: " << tmp_Kmer_length << endl;
		exit(1);
	}

	string return_mergedGenome_Kmer_path_with_KmerLength(int tmp_Kmer_length)
	{
		for(int tmp = 0; tmp < mergedGenome_Kmer_pathVec.size(); tmp++)
		{
			if(mergedGenome_Kmer_pathVec[tmp].first == tmp_Kmer_length)
				return mergedGenome_Kmer_pathVec[tmp].second;
		}
		cout << "no Kmer file was found for Kmer length: " << tmp_Kmer_length << endl;
		exit(1);
	}

	string return_mergedGenome_sortedKmer_path_with_KmerLength(int tmp_Kmer_length)
	{
		for(int tmp = 0; tmp < mergedGenome_sortedKmer_pathVec.size(); tmp++)
		{
			if(mergedGenome_sortedKmer_pathVec[tmp].first == tmp_Kmer_length)
				return mergedGenome_sortedKmer_pathVec[tmp].second;
		}
		cout << "no sortedKmer file was found for Kmer length: " << tmp_Kmer_length << endl;
		exit(1);
	}	

	int return_accession_id()
	{
		return accession_id;
	}

	int return_species_id()
	{
		return species_id;
	}

	int return_genus_id()
	{
		return genus_id;
	}

	int return_family_id()
	{
		return family_id;
	}		

	int return_order_id()
	{
		return order_id;
	}	

	int return_class_id()
	{
		return class_id;
	}		

	int return_phylum_id()
	{
		return phylum_id;
	}	

	string return_species_name()
	{
		return species_name;
	}

	string return_genus_name()
	{
		return genus_name;
	}

	string return_family_name()
	{
		return family_name;
	}		

	string return_order_name()
	{
		return order_name;
	}	

	string return_class_name()
	{
		return class_name;
	}		

	string return_phylum_name()
	{
		return phylum_name;
	}	

	string return_mergedGenome_fa_path()
	{
		return mergedGenome_fa_path;
	}

	void assign_mergedGenome_fa_path(string& tmp_mergedGenome_fa_path)
	{
		mergedGenome_fa_path = tmp_mergedGenome_fa_path;
	}

	void assign_mergedGenome_jf_path(string& tmp_mergedGenome_jf_path, int tmp_Kmer_length)
	{
		mergedGenome_jf_pathVec.push_back(pair<int,string>(tmp_Kmer_length, tmp_mergedGenome_jf_path));
	}

	void assign_mergedGenome_Kmer_path(string& tmp_mergedGenome_Kmer_path, int tmp_Kmer_length)
	{
		mergedGenome_Kmer_pathVec.push_back(pair<int,string>(tmp_Kmer_length, tmp_mergedGenome_Kmer_path));
	}

	void assign_mergedGenome_sortedKmer_path(string& tmp_mergedGenome_sortedKmer_path, int tmp_Kmer_length)
	{
		mergedGenome_sortedKmer_pathVec.push_back(pair<int,string>(tmp_Kmer_length, tmp_mergedGenome_sortedKmer_path));
	}	

	string return_chromosomeSeqFilePath(int tmp)
	{
		return chromosomeSeqFilePathVec[tmp];
	}

	int return_chromosomeSeqFileNum()
	{
		return chromosomeSeqFilePathVec.size();
	}

	void initiate_1stChromosomeSeq(string& chromosomeSeqFilePath_1st, int tmp_accession_id, 
		string& taxonomyInfoStr, Taxonomy_Info& taxonomyInfo)
	{
		chromosomeSeqFilePathVec.push_back(chromosomeSeqFilePath_1st);
		accession_id = tmp_accession_id;

		int tabLoc_1 = taxonomyInfoStr.find("\t");
		int tabLoc_2 = taxonomyInfoStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = taxonomyInfoStr.find("\t", tabLoc_2 + 1);
		int tabLoc_4 = taxonomyInfoStr.find("\t", tabLoc_3 + 1);
		int tabLoc_5 = taxonomyInfoStr.find("\t", tabLoc_4 + 1);

		string species_id_str = taxonomyInfoStr.substr(0, tabLoc_1);
		string genus_id_str = taxonomyInfoStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string family_id_str = taxonomyInfoStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		string order_id_str = taxonomyInfoStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
		string class_id_str = taxonomyInfoStr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
		string phylum_id_str = taxonomyInfoStr.substr(tabLoc_5 + 1);

		if(species_id_str != "UNKNOWN")
		{
			species_assignedOrNot_bool = true;
			species_id = atoi(species_id_str.c_str());
			species_name = taxonomyInfo.return_name_from_id_species(species_id);
		}

		if(genus_id_str != "UNKNOWN")
		{
			genus_assignedOrNot_bool = true;
			genus_id = atoi(genus_id_str.c_str());
			genus_name = taxonomyInfo.return_name_from_id_genus(genus_id);
		}

		if(family_id_str != "UNKNOWN")
		{
			family_assignedOrNot_bool = true;
			family_id = atoi(family_id_str.c_str());
			family_name = taxonomyInfo.return_name_from_id_family(family_id);
		}

		if(order_id_str != "UNKNOWN")
		{
			order_assignedOrNot_bool = true;
			order_id = atoi(order_id_str.c_str());
			order_name = taxonomyInfo.return_name_from_id_order(order_id);
		}

		if(class_id_str != "UNKNOWN")
		{
			class_assignedOrNot_bool = true;
			class_id = atoi(class_id_str.c_str());
			class_name = taxonomyInfo.return_name_from_id_class(class_id);
		}								

		if(phylum_id_str != "UNKNOWN")
		{
			phylum_assignedOrNot_bool = true;
			phylum_id = atoi(phylum_id_str.c_str());
			phylum_name = taxonomyInfo.return_name_from_id_phylum(phylum_id);
		}			
	}

	void initiate_1stChromosomeSeq(string& chromosomeSeqFilePath_1st, int tmp_accession_id, 
		string& taxonomyInfoStr)//, Taxonomy_Info& taxonomyInfo)
	{
		chromosomeSeqFilePathVec.push_back(chromosomeSeqFilePath_1st);
		accession_id = tmp_accession_id;

		int tabLoc_1 = taxonomyInfoStr.find("\t");
		int tabLoc_2 = taxonomyInfoStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = taxonomyInfoStr.find("\t", tabLoc_2 + 1);
		int tabLoc_4 = taxonomyInfoStr.find("\t", tabLoc_3 + 1);
		int tabLoc_5 = taxonomyInfoStr.find("\t", tabLoc_4 + 1);

		string species_id_str = taxonomyInfoStr.substr(0, tabLoc_1);
		string genus_id_str = taxonomyInfoStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string family_id_str = taxonomyInfoStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		string order_id_str = taxonomyInfoStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
		string class_id_str = taxonomyInfoStr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
		string phylum_id_str = taxonomyInfoStr.substr(tabLoc_5 + 1);

		if(species_id_str != "UNKNOWN")
		{
			species_assignedOrNot_bool = true;
			species_id = atoi(species_id_str.c_str());
			species_name = species_id_str; //taxonomyInfo.return_name_from_id_species(species_id);
		}

		if(genus_id_str != "UNKNOWN")
		{
			genus_assignedOrNot_bool = true;
			genus_id = atoi(genus_id_str.c_str());
			genus_name = genus_id_str; // taxonomyInfo.return_name_from_id_genus(genus_id);
		}

		if(family_id_str != "UNKNOWN")
		{
			family_assignedOrNot_bool = true;
			family_id = atoi(family_id_str.c_str());
			family_name = family_id_str; // taxonomyInfo.return_name_from_id_family(family_id);
		}

		if(order_id_str != "UNKNOWN")
		{
			order_assignedOrNot_bool = true;
			order_id = atoi(order_id_str.c_str());
			order_name = order_id_str; // taxonomyInfo.return_name_from_id_order(order_id);
		}

		if(class_id_str != "UNKNOWN")
		{
			class_assignedOrNot_bool = true;
			class_id = atoi(class_id_str.c_str());
			class_name = class_id_str; // taxonomyInfo.return_name_from_id_class(class_id);
		}								

		if(phylum_id_str != "UNKNOWN")
		{
			phylum_assignedOrNot_bool = true;
			phylum_id = atoi(phylum_id_str.c_str());
			phylum_name = phylum_id_str; // taxonomyInfo.return_name_from_id_phylum(phylum_id);
		}			
	}

	void update_newChromosomeSeq(string& chromosomeSeqFilePath_new)
	{
		chromosomeSeqFilePathVec.push_back(chromosomeSeqFilePath_new);
	}

};
#endif