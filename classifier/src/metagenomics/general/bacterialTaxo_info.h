// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#ifndef BACTERIALTAXO_INFO_H
#define BACTERIALTAXO_IFNO_H
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
using namespace std;

class BacterialTaxo_Info
{
private:
	vector<int> speciesIdVec;
	vector<int> genusIdVec;
	vector<int> familyIdVec;
	vector<int> orderIdVec;
	vector<int> classIdVec;
	vector<int> phylumIdVec;	

	vector<string> speciesNameVec;
	vector<string> genusNameVec;
	vector<string> familyNameVec;
	vector<string> orderNameVec;
	vector<string> classNameVec;
	vector<string> phylumNameVec;	

	vector<int> reissuedSpeciesIdVec;
	vector<int> reissuedGenusIdVec;
	vector<int> reissuedFamilyIdVec;
	vector<int> reissuedOrderIdVec;
	vector<int> reissuedClassIdVec;
	vector<int> reissuedPhylumIdVec;	

	// reissued id 2 taxoId and name
	vector< pair<int, string> > speciesIdNamePairVec;
	vector< pair<int, string> > genusIdNamePairVec;
	vector< pair<int, string> > familyIdNamePairVec;
	vector< pair<int, string> > orderIdNamePairVec;
	vector< pair<int, string> > classIdNamePairVec;
	vector< pair<int, string> > phylumIdNamePairVec;		

	// total taxo num in each level
	int species_num;
	int genus_num;
	int family_num;
	int order_num;
	int class_num;
	int phylum_num;

	int* speciesReissuedId_2_higherRankReissuedId_array;
	int* genusReissuedId_2_higherRankReissuedId_array;
	int* familyReissuedId_2_higherRankReissuedId_array;
	int* orderReissuedId_2_higherRankReissuedId_array;
	int* classReissuedId_2_higherRankReissuedId_array;
public:
	BacterialTaxo_Info()
	{}

	int get_rawTaxoId_from_reissuedTaxoId(int tmpReissuedId, int tmpRank)
	{
		if(tmpRank == 3) // phylum
			return reissuedPhylumIdVec[tmpReissuedId];
		else if(tmpRank == 4)
			return reissuedClassIdVec[tmpReissuedId];
		else if(tmpRank == 5)
			return reissuedOrderIdVec[tmpReissuedId];
		else if(tmpRank == 6)
			return reissuedFamilyIdVec[tmpReissuedId];
		else if(tmpRank == 7)
			return reissuedGenusIdVec[tmpReissuedId];
		else if(tmpRank == 8)
			return reissuedSpeciesIdVec[tmpReissuedId];
		else
		{
			cout << "invalid tmpRank in get_rawTaxoId_from_reissuedTaxoId, tmpRank: " << tmpRank << endl;
			exit(1);
		}
	}

	void initiate_lowerRankTaxoReissuedId_to_higherRankTaxoReissuedId_array()
	{
		speciesReissuedId_2_higherRankReissuedId_array = new int[5 * species_num];
		for(int tmp = 0; tmp < 5 * species_num; tmp++)
			speciesReissuedId_2_higherRankReissuedId_array[tmp] = -2;
		genusReissuedId_2_higherRankReissuedId_array = new int[4 * genus_num];
		for(int tmp = 0; tmp < 4 * genus_num; tmp++)
			genusReissuedId_2_higherRankReissuedId_array[tmp] = -2;
		familyReissuedId_2_higherRankReissuedId_array = new int[3 * family_num];
		for(int tmp = 0; tmp < 3 * family_num; tmp++)
			familyReissuedId_2_higherRankReissuedId_array[tmp] = -2;
		orderReissuedId_2_higherRankReissuedId_array = new int[2 * order_num];
		for(int tmp = 0; tmp < 2 * order_num; tmp++)
			orderReissuedId_2_higherRankReissuedId_array[tmp] = -2;
		classReissuedId_2_higherRankReissuedId_array = new int[class_num];
		for(int tmp = 0; tmp < class_num; tmp++)
			classReissuedId_2_higherRankReissuedId_array[tmp] = -2;
		for(int tmpIndex = 0; tmpIndex < reissuedSpeciesIdVec.size(); tmpIndex ++)
		{
			int tmpReissuedId_species = reissuedSpeciesIdVec[tmpIndex];
			int tmpReissuedId_genus = reissuedGenusIdVec[tmpIndex];
			int tmpReissuedId_family = reissuedFamilyIdVec[tmpIndex];
			int tmpReissuedId_order = reissuedOrderIdVec[tmpIndex];
			int tmpReissuedId_class = reissuedClassIdVec[tmpIndex];
			int tmpReissuedId_phylum = reissuedPhylumIdVec[tmpIndex];

			// species array
			speciesReissuedId_2_higherRankReissuedId_array[0 * species_num + tmpReissuedId_species] = tmpReissuedId_genus;
			speciesReissuedId_2_higherRankReissuedId_array[1 * species_num + tmpReissuedId_species] = tmpReissuedId_family;
			speciesReissuedId_2_higherRankReissuedId_array[2 * species_num + tmpReissuedId_species] = tmpReissuedId_order;
			speciesReissuedId_2_higherRankReissuedId_array[3 * species_num + tmpReissuedId_species] = tmpReissuedId_class;
			speciesReissuedId_2_higherRankReissuedId_array[4 * species_num + tmpReissuedId_species] = tmpReissuedId_phylum;

			// genus array
			genusReissuedId_2_higherRankReissuedId_array[0 * genus_num + tmpReissuedId_genus] = tmpReissuedId_family;
			genusReissuedId_2_higherRankReissuedId_array[1 * genus_num + tmpReissuedId_genus] = tmpReissuedId_order;
			genusReissuedId_2_higherRankReissuedId_array[2 * genus_num + tmpReissuedId_genus] = tmpReissuedId_class;
			genusReissuedId_2_higherRankReissuedId_array[3 * genus_num + tmpReissuedId_genus] = tmpReissuedId_phylum;

			// family array
			familyReissuedId_2_higherRankReissuedId_array[0 * family_num + tmpReissuedId_family] = tmpReissuedId_order;
			familyReissuedId_2_higherRankReissuedId_array[1 * family_num + tmpReissuedId_family] = tmpReissuedId_class;
			familyReissuedId_2_higherRankReissuedId_array[2 * family_num + tmpReissuedId_family] = tmpReissuedId_phylum;

			// order array
			orderReissuedId_2_higherRankReissuedId_array[0 * order_num + tmpReissuedId_order] = tmpReissuedId_class;
			orderReissuedId_2_higherRankReissuedId_array[1 * order_num + tmpReissuedId_order] = tmpReissuedId_phylum;
			
			// class array
			classReissuedId_2_higherRankReissuedId_array[0 * class_num + tmpReissuedId_class] = tmpReissuedId_phylum;
		}
	}

	int get_higherRankTaxoReissuedId_from_lowerRankTaxoReissuedId(
		//int& tmp_output_ReissuedId_higherRank, 
		int tmp_input_TaxoRank_higherRank, 
		int tmp_input_ReissuedId_lowerRank, int tmp_input_TaxoRank_lowerRank)
	{
		int rankDiff = tmp_input_TaxoRank_lowerRank - tmp_input_TaxoRank_higherRank;
		switch(tmp_input_TaxoRank_lowerRank)
		{
			case 8: // species
				return speciesReissuedId_2_higherRankReissuedId_array[(rankDiff-1) * species_num + tmp_input_ReissuedId_lowerRank];
			case 7: // genus
				return genusReissuedId_2_higherRankReissuedId_array[(rankDiff-1) * genus_num + tmp_input_ReissuedId_lowerRank];
			case 6: // family
				return familyReissuedId_2_higherRankReissuedId_array[(rankDiff-1) * family_num + tmp_input_ReissuedId_lowerRank];
			case 5: // order
				return orderReissuedId_2_higherRankReissuedId_array[(rankDiff-1) * order_num + tmp_input_ReissuedId_lowerRank];
			case 4: // class
				return classReissuedId_2_higherRankReissuedId_array[(rankDiff-1) * class_num + tmp_input_ReissuedId_lowerRank];
			default:
				cout << "invalid rank in get_higherRankTaxoReissuedId_from_lowerRankTaxoReissuedId, tmp_input_TaxoRank_lowerRank: " 
					<< tmp_input_TaxoRank_lowerRank << endl;
				exit(1);
		}
	}

	void get_reissuedId_allRank(int index, int& tmpReissuedId_species, int& tmpReissuedId_genus,
		int& tmpReissuedId_family, int& tmpReissuedId_order, int& tmpReissuedId_class, int& tmpReissuedId_phylum)
	{
		tmpReissuedId_species = reissuedSpeciesIdVec[index];
		tmpReissuedId_genus = reissuedGenusIdVec[index];
		tmpReissuedId_family = reissuedFamilyIdVec[index];
		tmpReissuedId_order = reissuedOrderIdVec[index];
		tmpReissuedId_class = reissuedClassIdVec[index];
		tmpReissuedId_phylum = reissuedPhylumIdVec[index];
	}

	void print2Vec_species_genus_family_order_class_phylum(
		vector<int>& tmpSpeciesIdVec, vector<int>& tmpGenusIdVec, vector<int>& tmpFamilyIdVec, 
		vector<int>& tmpOrderIdVec, vector<int>& tmpClassIdVec, vector<int>& tmpPhylumIdVec)
	{
		for(int tmp = 0; tmp < speciesIdVec.size(); tmp++)
		{
			tmpSpeciesIdVec.push_back(speciesIdVec[tmp]);
			tmpGenusIdVec.push_back(genusIdVec[tmp]);
			tmpFamilyIdVec.push_back(familyIdVec[tmp]);
			tmpOrderIdVec.push_back(orderIdVec[tmp]);
			tmpClassIdVec.push_back(classIdVec[tmp]);
			tmpPhylumIdVec.push_back(phylumIdVec[tmp]);
		}
	}

	string return_taxo_name(int tmpRank, int tmp_reissued_taxo_id)
	{
		if(tmp_reissued_taxo_id < 0)
			return "INVALID";
		else
		{}
		if(tmpRank == 3) // phylum
			return phylumIdNamePairVec[tmp_reissued_taxo_id].second;
		else if(tmpRank == 4) // class
			return classIdNamePairVec[tmp_reissued_taxo_id].second;
		else if(tmpRank == 5) // order
			return orderIdNamePairVec[tmp_reissued_taxo_id].second;
		else if(tmpRank == 6) // family
			return familyIdNamePairVec[tmp_reissued_taxo_id].second;
		else if(tmpRank == 7) // genus
			return genusIdNamePairVec[tmp_reissued_taxo_id].second;
		else if(tmpRank == 8) // species
			return speciesIdNamePairVec[tmp_reissued_taxo_id].second;
		else
		{
			cout << "error! invalid rank: " << tmpRank << endl;
			exit(1);
		}			
	}

	int return_taxo_Id(int tmpRank, int tmp_reissued_taxo_id)
	{
		if(tmp_reissued_taxo_id < 0)
			return tmp_reissued_taxo_id;
		else
		{}
		if(tmpRank == 3) // phylum
			return phylumIdNamePairVec[tmp_reissued_taxo_id].first;
		else if(tmpRank == 4) // class
			return classIdNamePairVec[tmp_reissued_taxo_id].first;
		else if(tmpRank == 5) // order
			return orderIdNamePairVec[tmp_reissued_taxo_id].first;
		else if(tmpRank == 6) // family
			return familyIdNamePairVec[tmp_reissued_taxo_id].first;
		else if(tmpRank == 7) // genus
			return genusIdNamePairVec[tmp_reissued_taxo_id].first;
		else if(tmpRank == 8) // species
			return speciesIdNamePairVec[tmp_reissued_taxo_id].first;
		else
		{
			cout << "error! invalid rank: " << tmpRank << endl;
			exit(1);
		}			
	}

	void get_taxoIdVec_from_taxoReissuedIdVec(vector<int>& tmpTaxoIdVec, vector<int>& tmpTaxoReissuedIdVec, int specific_rank)
	{
		for(int tmp = 0; tmp < tmpTaxoReissuedIdVec.size(); tmp++)
		{
			int tmpReissuedId = tmpTaxoReissuedIdVec[tmp];
			int tmpRawId = this->return_taxo_Id(specific_rank, tmpReissuedId);
			tmpTaxoIdVec.push_back(tmpRawId);
		}
	}

	int return_taxo_reissued_id(int tmpRank, int tmp_reissued_species_id)
	{
		if(tmpRank == 3) // phylum
			return reissuedPhylumIdVec[tmp_reissued_species_id];
		else if(tmpRank == 4) // class
			return reissuedClassIdVec[tmp_reissued_species_id];
		else if(tmpRank == 5) // order
			return reissuedOrderIdVec[tmp_reissued_species_id];
		else if(tmpRank == 6) // family
			return reissuedFamilyIdVec[tmp_reissued_species_id];
		else if(tmpRank == 7) // genus
			return reissuedGenusIdVec[tmp_reissued_species_id];
		else if(tmpRank == 8) // species
			return reissuedSpeciesIdVec[tmp_reissued_species_id];
		else
		{
			cout << "error! invalid rank: " << tmpRank << endl;
			exit(1);
		}			
	}

	int return_taxo_num(int tmpRank)
	{
		if(tmpRank == 3) // phylum
			return phylum_num;
		else if(tmpRank == 4) // class
			return class_num;
		else if(tmpRank == 5) // order
			return order_num;
		else if(tmpRank == 6) // family
			return family_num;
		else if(tmpRank == 7) // genus
			return genus_num;
		else if(tmpRank == 8) // species
			return species_num;
		else
		{
			cout << "error! invalid rank: " << tmpRank << endl;
			exit(1);
		}				
	}

	int return_taxo_num_total()
	{
		return species_num + genus_num + family_num + order_num + class_num + phylum_num;
	}

	int return_species_id(int tmpIndex)
	{
		return speciesIdVec[tmpIndex];
	}

	int return_species_num()
	{
		return species_num;
	}

	int return_genus_num()
	{
		return genus_num;
	}

	int return_family_num()
	{
		return family_num;
	}		

	int return_order_num()
	{
		return order_num;
	}

	int return_class_num()
	{
		return class_num;
	}	

	int return_phylum_num()
	{
		return phylum_num;
	}	

	int returnReissuedId_specificTaxoLevel_from_reissuedSpeciesId(int rank, int reissuedSpeciesId)
	{
		if(rank == 3) // phylum
			return reissuedPhylumIdVec[reissuedSpeciesId];
		else if(rank == 4) // class
			return reissuedClassIdVec[reissuedSpeciesId];
		else if(rank == 5) // order
			return reissuedOrderIdVec[reissuedSpeciesId];
		else if(rank == 6) // family
			return reissuedFamilyIdVec[reissuedSpeciesId];
		else if(rank == 7) // genus
			return reissuedGenusIdVec[reissuedSpeciesId];
		else if(rank == 8) // species
			return reissuedSpeciesId;
		else
		{
			cout << "error! invalid rank: " << rank << endl;
			exit(1);
		}			
	}

	// int get_lothelloKmerSetId_specificTaxoLevel_from_reissuedSpeciesId(int rank, int reissuedSpeciesId)
	// {
	// 	if(rank == 3) // phylum
	// 	{
	// 		int reissuedPhylumId = reissuedPhylumIdVec[reissuedSpeciesId];
	// 		return species_num + genus_num + family_num + order_num + class_num + reissuedPhylumId + 1;
	// 	}
	// 	else if(rank == 4) // class
	// 	{
	// 		int reissuedClassId = reissuedClassIdVec[reissuedSpeciesId];
	// 		return species_num + genus_num + family_num + order_num + reissuedClassId + 1;
	// 	}
	// 	else if(rank == 5) // order
	// 	{
	// 		int reissuedOrderId = reissuedOrderIdVec[reissuedSpeciesId];
	// 		return species_num + genus_num + family_num + reissuedOrderId + 1;
	// 	}
	// 	else if(rank == 6) // family
	// 	{
	// 		int reissuedFamilyId = reissuedFamilyIdVec[reissuedSpeciesId];
	// 		return species_num + genus_num + reissuedFamilyId + 1;
	// 	}
	// 	else if(rank == 7) // genus
	// 	{
	// 		int reissuedGenusId = reissuedGenusIdVec[reissuedSpeciesId];
	// 		return species_num + reissuedGenusId + 1;
	// 	}
	// 	else if(rank == 8) // species
	// 		return reissuedSpeciesId + 1;
	// 	else
	// 	{
	// 		cout << "error! invalid rank: " << rank << endl;
	// 		exit(1);
	// 	}						
	// }

	// int get_reissuedSpecificTaxoId_from_lothelloKmerSetId(int rank, int lothelloKmerSetId)
	// {
	// 	if(rank == 3) // phylum
	// 		return lothelloKmerSetId - (species_num + genus_num + family_num + order_num + class_num + 1);
	// 	else if(rank == 4) // class
	// 		return lothelloKmerSetId - (species_num + genus_num + family_num + order_num + 1);
	// 	else if(rank == 5) // order
	// 		return lothelloKmerSetId - (species_num + genus_num + family_num + 1);
	// 	else if(rank == 6) // family
	// 		return lothelloKmerSetId - (species_num + genus_num + 1);
	// 	else if(rank == 7) // genus
	// 		return lothelloKmerSetId - (species_num + 1);
	// 	else if(rank == 8) // species
	// 		return lothelloKmerSetId - 1;
	// 	else
	// 	{
	// 		cout << "error! invalid rank: " << rank << endl;
	// 		exit(1);
	// 	}			
	// }

	// void get_LCA(int rank_1, int reissuedTaxoId_1, int rank_2, int reissuedTaxoId_2,
	// 	int& rank_LCA, int& reissuedTaxoId_LCA)
	// {
	// 	int rank_high, reissuedTaxoId_high, rank_low, reissuedTaxoId_low; 
	// 	if(rank_1 >= rank_2)
	// 	{
	// 		rank_high = rank_1;
	// 		reissuedTaxoId_high = reissuedTaxoId_1;
	// 		rank_low = rank_2;
	// 		reissuedTaxoId_low = reissuedTaxoId_2;			
	// 	}
	// 	else
	// 	{
	// 		rank_high = rank_2;
	// 		reissuedTaxoId_high = reissuedTaxoId_2;
	// 		rank_low = rank_1;
	// 		reissuedTaxoId_low = reissuedTaxoId_1;			
	// 	}
	// 	for(int tmpRank = rank_high; tmpRank >= 3; tmpRank --)
	// 	{

	// 	}
	// }

	// void update_lothelloKmerSetId_taxoSpecificRepetitiveFlag(
	// 	int& taxoSpecificRepetitiveFlag_update, uint16_t& lothelloKmerSetId_update,
	// 	int taxoSpecificRepetitiveFlag_ori, uint16_t lothelloKmerSetId_ori,
	// 	int toInputReissuedSpeciesId)
	// 	// toReturnTaxoSpecificRepetitiveFlag: 0-not assigned; 3-phylum specific; 4-class specific;
	// 	// 5-order specific; 6-family specific; 7-genus specific; 8--species specific; 
	// 	// 9--genome specific (not applicable now); 10 -- repetitive
	// {
	// 	if(taxoSpecificRepetitiveFlag_ori == 0) 
	// 	{
	// 		taxoSpecificRepetitiveFlag_update = 8;
	// 		lothelloKmerSetId_update = this->get_lothelloKmerSetId_specificTaxoLevel_from_reissuedSpeciesId(
	// 			8, toInputReissuedSpeciesId);
	// 	}
	// 	else if(toInputTaxoSpecificRepetitiveFlag == 3)// phylum specific before
	// 	{
	// 		int tmp_lothelloKmerSetId_phylum = this->get_lothelloKmerSetId_specificTaxoLevel_from_reissuedSpeciesId(
	// 			3, toInputReissuedSpeciesId);
	// 		if(tmp_lothelloKmerSetId_phylum == (int)lothelloKmerSetId_ori)
	// 		{}
	// 		else
	// 			taxoSpecificRepetitiveFlag_update = 10;
	// 	}
	// 	else if(toInputTaxoSpecificRepetitiveFlag == 4)// class specific before
	// 	{}
	// 	else if(toInputTaxoSpecificRepetitiveFlag == 5)// order specific before
	// 	{}
	// 	else if(toInputTaxoSpecificRepetitiveFlag == 6)// family specific before
	// 	{}
	// 	else if(toInputTaxoSpecificRepetitiveFlag == 7)// genus specific before
	// 	{}
	// 	else if(toInputTaxoSpecificRepetitiveFlag == 8)// species specific before
	// 	{
	// 		int tmp_lothelloKmerSetId_species = this->get_lothelloKmerSetId_specificTaxoLevel_from_reissuedSpeciesId(
	// 			8, toInputReissuedSpeciesId);
	// 		if(tmp_lothelloKmerSetId_species == (int)lothelloKmerSetId_ori)
	// 		{}
	// 		else
	// 		{
	// 			int lothelloKmerSetId_genus_ori = 
	// 		}
	// 	}
	// 	else if(toInputTaxoSpecificRepetitiveFlag == 10)// phylum repetitive before
	// 		toReturnTaxoSpecificRepetitiveFlag = 10;
	// 	else
	// 	{
	// 		cout << "error! invalid toInputTaxoSpecificRepetitiveFlag: " << toInputTaxoSpecificRepetitiveFlag << endl;
	// 		exit(1);
	// 	}				
	// }

	void print(string& outputDir)
	{
		outputDir += "/";
		string mkdir_cmd = "mkdir " + outputDir;
		system(mkdir_cmd.c_str());
		string overall_file = outputDir + "overall.txt";
		string species_file = outputDir + "species.txt";
		string genus_file = outputDir + "genus.txt";
		string family_file = outputDir + "family.txt";
		string order_file = outputDir + "order.txt";
		string class_file = outputDir + "class.txt";
		string phylum_file = outputDir + "phylum.txt";
		this->print(overall_file, species_file, genus_file, family_file, 
			order_file, class_file, phylum_file);
	}

	void print(string& overall_file, string& species_file, string& genus_file, 
		string& family_file, string& order_file, string& class_file, string& phylum_file)
	{
		ofstream overall_ofs(overall_file.c_str());
		overall_ofs << "Species_index\tSpecies_ID\tSpecies_name\tGenus_index\tGenus_ID\tGenus_name\t";
		overall_ofs << "Family_index\tFamily_ID\tFamily_name\tOrder_index\tOrder_ID\tOrder_name\t";
		overall_ofs << "Class_index\tClass_ID\tClass_name\tPhylum_index\tPhylum_ID\tPhylum_name" << endl;
		for(int tmp = 0; tmp < speciesIdVec.size(); tmp++)
			overall_ofs << reissuedSpeciesIdVec[tmp] << "\t" << speciesIdVec[tmp] << "\t" << speciesNameVec[tmp] 
				<< "\t" << reissuedGenusIdVec[tmp] << "\t" << genusIdVec[tmp] << "\t" << genusNameVec[tmp]
				<< "\t" << reissuedFamilyIdVec[tmp] << "\t" << familyIdVec[tmp] << "\t" << familyNameVec[tmp]
				<< "\t" << reissuedOrderIdVec[tmp] << "\t" << orderIdVec[tmp] << "\t" << orderNameVec[tmp]
				<< "\t" << reissuedClassIdVec[tmp] << "\t" << classIdVec[tmp] << "\t" << classNameVec[tmp]
				<< "\t" << reissuedPhylumIdVec[tmp] << "\t" << phylumIdVec[tmp] << "\t" << phylumNameVec[tmp] << endl;
		overall_ofs.close();

		ofstream species_ofs(species_file.c_str());
		species_ofs << "Species_index\tSpecies_ID\tSpecies_name" << endl;
		for(int tmp = 0; tmp < speciesIdNamePairVec.size(); tmp++)
			species_ofs << tmp << "\t" << speciesIdNamePairVec[tmp].first << "\t" << speciesIdNamePairVec[tmp].second << endl;
		species_ofs.close();
		ofstream genus_ofs(genus_file.c_str());
		genus_ofs << "Genus_index\tGenus_ID\tGenus_name" << endl;
		for(int tmp = 0; tmp < genusIdNamePairVec.size(); tmp++)
			genus_ofs << tmp << "\t" << genusIdNamePairVec[tmp].first << "\t" << genusIdNamePairVec[tmp].second << endl;
		genus_ofs.close();
		ofstream family_ofs(family_file.c_str());
		family_ofs << "Family_index\tFamily_ID\tFamily_name" << endl;
		for(int tmp = 0; tmp < familyIdNamePairVec.size(); tmp++)
			family_ofs << tmp << "\t" << familyIdNamePairVec[tmp].first << "\t" << familyIdNamePairVec[tmp].second << endl;
		family_ofs.close();
		ofstream order_ofs(order_file.c_str());
		order_ofs << "Order_index\tOrder_ID\tOrder_name" << endl;
		for(int tmp = 0; tmp < orderIdNamePairVec.size(); tmp++)
			order_ofs << tmp << "\t" << orderIdNamePairVec[tmp].first << "\t" << orderIdNamePairVec[tmp].second << endl;		
		order_ofs.close();
		ofstream class_ofs(class_file.c_str());
		class_ofs << "Class_index\tClass_ID\tClass_name" << endl;
		for(int tmp = 0; tmp < classIdNamePairVec.size(); tmp++)
			class_ofs << tmp << "\t" << classIdNamePairVec[tmp].first << "\t" << classIdNamePairVec[tmp].second << endl;		
		class_ofs.close();
		ofstream phylum_ofs(phylum_file.c_str());
		phylum_ofs << "Phylum_index\tPhylum_ID\tPhylum_name" << endl;
		for(int tmp = 0; tmp < phylumIdNamePairVec.size(); tmp++)
			phylum_ofs << tmp << "\t" << phylumIdNamePairVec[tmp].first << "\t" << phylumIdNamePairVec[tmp].second << endl;
		phylum_ofs.close();
	}

	void reissueTaxoIdName_all()
	{
		this->reissueTaxoIdName_specific(speciesIdVec, speciesNameVec, reissuedSpeciesIdVec, speciesIdNamePairVec);
		this->reissueTaxoIdName_specific(genusIdVec, genusNameVec, reissuedGenusIdVec, genusIdNamePairVec);
		this->reissueTaxoIdName_specific(familyIdVec, familyNameVec, reissuedFamilyIdVec, familyIdNamePairVec);
		this->reissueTaxoIdName_specific(orderIdVec, orderNameVec, reissuedOrderIdVec, orderIdNamePairVec);
		this->reissueTaxoIdName_specific(classIdVec, classNameVec, reissuedClassIdVec, classIdNamePairVec);
		this->reissueTaxoIdName_specific(phylumIdVec, phylumNameVec, reissuedPhylumIdVec, phylumIdNamePairVec);											
		species_num = speciesIdNamePairVec.size();
		genus_num = genusIdNamePairVec.size();
		family_num = familyIdNamePairVec.size();
		order_num = orderIdNamePairVec.size();
		class_num = classIdNamePairVec.size();
		phylum_num = phylumIdNamePairVec.size();
	}

	void reissueTaxoIdName_specific(vector<int>& rawTaxoIdVec, vector<string>& rawTaxoNameVec,
		vector<int>& reissuedTaxoIdVec, vector< pair<int, string> >& reissuedTaxoIdNamePairVec)
	{
		for(int tmp = 0; tmp < rawTaxoIdVec.size(); tmp++)
		{
			int tmpRawTaxoId = rawTaxoIdVec[tmp];
			string tmpRawTaxoName = rawTaxoNameVec[tmp];
			int currentTaxoIdNamePairVecSize = reissuedTaxoIdNamePairVec.size();
			bool tmpAlreadyExist_bool = false;
			for(int tmp2 = 0; tmp2 < currentTaxoIdNamePairVecSize; tmp2++)
			{
				if(reissuedTaxoIdNamePairVec[tmp2].first == tmpRawTaxoId)
				{
					reissuedTaxoIdVec.push_back(tmp2);
					tmpAlreadyExist_bool = true;
					break;
				}
			}
			if(!tmpAlreadyExist_bool)
			{
				reissuedTaxoIdVec.push_back(currentTaxoIdNamePairVecSize);
				reissuedTaxoIdNamePairVec.push_back(pair<int,string>(tmpRawTaxoId, tmpRawTaxoName));
			}
		}
	}

	void initiate_bacterialTaxoFile_NCBIfullTaxoId2NameFile(
		string& bacterialTaxoFile, string& NCBIfullTaxoId2NameFile)
	{
		cout << "start to initiate_bacterialTaxoFile_NCBIfullTaxoId2NameFile" << endl;
		ifstream bacterialTaxo_ifs(bacterialTaxoFile.c_str());
		int tmpReissuedId = 0;
		while(!bacterialTaxo_ifs.eof())
		{
			string tmpStr;
			getline(bacterialTaxo_ifs, tmpStr);
			if(tmpStr == "")
				break;
			vector<string> tmpTaxoIdStrVec;
			parseStr2fieldVec(tmpTaxoIdStrVec, tmpStr);
			int tmpSpeciesId = atoi(tmpTaxoIdStrVec[0].c_str());
			int tmpGenusId = atoi(tmpTaxoIdStrVec[1].c_str());
			int tmpFamilyId = atoi(tmpTaxoIdStrVec[2].c_str());
			int tmpOrderId = atoi(tmpTaxoIdStrVec[3].c_str());
			int tmpClassId = atoi(tmpTaxoIdStrVec[4].c_str());
			int tmpPhylumId = atoi(tmpTaxoIdStrVec[5].c_str());
			speciesIdVec.push_back(tmpSpeciesId);
			genusIdVec.push_back(tmpGenusId);
			familyIdVec.push_back(tmpFamilyId);	
			orderIdVec.push_back(tmpOrderId);
			classIdVec.push_back(tmpClassId);
			phylumIdVec.push_back(tmpPhylumId);							
			tmpReissuedId ++;
		}
		bacterialTaxo_ifs.close();
		cout << "start to initiate NCBIfullTaxoId2Name" << endl;
		this->initiateNameVec_idVec_NCBIfullTaxoId2NameFile(NCBIfullTaxoId2NameFile);
	}

	void initiateNameVec_idVec_NCBIfullTaxoId2NameFile(string& NCBIfullTaxoId2NameFile)
	{
		NCBIfullTaxoID2Name_Info fullTaxoId2NameInfo;
		cout << "start to do initiateNameVec_idVec_NCBIfullTaxoId2NameFile" << endl;
		fullTaxoId2NameInfo.initiate_taxoID2NameFile(NCBIfullTaxoId2NameFile);
		int NCBIfullTaxoIdMax = fullTaxoId2NameInfo.return_NCBIfullTaxoIdMax();
		int speciesIdVecSize = speciesIdVec.size();
		cout << "speciesIdVecSize: " << speciesIdVecSize << endl;
		for(int tmp = 0; tmp < speciesIdVecSize; tmp++)
		{
			//cout << "tmp: " << tmp << endl;
			int tmpSpeciesId = speciesIdVec[tmp];
			int tmpGenusId = genusIdVec[tmp];
			int tmpFamilyId = familyIdVec[tmp];
			int tmpOrderId = orderIdVec[tmp];
			int tmpClassId = classIdVec[tmp];
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

			if((tmpFamilyId < 0)||(tmpFamilyId > NCBIfullTaxoIdMax))
				familyNameVec.push_back("NULL");
			else
				familyNameVec.push_back(fullTaxoId2NameInfo.return_taxoName(tmpFamilyId));

			if((tmpOrderId < 0)||(tmpOrderId > NCBIfullTaxoIdMax))
				orderNameVec.push_back("NULL");
			else
				orderNameVec.push_back(fullTaxoId2NameInfo.return_taxoName(tmpOrderId));

			if((tmpClassId < 0)||(tmpClassId > NCBIfullTaxoIdMax))
				classNameVec.push_back("NULL");
			else
				classNameVec.push_back(fullTaxoId2NameInfo.return_taxoName(tmpClassId));											

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