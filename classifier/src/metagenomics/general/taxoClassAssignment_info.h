// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#ifndef TAXOCLASSASSIGNMENT_INFO_H
#define TAXOCLASSASSIGNMENT_IFNO_H
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

class TaxoClassAssignment_Info
{
private:
	int taxo_class_num;
	int rank;
	string taxo_rank_name;

	vector<int> rawKmerClassIndex2taxoClassIndexVec;

	int genome_num;

	//////////////////// assignment on multiTaxoRank
	int total_taxo_set_num;
	int total_group_num;
	vector< pair<int,int> > groupId_range_eachRank_vec; // size = 6, groupId_range_eachRank_vec[8-taxoRank]
	vector< pair<int,int> > groupId_to_taxoReissuedIdRankPair_vec; // size = total_taxo_num + 2, <taxoReissuedId, taxoRank>
	//vector< vector<int> > groupIdPair_inLine_vecVec; // 
	char* groupIdPair_inLineOrNot_array; 
public:
	TaxoClassAssignment_Info()
	{
		//invalid_taxo_class_index = 0;
	}

	void get_groupIdRange_forSpecificRank(int& groupId_min, int& groupId_max, int specificRank)
	{
		groupId_min = groupId_range_eachRank_vec[8 - specificRank].first;
		groupId_max = groupId_range_eachRank_vec[8 - specificRank].second;
	}

	int return_groupId_forSpecificRank_min(int specificRank)
	{
		return groupId_range_eachRank_vec[8 - specificRank].first;
	}

	int return_groupId_forSpecificRank_max(int specificRank)
	{
		return groupId_range_eachRank_vec[8 - specificRank].second;
	}	

	void get_assignment_taxoId_allTaxoRank(
		#ifdef ASSIGN_INFO
		ofstream& taxo_assignment_ofs_assignInfo_detail,
		#endif
		// INPUT parameters
		BacterialTaxo_Info& bacterialTaxoInfo, vector<MaximalSpecificKmerWindow_Info>& tmpMSKWinfoVec, 
		//int assignment_windowSize_min,
		double assignment_windowSizeSquareSum_min,
		vector<int>& tmpInvalidKmerCountVec, vector<int>& tmpRepetitiveKmerCountVec, vector<int>& tmpDiscriminativeKmerCountVec,
		// OUTPUT parameters
		int& tmp_assignment_taxoId_species, int& tmp_assignment_taxoId_genus, int& tmp_assignment_taxoId_family,
		int& tmp_assignment_taxoId_order, int& tmp_assignment_taxoId_class, int& tmp_assignment_taxoId_phylum)
	{
		#ifdef ASSIGN_INFO
		taxo_assignment_ofs_assignInfo_detail << "WindowSizeSquareSum-based approach results: " << endl;
		#endif		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////// WINSOW-SIZE based  //////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// vector<int> tmp_assignment_reissuedIdVec_specificRank;
		// for(int tmpRank = 8; tmpRank >= 3; tmpRank--)
		// {
		// 	int tmpRankIndex = 8 - tmpRank;

		// 	int tmp_reissuedId_best, tmp_windowSize_best, tmp_reissuedId_secondBest, tmp_windowSize_secondBest;
		// 	tmpMSKWinfoVec[tmpRankIndex].get_taxoReissuedId_maximalWindow_best_secondBest(tmp_reissuedId_best, 
		// 		tmp_windowSize_best, tmp_reissuedId_secondBest, tmp_windowSize_secondBest, assignment_windowSize_min);

		// 	#ifdef ASSIGN_INFO
		// 	taxo_assignment_ofs_assignInfo_detail << "Rank: " << tmpRank << endl;
		// 	taxo_assignment_ofs_assignInfo_detail 
		// 		<< "Best: " << bacterialTaxoInfo.return_taxo_Id(tmpRank, tmp_reissuedId_best) 
		// 		<< "-" << tmp_windowSize_best 
		// 		<< " SecondBest: " << bacterialTaxoInfo.return_taxo_Id(tmpRank, tmp_reissuedId_secondBest)
		// 		<< "-" << tmp_windowSize_secondBest << endl;
		// 	#endif

		// 	if(tmp_windowSize_best == tmp_windowSize_secondBest)
		// 		tmp_assignment_reissuedIdVec_specificRank.push_back(-2);
		// 	else
		// 		tmp_assignment_reissuedIdVec_specificRank.push_back(tmp_reissuedId_best);
		// }
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////   WINDOW-SIZE-SQUARESUM based  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		vector<int> tmp_assignment_reissuedIdVec_specificRank;
		for(int tmpRank = 8; tmpRank >= 3; tmpRank--)
		{
			int tmpRankIndex = 8 - tmpRank;
			int tmp_reissuedId_best, tmp_windowSizeSquareSum_best, tmp_reissuedId_secondBest, tmp_windowSizeSquareSum_secondBest;
			tmpMSKWinfoVec[tmpRankIndex].get_taxoReissuedId_maximalWindowSquareSum_best_secondBest(
				tmp_reissuedId_best, tmp_windowSizeSquareSum_best, 
				tmp_reissuedId_secondBest, tmp_windowSizeSquareSum_secondBest, assignment_windowSizeSquareSum_min);
			
			#ifdef ASSIGN_INFO
			taxo_assignment_ofs_assignInfo_detail << "Rank: " << tmpRank;
			//taxo_assignment_ofs_assignInfo_detail << "Map:" << endl << tmpMSKWinfoVec[tmpRankIndex].print_taxoId2windowMap_2_string() << endl;
			taxo_assignment_ofs_assignInfo_detail
				<< " Best: " << bacterialTaxoInfo.return_taxo_Id(tmpRank, tmp_reissuedId_best) 
				<< "-" << tmp_windowSizeSquareSum_best 
				<< " SecondBest: " << bacterialTaxoInfo.return_taxo_Id(tmpRank, tmp_reissuedId_secondBest)
				<< "-" << tmp_windowSizeSquareSum_secondBest << endl;
			#endif

			if(tmp_windowSizeSquareSum_best == tmp_windowSizeSquareSum_secondBest)
				tmp_assignment_reissuedIdVec_specificRank.push_back(-2);
			else
				tmp_assignment_reissuedIdVec_specificRank.push_back(tmp_reissuedId_best);
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////    check consistency between taxoId assignment on different taxo levels  ////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// this->update_assignment_reissuedId_checkAssignmentConsistency_eachTaxoRankPair(
		// 	tmp_assignment_reissuedIdVec_specificRank[0], tmp_assignment_reissuedIdVec_specificRank[1],
		// 	tmp_assignment_reissuedIdVec_specificRank[2], tmp_assignment_reissuedIdVec_specificRank[3],
		// 	tmp_assignment_reissuedIdVec_specificRank[4], tmp_assignment_reissuedIdVec_specificRank[5]);

		this->update_assignment_reissuedId_checkAssignmentConsistency_topDown(
			tmp_assignment_reissuedIdVec_specificRank[0], tmp_assignment_reissuedIdVec_specificRank[1],
			tmp_assignment_reissuedIdVec_specificRank[2], tmp_assignment_reissuedIdVec_specificRank[3],
			tmp_assignment_reissuedIdVec_specificRank[4], tmp_assignment_reissuedIdVec_specificRank[5]);
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////// for either windowSize-based or windowSizeSquareSum-based approach /////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		tmp_assignment_taxoId_species = bacterialTaxoInfo.return_taxo_Id(8, tmp_assignment_reissuedIdVec_specificRank[0]);
		tmp_assignment_taxoId_genus = bacterialTaxoInfo.return_taxo_Id(7, tmp_assignment_reissuedIdVec_specificRank[1]);
		tmp_assignment_taxoId_family = bacterialTaxoInfo.return_taxo_Id(6, tmp_assignment_reissuedIdVec_specificRank[2]);
		tmp_assignment_taxoId_order = bacterialTaxoInfo.return_taxo_Id(5, tmp_assignment_reissuedIdVec_specificRank[3]);
		tmp_assignment_taxoId_class = bacterialTaxoInfo.return_taxo_Id(4, tmp_assignment_reissuedIdVec_specificRank[4]);
		tmp_assignment_taxoId_phylum = bacterialTaxoInfo.return_taxo_Id(3, tmp_assignment_reissuedIdVec_specificRank[5]);
	}

	void update_assignment_reissuedId_checkAssignmentConsistency_topDown(
		int& assignment_reissuedId_species, int& assignment_reissuedId_genus, int& assignment_reissuedId_family, 
		int& assignment_reissuedId_order, int& assignment_reissuedId_class, int& assignment_reissuedId_phylum)
	{
		if(assignment_reissuedId_species >= 0)
		{
			if(((assignment_reissuedId_genus >= 0)
					&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_species, 8, assignment_reissuedId_genus, 7))))
				||((assignment_reissuedId_family >= 0)
					&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_species, 8, assignment_reissuedId_family, 6))))
				||((assignment_reissuedId_order >= 0)
					&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_species, 8, assignment_reissuedId_order, 5))))
				||((assignment_reissuedId_class >= 0)
					&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_species, 8, assignment_reissuedId_class, 4))))
				||((assignment_reissuedId_phylum >= 0)
					&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_species, 8, assignment_reissuedId_phylum, 3)))))
				assignment_reissuedId_species = -1;
		}
		if(assignment_reissuedId_genus >= 0)
		{
			if(((assignment_reissuedId_family >= 0)
					&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_genus, 7, assignment_reissuedId_family, 6))))
				||((assignment_reissuedId_order >= 0)
					&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_genus, 7, assignment_reissuedId_order, 5))))
				||((assignment_reissuedId_class >= 0)
					&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_genus, 7, assignment_reissuedId_class, 4))))
				||((assignment_reissuedId_phylum >= 0)
					&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_genus, 7, assignment_reissuedId_phylum, 3)))))					
				assignment_reissuedId_genus = -1;
		}
		if(assignment_reissuedId_family >= 0)
		{
			if(((assignment_reissuedId_order >= 0)
					&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_family, 6, assignment_reissuedId_order, 5))))
				||((assignment_reissuedId_class >= 0)
					&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_family, 6, assignment_reissuedId_class, 4))))
				||((assignment_reissuedId_phylum >= 0)
					&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_family, 6, assignment_reissuedId_phylum, 3)))))
				assignment_reissuedId_family = -1;
		}
		if(assignment_reissuedId_order >= 0)
		{
			if(((assignment_reissuedId_class >= 0)
					&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_order, 5, assignment_reissuedId_class, 4))))
				||((assignment_reissuedId_phylum >= 0)
					&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_order, 5, assignment_reissuedId_phylum, 3)))))		
				assignment_reissuedId_order = -1;
		}
		if(assignment_reissuedId_class >= 0)
		{
			if((assignment_reissuedId_phylum >= 0)
					&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_class, 4, assignment_reissuedId_phylum, 3))))
				assignment_reissuedId_class = -1;
		}			
	}

	void update_assignment_reissuedId_checkAssignmentConsistency_eachTaxoRankPair(
		int& assignment_reissuedId_species, int& assignment_reissuedId_genus, int& assignment_reissuedId_family, 
		int& assignment_reissuedId_order, int& assignment_reissuedId_class, int& assignment_reissuedId_phylum)
	{
		bool consistent_or_not_eachTaxoRankPair_bool = this->checkAssignmentConsistency_eachTaxoRankPair(
			assignment_reissuedId_species, assignment_reissuedId_genus, assignment_reissuedId_family, 
			assignment_reissuedId_order, assignment_reissuedId_class, assignment_reissuedId_phylum);
		if(!consistent_or_not_eachTaxoRankPair_bool)
		{
			assignment_reissuedId_species = -1;
			assignment_reissuedId_genus = -1;
			assignment_reissuedId_family = -1;
			assignment_reissuedId_order = -1;
			assignment_reissuedId_class = -1;
			assignment_reissuedId_phylum = -1;
		}
	}

	bool checkAssignmentConsistency_eachTaxoRankPair(int assignment_reissuedId_species, int assignment_reissuedId_genus,
		int assignment_reissuedId_family, int assignment_reissuedId_order, int assignment_reissuedId_class, int assignment_reissuedId_phylum)
	{
		if((assignment_reissuedId_species >= 0)&&(assignment_reissuedId_genus >= 0)
			&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_species, 8, assignment_reissuedId_genus, 7))))
			return false;
		if((assignment_reissuedId_species >= 0)&&(assignment_reissuedId_family >= 0)
			&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_species, 8, assignment_reissuedId_family, 6))))
			return false;
		if((assignment_reissuedId_species >= 0)&&(assignment_reissuedId_order >= 0)
			&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_species, 8, assignment_reissuedId_order, 5))))
			return false;
		if((assignment_reissuedId_species >= 0)&&(assignment_reissuedId_class >= 0)
			&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_species, 8, assignment_reissuedId_class, 4))))
			return false;
		if((assignment_reissuedId_species >= 0)&&(assignment_reissuedId_phylum >= 0)
			&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_species, 8, assignment_reissuedId_phylum, 3))))
			return false;


		if((assignment_reissuedId_genus >= 0)&&(assignment_reissuedId_family >= 0)
			&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_genus, 7, assignment_reissuedId_family, 6))))
			return false;
		if((assignment_reissuedId_genus >= 0)&&(assignment_reissuedId_order >= 0)
			&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_genus, 7, assignment_reissuedId_order, 5))))
			return false;
		if((assignment_reissuedId_genus >= 0)&&(assignment_reissuedId_class >= 0)
			&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_genus, 7, assignment_reissuedId_class, 4))))
			return false;
		if((assignment_reissuedId_genus >= 0)&&(assignment_reissuedId_phylum >= 0)
			&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_genus, 7, assignment_reissuedId_phylum, 3))))
			return false;					


		if((assignment_reissuedId_family >= 0)&&(assignment_reissuedId_order >= 0)
			&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_family, 6, assignment_reissuedId_order, 5))))
			return false;
		if((assignment_reissuedId_family >= 0)&&(assignment_reissuedId_class >= 0)
			&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_family, 6, assignment_reissuedId_class, 4))))
			return false;
		if((assignment_reissuedId_family >= 0)&&(assignment_reissuedId_phylum >= 0)
			&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_family, 6, assignment_reissuedId_phylum, 3))))
			return false;		

		if((assignment_reissuedId_order >= 0)&&(assignment_reissuedId_class >= 0)
			&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_order, 5, assignment_reissuedId_class, 4))))
			return false;
		if((assignment_reissuedId_order >= 0)&&(assignment_reissuedId_phylum >= 0)
			&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_order, 5, assignment_reissuedId_phylum, 3))))
			return false;


		if((assignment_reissuedId_class >= 0)&&(assignment_reissuedId_phylum >= 0)
			&&(!(this->taxoIdRankPairPair_inLineOrNot(assignment_reissuedId_class, 4, assignment_reissuedId_phylum, 3))))
			return false;
		return true;										
	}

	void filter_KmerValueVec(vector<valueT>& tmpKmerValueVec_filtered, vector<valueT>& tmpKmerValueVec)
	{
		int tmpKmerValueVecSize = tmpKmerValueVec.size();
		if(tmpKmerValueVecSize == 0)
			return;
		else if(tmpKmerValueVecSize == 1)
		{	
			tmpKmerValueVec_filtered.push_back(tmpKmerValueVec[0]);
			return;
		}
		else if(tmpKmerValueVecSize == 2)
		{
			#ifndef NO_FILTERING_AT_BOTH_ENDS
			int tmp1stValue = tmpKmerValueVec[0];
			int tmp2ndValue = tmpKmerValueVec[1];
			if(this->groupIdPair_inLineOrNot(tmp1stValue, tmp2ndValue))
			{
				tmpKmerValueVec_filtered.push_back(tmp1stValue);
				tmpKmerValueVec_filtered.push_back(tmp2ndValue);
			}
			else
			{
				tmpKmerValueVec_filtered.push_back(0);
				tmpKmerValueVec_filtered.push_back(0);
			}
			#else
			tmpKmerValueVec_filtered.push_back(tmpKmerValueVec[0]);
			tmpKmerValueVec_filtered.push_back(tmpKmerValueVec[1]);
			#endif
			return;
		}
		// the 1st KmerValue
		#ifndef NO_FILTERING_AT_BOTH_ENDS
		int tmp1stValue = tmpKmerValueVec[0];
		int tmp2ndValue = tmpKmerValueVec[1];
		if(this->groupIdPair_inLineOrNot(tmp1stValue, tmp2ndValue))
			tmpKmerValueVec_filtered.push_back(tmp1stValue);
		else
			tmpKmerValueVec_filtered.push_back(0);
		#else
		tmpKmerValueVec_filtered.push_back(tmpKmerValueVec[0]);
		#endif

		// KmerValue except the 1st and the last one
		for(int tmp = 1; tmp <= tmpKmerValueVecSize - 2; tmp++)
		{
			int tmpThisValue = tmpKmerValueVec[tmp];
			int tmpLastValue = tmpKmerValueVec[tmp-1];
			int tmpNextValue = tmpKmerValueVec[tmp+1];
			if((this->groupIdPair_inLineOrNot(tmpLastValue, tmpThisValue))
				||(this->groupIdPair_inLineOrNot(tmpThisValue, tmpNextValue)))
				tmpKmerValueVec_filtered.push_back(tmpThisValue);
			else
				tmpKmerValueVec_filtered.push_back(0);
		}

		// the last KmerValue
		#ifndef NO_FILTERING_AT_BOTH_ENDS
		int penultimateKmerValue = tmpKmerValueVec[tmpKmerValueVecSize - 2];
		int lastKmerValue = tmpKmerValueVec[tmpKmerValueVecSize - 1];
		if(this->groupIdPair_inLineOrNot(penultimateKmerValue, lastKmerValue))
			tmpKmerValueVec_filtered.push_back(tmpKmerValueVec[tmpKmerValueVecSize - 1]);
		else
			tmpKmerValueVec_filtered.push_back(0);
		#else
		tmpKmerValueVec_filtered.push_back(tmpKmerValueVec[tmpKmerValueVecSize - 1]);
		#endif
	}

	void get_taxoId2windowVecMapVec_specific_allTaxoRank(
		vector<MaximalSpecificKmerWindow_Info>& maximalSpecificKmerWindowInfoVec,
		vector<valueT>& tmpKmerValueVec, int Kmer_length, BacterialTaxo_Info& bacterialTaxoInfo,
		vector<int>& invalidKmerCountVec, vector<int>& repetitiveKmerCountVec, vector<int>& discriminativeKmerCountVec)
	{
		//cout << "start to do get_taxoId2windowVecMapVec_specific_allTaxoRank(" << endl;
		int taxoRank_highest = 3;
		int taxoRank_lowest = 8;
		//cout << "total_taxo_set_num + 1: " << total_taxo_set_num + 1 << endl;
		for(int specificRank = taxoRank_lowest; specificRank >= taxoRank_highest; specificRank--)
		{
			//cout << "tmpSpecificRank: " << specificRank << endl;
			MaximalSpecificKmerWindow_Info tmpMSKWinfo;
			int tmpInvalidKmerCount = 0, tmpRepetitiveKmerCount = 0, tmpDiscriminativeKmerCount = 0;
			if(tmpKmerValueVec.size() == 0)
			{
				maximalSpecificKmerWindowInfoVec.push_back(tmpMSKWinfo);
				invalidKmerCountVec.push_back(0);
				repetitiveKmerCountVec.push_back(0);
				discriminativeKmerCountVec.push_back(0);
				continue;
			}
			//map<int,int> groupId2countMap;
			int currentWindow_taxoReissuedId = -1;
			int currentWindow_1stKmerIndexInVec = -1;
			int currentWindow_KmerCount = -1;	
			int specificRank_groupId_max = this->return_groupId_forSpecificRank_max(specificRank);
			// all the kmers except the last one
			if(tmpKmerValueVec.size() > 1)
			{	
				for(int tmp = 0; tmp < tmpKmerValueVec.size() - 1; tmp ++)
				{
					valueT tmpValue = tmpKmerValueVec[tmp];
					//cout << "tmpValue: " << tmpValue << endl;
					if((tmpValue < 1)||(tmpValue > total_taxo_set_num + 1)) // invalid, insert 2 MaximalSpecificKmerWindow_Info
					{
						//cout << "invalid!" << endl;
						tmpInvalidKmerCount ++;
						//insert 2 MaximalSpecificKmerWindow_Info
						#ifdef W2
						if((currentWindow_taxoReissuedId >= 0)&&(currentWindow_KmerCount > 1))	
							tmpMSKWinfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
								currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
						#else
						if((currentWindow_taxoReissuedId >= 0))//&&(currentWindow_KmerCount > 1))	
							tmpMSKWinfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
								currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
						#endif
						currentWindow_taxoReissuedId = -1;
						currentWindow_1stKmerIndexInVec = -1;
						currentWindow_KmerCount = -1;
					}
					else if(tmpValue > specificRank_groupId_max) // specific-taxo-level-repetitive,
					{
						//cout << "repetitive! " << endl;
						tmpRepetitiveKmerCount ++;
						//insert 2 MaximalSpecificKmerWindow_Info
						#ifdef W2
						if((currentWindow_taxoReissuedId >= 0)&&(currentWindow_KmerCount > 1))
							tmpMSKWinfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
								currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
						#else
						if((currentWindow_taxoReissuedId >= 0))//&&(currentWindow_KmerCount > 1))
							tmpMSKWinfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
								currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
						#endif
						currentWindow_taxoReissuedId = -1;
						currentWindow_1stKmerIndexInVec = -1;
						currentWindow_KmerCount = -1;	
					}
					else // specific-taxo-level-discriminative
					{
						//cout << "discriminative!" << endl;
						tmpDiscriminativeKmerCount ++;
						int tmpGroupId = tmpValue;
						int tmpTaxoRank, tmpTaxoReissuedId;
						this->get_taxoReissuedIdRank_from_groupId(tmpGroupId, tmpTaxoReissuedId, tmpTaxoRank);
						//cout << "tmpTaxoRank: " << tmpTaxoRank << endl;
						//cout << "tmpTaxoReissuedId: " << tmpTaxoReissuedId << endl;
						int tmpTaxoReissuedId_thisRank;
						if(tmpTaxoRank > specificRank) // lower rank
							tmpTaxoReissuedId_thisRank 
								= bacterialTaxoInfo.get_higherRankTaxoReissuedId_from_lowerRankTaxoReissuedId(
									specificRank, tmpTaxoReissuedId, tmpTaxoRank);
						else // this rank
							tmpTaxoReissuedId_thisRank = tmpTaxoReissuedId;
						//cout << "tmpTaxoReissuedId_thisRank: " << tmpTaxoReissuedId_thisRank << endl;
						if(tmpTaxoReissuedId_thisRank != currentWindow_taxoReissuedId)
						{
							//cout << " insert 2 MaximalSpecificKmerWindow_Info, restart a new window" << endl;
							#ifdef W2
							if((currentWindow_taxoReissuedId >= 0)&&(currentWindow_KmerCount > 1))
								tmpMSKWinfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
									currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
							#else
							if((currentWindow_taxoReissuedId >= 0))//&&(currentWindow_KmerCount > 1))
								tmpMSKWinfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
									currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
							#endif
							currentWindow_taxoReissuedId = tmpTaxoReissuedId_thisRank;
							currentWindow_1stKmerIndexInVec = tmp;
							currentWindow_KmerCount = 1;						
						}
						else //tmpTaxoReissuedId_thisRank == currentWindow_taxoReissuedId, update current window
							currentWindow_KmerCount ++;
					}
					//cout << "tmpMap:" << endl << tmpMSKWinfo.print_taxoId2windowMap_2_string() << endl;
				}
			}
			// last Kmer
			valueT lastValue = tmpKmerValueVec[tmpKmerValueVec.size() - 1];
			//cout << "lastValue: " << lastValue << endl;
			//cout << "currentWindow_taxoReissuedId: " << currentWindow_taxoReissuedId << endl;
			if((lastValue < 1)||(lastValue > specificRank_groupId_max))// invalid or repetitive Kmer
			{
				//cout << "invalid or repetitive Kmer" << endl;
				if((lastValue < 1)||(lastValue > total_taxo_set_num + 1))
					tmpInvalidKmerCount ++;
				else
					tmpRepetitiveKmerCount ++;
				//insert 2 MaximalSpecificKmerWindow_Info
				#ifdef W2
				if((currentWindow_taxoReissuedId >= 0)&&(currentWindow_KmerCount > 1))
					tmpMSKWinfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
						currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
				#else
				if((currentWindow_taxoReissuedId >= 0))//&&(currentWindow_KmerCount > 1))
					tmpMSKWinfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
						currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
				#endif
			}
			else // discriminative Kmer
			{	
				//cout << "discriminative Kmer" << endl;
				tmpDiscriminativeKmerCount ++;
				int lastGroupId = lastValue;
				int lastTaxoRank, lastTaxoReissuedId;
				this->get_taxoReissuedIdRank_from_groupId(lastGroupId, lastTaxoReissuedId, lastTaxoRank);			
				
				int lastTaxoReissuedId_thisRank;
				if(lastTaxoRank > specificRank) // lower rank
					lastTaxoReissuedId_thisRank 
						= bacterialTaxoInfo.get_higherRankTaxoReissuedId_from_lowerRankTaxoReissuedId(
							specificRank, lastTaxoReissuedId, lastTaxoRank);
				else // this rank
					lastTaxoReissuedId_thisRank = lastTaxoReissuedId;
				if(lastTaxoReissuedId_thisRank != currentWindow_taxoReissuedId)
				{
					// insert 2 MaximalSpecificKmerWindow_Info, restart a new window
					#ifdef W2
					if((currentWindow_taxoReissuedId >= 0)&&(currentWindow_KmerCount > 1))
						tmpMSKWinfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
							currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
					#else
					if((currentWindow_taxoReissuedId >= 0))//&&(currentWindow_KmerCount > 1))
						tmpMSKWinfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
							currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);	
					#endif
					currentWindow_taxoReissuedId = lastTaxoReissuedId_thisRank;
					currentWindow_1stKmerIndexInVec = tmpKmerValueVec.size() - 1;
					currentWindow_KmerCount = 1;	
					#ifdef W2
					if((currentWindow_taxoReissuedId >= 0)&&(currentWindow_KmerCount > 1))
						tmpMSKWinfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
							currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
					#else					
					if((currentWindow_taxoReissuedId >= 0))//&&(currentWindow_KmerCount > 1)) // actually do not need checking, as lastTaxoReissuedId_thisRank >= 0
						tmpMSKWinfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
							currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
					#endif
				}
				else //lastTaxoReissuedId_thisRank == currentWindow_taxoReissuedId, update current window and insert
				{
					currentWindow_KmerCount ++;

					#ifdef W2
					if((currentWindow_taxoReissuedId >= 0)&&(currentWindow_KmerCount > 1))
						tmpMSKWinfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
							currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
					#else
					if((currentWindow_taxoReissuedId >= 0))//&&(currentWindow_KmerCount > 1))
						tmpMSKWinfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
							currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
					#endif
				}
			}	
			maximalSpecificKmerWindowInfoVec.push_back(tmpMSKWinfo);
			invalidKmerCountVec.push_back(tmpInvalidKmerCount);
			repetitiveKmerCountVec.push_back(tmpRepetitiveKmerCount);
			discriminativeKmerCountVec.push_back(tmpDiscriminativeKmerCount);
			//cout << "tmpMap:" << endl << tmpMSKWinfo.print_taxoId2windowMap_2_string() << endl;
		}
	}

	/*
	void get_taxoId2windowVecMap_specific(MaximalSpecificKmerWindow_Info& maximalSpecificKmerWindowInfo, 
		vector<valueT>& tmpKmerValueVec, 
		int Kmer_length, int specificRank, BacterialTaxo_Info& bacterialTaxoInfo,
		int& invalidKmerCount, int& repetitiveKmerCount, int& discriminativeKmerCount)
	{
		int tmpInvalidKmerCount = 0, tmpRepetitiveKmerCount = 0, tmpDiscriminativeKmerCount = 0;
		//map<int,int> groupId2countMap;
		int currentWindow_taxoReissuedId = -1;
		int currentWindow_1stKmerIndexInVec = -1;
		int currentWindow_KmerCount = -1;

		if(tmpKmerValueVec.size() == 0)
			return;

		// all the kmers except the last one
		for(int tmp = 0; tmp < tmpKmerValueVec.size() - 1; tmp ++)
		{
			valueT tmpValue = tmpKmerValueVec[tmp];
			if((tmpValue < 1)||(tmpValue > total_taxo_set_num + 1)) // invalid, insert 2 MaximalSpecificKmerWindow_Info
			{
				tmpInvalidKmerCount ++;
				//insert 2 MaximalSpecificKmerWindow_Info
				if(currentWindow_taxoReissuedId >= 0)
					maximalSpecificKmerWindowInfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
						currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
				currentWindow_taxoReissuedId = -1;
				currentWindow_1stKmerIndexInVec = -1;
				currentWindow_KmerCount = -1;
			}
			else if(tmpValue == total_taxo_set_num + 1) // phylum-repetitive,
			{
				tmpRepetitiveKmerCount ++;
				//insert 2 MaximalSpecificKmerWindow_Info
				if(currentWindow_taxoReissuedId >= 0)
					maximalSpecificKmerWindowInfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
						currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
				currentWindow_taxoReissuedId = -1;
				currentWindow_1stKmerIndexInVec = -1;
				currentWindow_KmerCount = -1;	
			}
			else // branch-specific
			{
				tmpDiscriminativeKmerCount ++;
				int tmpGroupId = tmpValue;
				int tmpTaxoRank, tmpTaxoReissuedId;
				this->get_taxoReissuedIdRank_from_groupId(tmpGroupId, tmpTaxoReissuedId, tmpTaxoRank);				
				if(tmpTaxoRank < specificRank)//higher rank, doesn't count for get_taxoId2windowVecMap_specific
				{
					//insert 2 MaximalSpecificKmerWindow_Info
					if(currentWindow_taxoReissuedId >= 0)
						maximalSpecificKmerWindowInfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
							currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
					currentWindow_taxoReissuedId = -1;
					currentWindow_1stKmerIndexInVec = -1;
					currentWindow_KmerCount = -1;
				}
				else// lower rank or this rank
				{
					int tmpTaxoReissuedId_thisRank;
					if(tmpTaxoRank > specificRank) // lower rank
						tmpTaxoReissuedId_thisRank 
							= bacterialTaxoInfo.get_higherRankTaxoReissuedId_from_lowerRankTaxoReissuedId(
								specificRank, tmpTaxoReissuedId, tmpTaxoRank);
					else // this rank
						tmpTaxoReissuedId_thisRank = tmpTaxoReissuedId;
					if(tmpTaxoReissuedId_thisRank != currentWindow_taxoReissuedId)
					{
						// insert 2 MaximalSpecificKmerWindow_Info, restart a new window
						if(currentWindow_taxoReissuedId >= 0)
							maximalSpecificKmerWindowInfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
								currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
						currentWindow_taxoReissuedId = tmpTaxoReissuedId_thisRank;
						currentWindow_1stKmerIndexInVec = tmp;
						currentWindow_KmerCount = 1;						
					}
					else //tmpTaxoReissuedId_thisRank == currentWindow_taxoReissuedId, update current window
						currentWindow_KmerCount ++;
				}
			}
		}

		// the last Kmer
		valueT lastValue = tmpKmerValueVec[tmpKmerValueVec.size() - 1];
		if((lastValue < 1)||(lastValue >= total_taxo_set_num + 1))// invalid
		{
			if(lastValue == total_taxo_set_num + 1)
				tmpRepetitiveKmerCount ++;
			else
				tmpInvalidKmerCount ++;			
			//insert 2 MaximalSpecificKmerWindow_Info
			if(currentWindow_taxoReissuedId >= 0)
				maximalSpecificKmerWindowInfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
					currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
		}
		else
		{
			tmpDiscriminativeKmerCount ++;
			int lastGroupId = lastValue;
			int lastTaxoRank, lastTaxoReissuedId;
			this->get_taxoReissuedIdRank_from_groupId(lastGroupId, lastTaxoReissuedId, lastTaxoRank);				
			if(lastTaxoRank < specificRank)//higher rank, doesn't count for get_taxoId2windowVecMap_specific
			{
				//insert 2 MaximalSpecificKmerWindow_Info
				if(currentWindow_taxoReissuedId >= 0)
					maximalSpecificKmerWindowInfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
						currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
			}
			else// lower rank or this rank
			{
				int lastTaxoReissuedId_thisRank;
				if(lastTaxoRank > specificRank) // lower rank
					lastTaxoReissuedId_thisRank 
						= bacterialTaxoInfo.get_higherRankTaxoReissuedId_from_lowerRankTaxoReissuedId(
							specificRank, lastTaxoReissuedId, lastTaxoRank);
				else // this rank
					lastTaxoReissuedId_thisRank = lastTaxoReissuedId;
				if(lastTaxoReissuedId_thisRank != currentWindow_taxoReissuedId)
				{
					// insert 2 MaximalSpecificKmerWindow_Info, restart a new window
					if(currentWindow_taxoReissuedId >= 0)
						maximalSpecificKmerWindowInfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
							currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);	
					currentWindow_taxoReissuedId = lastTaxoReissuedId_thisRank;
					currentWindow_1stKmerIndexInVec = tmpKmerValueVec.size() - 1;
					currentWindow_KmerCount = 1;	
					maximalSpecificKmerWindowInfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
						currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
				}
				else //lastTaxoReissuedId_thisRank == currentWindow_taxoReissuedId, update current window and insert
				{
					currentWindow_KmerCount ++;
					if(currentWindow_taxoReissuedId >= 0)
						maximalSpecificKmerWindowInfo.insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(
							currentWindow_taxoReissuedId, currentWindow_1stKmerIndexInVec, currentWindow_KmerCount);
				}
			}
		}

		invalidKmerCount = tmpInvalidKmerCount;
		repetitiveKmerCount = tmpRepetitiveKmerCount;
		discriminativeKmerCount = tmpDiscriminativeKmerCount;
	}*/

	void convertGroupIdVec2reissuedIdVec(vector<valueT>& tmpKmerValueVec, 
		vector<int>& specificRankTaxoReissuedIdVec, int specificRank, BacterialTaxo_Info& bacterialTaxoInfo)
	{
		for(int tmp = 0; tmp < tmpKmerValueVec.size(); tmp ++)
		{
			valueT tmpValue = tmpKmerValueVec[tmp];
			if((tmpValue < 1)||(tmpValue > total_taxo_set_num + 1)) // invalid
				specificRankTaxoReissuedIdVec.push_back(-3);
			else if(tmpValue == total_taxo_set_num + 1)
				specificRankTaxoReissuedIdVec.push_back(-2);
			else
			{
				int tmpGroupId = tmpValue;
				int tmpTaxoRank, tmpTaxoReissuedId;
				this->get_taxoReissuedIdRank_from_groupId(tmpGroupId, tmpTaxoReissuedId, tmpTaxoRank);					
				if(tmpTaxoRank < specificRank)//higher rank, doesn't count for get_taxoId2windowVecMap_specific
					specificRankTaxoReissuedIdVec.push_back(-2);
				else// lower rank or this rank
				{
					int tmpTaxoReissuedId_thisRank;
					if(tmpTaxoRank > specificRank) // lower rank
						tmpTaxoReissuedId_thisRank 
							= bacterialTaxoInfo.get_higherRankTaxoReissuedId_from_lowerRankTaxoReissuedId(
								specificRank, tmpTaxoReissuedId, tmpTaxoRank);
					else // this rank
						tmpTaxoReissuedId_thisRank = tmpTaxoReissuedId;
					specificRankTaxoReissuedIdVec.push_back(tmpTaxoReissuedId_thisRank);
				}
			}
		}		
	}

	void get_taxoReissuedIdRank_from_groupId(int tmpGroupId, int& tmpTaxoReissuedId, int& tmpTaxoRank)
	{
		//tmpTaxoReissuedId = -1;
		//tmpTaxoRank = -1;
		// for(int tmpRank = 8; tmpRank >= 3; tmpRank --)
		// {
		// 	int tmpIndex = 8 - tmpRank;
		// 	int tmpRank_groupId_min = groupId_range_eachRank_vec[tmpIndex].first;
		// 	int tmpRank_groupId_max = groupId_range_eachRank_vec[tmpIndex].second;
		// 	if((tmpRank_groupId_min <= tmpGroupId)&&(tmpGroupId <= tmpRank_groupId_max))
		// 	{
		// 		tmpTaxoRank = tmpRank;
		// 		tmpTaxoReissuedId = tmpGroupId - tmpRank_groupId_min;
		// 		return;
		// 	}
		// }
		// cout << "current tmpGroupId not valid for any taxoRank, tmpGroupId: " << tmpGroupId << endl;
		// exit(1);
		if((tmpGroupId < 1)||(tmpGroupId > total_taxo_set_num))
		{
			cout << "invalid tmpGroupId for get_taxoReissuedIdRank_from_groupId, tmpGroupId: " << tmpGroupId << endl;
			exit(1);
		}
		tmpTaxoReissuedId = groupId_to_taxoReissuedIdRankPair_vec[tmpGroupId].first;
		tmpTaxoRank = groupId_to_taxoReissuedIdRankPair_vec[tmpGroupId].second;
	}

	int get_groupId_from_reissuedTaxoId_taxoRank(int tmpReissuedTaxoId, int tmpTaxoRank)
	{
		switch(tmpTaxoRank)
		{
			case 3: // phylum
				return (groupId_range_eachRank_vec[5].first + tmpReissuedTaxoId);
			case 4: // class
				return (groupId_range_eachRank_vec[4].first + tmpReissuedTaxoId);
			case 5: // order
				return (groupId_range_eachRank_vec[3].first + tmpReissuedTaxoId);
			case 6: // family
				return (groupId_range_eachRank_vec[2].first + tmpReissuedTaxoId);
			case 7: // genus
				return (groupId_range_eachRank_vec[1].first + tmpReissuedTaxoId);
			case 8: // species
				return (groupId_range_eachRank_vec[0].first + tmpReissuedTaxoId);
			default:
				cout << "Invalid tmpTaxoRank: " << tmpTaxoRank << endl;
		}
	}

	void initiate_groupIdPair_inLineOrNot_array(BacterialTaxo_Info& bacterialTaxoInfo)
	{
		groupIdPair_inLineOrNot_array = new char[total_group_num * total_group_num];
		for(int tmp = 1; tmp <= total_taxo_set_num + 1; tmp++)
		{
			groupIdPair_inLineOrNot_array[tmp * total_group_num + tmp] = 1;
			groupIdPair_inLineOrNot_array[tmp * total_group_num + total_taxo_set_num + 1] = 1;
			groupIdPair_inLineOrNot_array[(total_taxo_set_num + 1) * total_group_num + tmp] = 1;
		}
		int species_num_total = bacterialTaxoInfo.return_species_num();
		for(int tmp = 0; tmp < species_num_total; tmp++)
		{
			//vector< pair<int,int> > tmpReissuedTaxoIdRankPairVec;
			vector<int> tmpGroupIdVec;
			int tmp_reissuedId_species, tmp_reissuedId_genus, tmp_reissuedId_family,
				tmp_reissuedId_order, tmp_reissuedId_class, tmp_reissuedId_phylum;
			bacterialTaxoInfo.get_reissuedId_allRank(tmp, tmp_reissuedId_species, tmp_reissuedId_genus,
				tmp_reissuedId_family, tmp_reissuedId_order, tmp_reissuedId_class, tmp_reissuedId_phylum);
			int tmpGroupId_species = this->get_groupId_from_reissuedTaxoId_taxoRank(tmp_reissuedId_species, 8);
			int tmpGroupId_genus = this->get_groupId_from_reissuedTaxoId_taxoRank(tmp_reissuedId_genus, 7);
			int tmpGroupId_family = this->get_groupId_from_reissuedTaxoId_taxoRank(tmp_reissuedId_family, 6);
			int tmpGroupId_order = this->get_groupId_from_reissuedTaxoId_taxoRank(tmp_reissuedId_order, 5);
			int tmpGroupId_class = this->get_groupId_from_reissuedTaxoId_taxoRank(tmp_reissuedId_class, 4);
			int tmpGroupId_phylum = this->get_groupId_from_reissuedTaxoId_taxoRank(tmp_reissuedId_phylum, 3);
			tmpGroupIdVec.push_back(tmpGroupId_species);
			tmpGroupIdVec.push_back(tmpGroupId_genus);
			tmpGroupIdVec.push_back(tmpGroupId_family);
			tmpGroupIdVec.push_back(tmpGroupId_order);
			tmpGroupIdVec.push_back(tmpGroupId_class);
			tmpGroupIdVec.push_back(tmpGroupId_phylum);
			for(int tmp2 = 0; tmp2 < tmpGroupIdVec.size(); tmp2++)
			{
				int tmpGroupId_1 = tmpGroupIdVec[tmp2];
				for(int tmp3 = tmp2 + 1; tmp3 < tmpGroupIdVec.size(); tmp3++)
				{
					int tmpGroupId_2 = tmpGroupIdVec[tmp3];
					groupIdPair_inLineOrNot_array[tmpGroupId_1 * total_group_num + tmpGroupId_2] = 1;
					groupIdPair_inLineOrNot_array[tmpGroupId_2 * total_group_num + tmpGroupId_1] = 1;
				}
			}
		}
	}

	bool groupIdPair_inLineOrNot(int tmpGroupId_1, int tmpGroupId_2)
	{
		if((tmpGroupId_1 < 1)||(tmpGroupId_1 > total_taxo_set_num + 1)
			||(tmpGroupId_2 < 1)||(tmpGroupId_2 > total_taxo_set_num + 1))
			return false;
		return (groupIdPair_inLineOrNot_array[tmpGroupId_1 * total_group_num + tmpGroupId_2] == 1);
	}

	bool taxoIdRankPairPair_inLineOrNot(int tmpTaxoReissuedId_1, int tmpTaxoRank_1,
		int tmpTaxoReissuedId_2, int tmpTaxoRank_2)
	{
		int tmpGroupId_1 = this->get_groupId_from_reissuedTaxoId_taxoRank(
			tmpTaxoReissuedId_1, tmpTaxoRank_1);
		int tmpGroupId_2 = this->get_groupId_from_reissuedTaxoId_taxoRank(
			tmpTaxoReissuedId_2, tmpTaxoRank_2);
		return (groupIdPair_inLineOrNot_array[tmpGroupId_1 * total_group_num + tmpGroupId_2] == 1);
	}

	void initiate_multiTaxoRank(BacterialTaxo_Info& bacterialTaxoInfo)
	{
		total_taxo_set_num = bacterialTaxoInfo.return_taxo_num_total();
		total_group_num = total_taxo_set_num + 2; 
		int species_num_total = bacterialTaxoInfo.return_species_num();
		int genus_num_total = bacterialTaxoInfo.return_genus_num();
		int family_num_total = bacterialTaxoInfo.return_family_num();
		int order_num_total = bacterialTaxoInfo.return_order_num();
		int class_num_total = bacterialTaxoInfo.return_class_num();
		int phylum_num_total = bacterialTaxoInfo.return_phylum_num();
		// species, taxo_set_Id_range_specificRank[0], rank = 8
		groupId_range_eachRank_vec.push_back(pair<int,int>(
			1, 
			species_num_total));
		// genus, taxo_set_Id_range_specificRank[1], rank = 7
		groupId_range_eachRank_vec.push_back(pair<int,int>(
			species_num_total + 1, 
			species_num_total + genus_num_total));
		// family, taxo_set_Id_range_specificRank[2], rank = 6
		groupId_range_eachRank_vec.push_back(pair<int,int>(
			species_num_total + genus_num_total + 1, 
			species_num_total + genus_num_total + family_num_total));
		// order, taxo_set_Id_range_specificRank[3], rank = 5
		groupId_range_eachRank_vec.push_back(pair<int,int>(
			species_num_total + genus_num_total + family_num_total + 1, 
			species_num_total + genus_num_total + family_num_total + order_num_total));
		// class, taxo_set_Id_range_specificRank[4], rank = 4
		groupId_range_eachRank_vec.push_back(pair<int,int>(
			species_num_total + genus_num_total + family_num_total + order_num_total + 1, 
			species_num_total + genus_num_total + family_num_total + order_num_total + class_num_total));				
		// phylum, taxo_set_Id_range_specificRank[5], rank = 3
		groupId_range_eachRank_vec.push_back(pair<int,int>(
			species_num_total + genus_num_total + family_num_total + order_num_total + class_num_total + 1, 
			species_num_total + genus_num_total + family_num_total + order_num_total + class_num_total + phylum_num_total));	
		
		// initiate groupId_to_taxoSetIdRankPair_vec
		for(int tmp = 0; tmp <= total_taxo_set_num + 1; tmp++)
			groupId_to_taxoReissuedIdRankPair_vec.push_back(pair<int,int>(-1,-1));
		for(int tmp = 1; tmp <= species_num_total; tmp++) // species
		{
			groupId_to_taxoReissuedIdRankPair_vec[tmp].first = tmp - 1;
			groupId_to_taxoReissuedIdRankPair_vec[tmp].second = 8;
		}
		for(int tmp = species_num_total + 1; 
			tmp <= species_num_total + genus_num_total; tmp++) // genus
		{
			groupId_to_taxoReissuedIdRankPair_vec[tmp].first = tmp - species_num_total - 1;
			groupId_to_taxoReissuedIdRankPair_vec[tmp].second = 7;			
		}
		for(int tmp = species_num_total + genus_num_total + 1; 
			tmp <= species_num_total + genus_num_total + family_num_total; tmp++) // family
		{
			groupId_to_taxoReissuedIdRankPair_vec[tmp].first = tmp - species_num_total - genus_num_total - 1;
			groupId_to_taxoReissuedIdRankPair_vec[tmp].second = 6;
		}
		for(int tmp = species_num_total + genus_num_total + family_num_total + 1; 
			tmp <= species_num_total + genus_num_total + family_num_total + order_num_total; tmp++) // order
		{
			groupId_to_taxoReissuedIdRankPair_vec[tmp].first = tmp - species_num_total - genus_num_total - family_num_total - 1;
			groupId_to_taxoReissuedIdRankPair_vec[tmp].second = 5;
		}
		for(int tmp = species_num_total + genus_num_total + family_num_total + order_num_total + 1; 
			tmp <= species_num_total + genus_num_total + family_num_total + order_num_total + class_num_total; tmp++) // class
		{
			groupId_to_taxoReissuedIdRankPair_vec[tmp].first = tmp - species_num_total - genus_num_total - family_num_total - order_num_total - 1;
			groupId_to_taxoReissuedIdRankPair_vec[tmp].second = 4;
		}
		for(int tmp = species_num_total + genus_num_total + family_num_total + order_num_total + phylum_num_total + 1; 
			tmp <= species_num_total + genus_num_total + family_num_total + order_num_total + class_num_total + phylum_num_total; tmp++) // phylum
		{
			groupId_to_taxoReissuedIdRankPair_vec[tmp].first = tmp - species_num_total - genus_num_total - family_num_total - order_num_total - class_num_total - 1;
			groupId_to_taxoReissuedIdRankPair_vec[tmp].second = 3;
		}				
	}

	void initiate_singleTaxoRank(BacterialTaxo_Info& bacterialTaxoInfo, int tmpRank)
	{
		taxo_class_num = bacterialTaxoInfo.return_taxo_num(tmpRank);
		rank = tmpRank;
		if(rank == 3)
			taxo_rank_name = "Phylum";
		else if(rank == 4)
			taxo_rank_name = "Class";
		else if(rank == 5)
			taxo_rank_name = "Order";
		else if(rank == 6)
			taxo_rank_name = "Family";					
		else if(rank == 7)
			taxo_rank_name = "Genus";
		else if(rank == 8)
			taxo_rank_name = "Species";
		else
		{
			cout << "invalid taxo_rank: " << rank << endl;
			exit(1);
		}					
	}

	void initiate(ReissuedGenomeID2TaxoID_Info& inputReissuedGenomeId2TaxoIdInfo, int inputRank)
	{
		genome_num = inputReissuedGenomeId2TaxoIdInfo.return_genome_num();
		if(inputRank == 3)//phylum
		{
			taxo_class_num = inputReissuedGenomeId2TaxoIdInfo.return_phylum_num();			
			rank = 3;
			taxo_rank_name = "Phylum";
			for(int tmp = 0; tmp < genome_num; tmp++)
				rawKmerClassIndex2taxoClassIndexVec.push_back(
					inputReissuedGenomeId2TaxoIdInfo.return_reissuedPhylumId_from_indexInOriGenomeIdVec(tmp));	
		}
		else if(inputRank == 7)
		{
			taxo_class_num = inputReissuedGenomeId2TaxoIdInfo.return_genus_num();			
			rank = 7;
			taxo_rank_name = "Genus";
			for(int tmp = 0; tmp < genome_num; tmp++)
				rawKmerClassIndex2taxoClassIndexVec.push_back(
					inputReissuedGenomeId2TaxoIdInfo.return_reissuedGenusId_from_indexInOriGenomeIdVec(tmp));		
		}
		else if(inputRank == 8)
		{
			taxo_class_num = inputReissuedGenomeId2TaxoIdInfo.return_species_num();			
			rank = 8;
			taxo_rank_name = "Species";
			for(int tmp = 0; tmp < genome_num; tmp++)
				rawKmerClassIndex2taxoClassIndexVec.push_back(
					inputReissuedGenomeId2TaxoIdInfo.return_reissuedSpeciesId_from_indexInOriGenomeIdVec(tmp));			
		}
		else if(inputRank == 9)
		{
			taxo_class_num = inputReissuedGenomeId2TaxoIdInfo.return_genome_num();
			rank = 9;
			taxo_rank_name = "Genome";
			for(int tmp = 0; tmp < genome_num; tmp++)
				rawKmerClassIndex2taxoClassIndexVec.push_back(
					inputReissuedGenomeId2TaxoIdInfo.return_reissuedSpeciesId_from_indexInOriGenomeIdVec(tmp));			
		}
		else
		{
			cout << "invalid rank: " << inputRank << endl;
			exit(1);
		}
		//repetitive_taxo_class_index = taxo_class_num;// + 1;
		genome_num = inputReissuedGenomeId2TaxoIdInfo.return_genome_num();
		cout << "genome_num in TaxoClassAssignment_Info: " << genome_num << endl;
		cout << "taxo_class_num: " << taxo_class_num << endl;
	}

	/*
	void get_taxoReissuedId_discriminativeAndRepetitiveKmerCountPair_from_groupIdCountVec(
		vector< pair<int, pair<int,int> > >& tmpTaxoReissuedId_discriminativeAndRepetitiveKmerCountPair_vec,
		int tmpRepetitiveKmerCount,
		vector< pair<int,int> >& tmpGroupIdCountVec_lowerRank,
		vector< pair<int,int> >& tmpGroupIdCountVec_thisRank,
		vector< pair<int,int> >& tmpGroupIdCountVec_higherRank, int specificRank,
		BacterialTaxo_Info& bacterialTaxoInfo)
	{
		// tmpGroupIdCountVec_thisRank
		for(int tmp = 0; tmp < tmpGroupIdCountVec_thisRank.size(); tmp++)
		{
			int tmpGroupId = tmpGroupIdCountVec_thisRank[tmp].first;
			int tmpGroupCount = tmpGroupIdCountVec_thisRank[tmp].second;
			int tmpReissuedTaxoId, tmpTaxoRank;
			this->get_taxoReissuedIdRank_from_groupId(tmpGroupId, tmpReissuedTaxoId, tmpTaxoRank);
			if(tmpTaxoRank != specificRank)
			{
				cout << "(tmpTaxoRank != tmpRank) in get_taxoId_discriminativeAndRepetitiveKmerCountPair_from_groupIdCountVec" << endl;
				exit(1);
			}
			tmpTaxoReissuedId_discriminativeAndRepetitiveKmerCountPair_vec.push_back(pair<int, pair<int,int> >
				(tmpReissuedTaxoId, pair<int,int>(tmpGroupCount, 0)));
		}

		// tmpGroupIdCountVec_lowerRank
		for(int tmp = 0; tmp < tmpGroupIdCountVec_lowerRank.size(); tmp++)
		{
			int tmpGroupId = tmpGroupIdCountVec_lowerRank[tmp].first;
			int tmpGroupCount = tmpGroupIdCountVec_lowerRank[tmp].second;
			int tmpReissuedTaxoId_lowerRank, tmpTaxoRank_lowerRank;
			this->get_taxoReissuedIdRank_from_groupId(tmpGroupId, tmpReissuedTaxoId_lowerRank, tmpTaxoRank_lowerRank);
			int tmpRawTaxoId_lowerRank = bacterialTaxoInfo.get_rawTaxoId_from_reissuedTaxoId(tmpReissuedTaxoId_lowerRank, tmpTaxoRank_lowerRank);
			if(tmpRawTaxoId_lowerRank != -1)
			{			
				int tmpReissuedTaxoId = bacterialTaxoInfo.get_higherRankTaxoReissuedId_from_lowerRankTaxoReissuedId(
					//tmpReissuedTaxoId, 
					specificRank, tmpReissuedTaxoId_lowerRank, tmpTaxoRank_lowerRank);
				bool tmp_taxoId_exists_bool = false;
				for(int tmp2 = 0; tmp2 < tmpTaxoReissuedId_discriminativeAndRepetitiveKmerCountPair_vec.size(); tmp2++)
				{
					if(tmpReissuedTaxoId == tmpTaxoReissuedId_discriminativeAndRepetitiveKmerCountPair_vec[tmp2].first)
					{
						(tmpTaxoReissuedId_discriminativeAndRepetitiveKmerCountPair_vec[tmp2].second).first += tmpGroupCount;
						tmp_taxoId_exists_bool = true;
						break;
					}
				}
				if(!tmp_taxoId_exists_bool)
					tmpTaxoReissuedId_discriminativeAndRepetitiveKmerCountPair_vec.push_back(pair<int, pair<int,int> >
						(tmpReissuedTaxoId, pair<int,int>(tmpGroupCount, 0)));
			}
		}		

		// tmpGroupIdCountVec_higherRank
		for(int tmp = 0; tmp < tmpGroupIdCountVec_higherRank.size(); tmp++)
		{
			int tmpGroupId = tmpGroupIdCountVec_higherRank[tmp].first;
			int tmpGroupCount = tmpGroupIdCountVec_higherRank[tmp].second;
			// int tmpReissuedTaxoId_higherRank, tmpTaxoRank_higherRank;
			// this->get_taxoReissuedIdRank_from_groupId(tmpGroupId, tmpReissuedTaxoId_higherRank, tmpTaxoRank_higherRank);
			for(int tmp2 = 0; tmp2 < tmpTaxoReissuedId_discriminativeAndRepetitiveKmerCountPair_vec.size(); tmp2++)
			{
				int tmpTaxoReissuedIdInVec = tmpTaxoReissuedId_discriminativeAndRepetitiveKmerCountPair_vec[tmp2].first;
				int tmpGroupIdInVec = this->get_groupId_from_reissuedTaxoId_taxoRank(tmpTaxoReissuedIdInVec, specificRank);
				bool tmpTaxoIdInLineOrNot_bool = this->groupIdPair_inLineOrNot(tmpGroupId, tmpGroupIdInVec);
				if(tmpTaxoIdInLineOrNot_bool)
					((tmpTaxoReissuedId_discriminativeAndRepetitiveKmerCountPair_vec[tmp2].second).second) += tmpGroupCount;
			}
		}

		// tmpRepetitiveKmerCount
		for(int tmp = 0; tmp < tmpTaxoReissuedId_discriminativeAndRepetitiveKmerCountPair_vec.size(); tmp++)
			((tmpTaxoReissuedId_discriminativeAndRepetitiveKmerCountPair_vec[tmp].second).second) += tmpRepetitiveKmerCount;
	}

	void get_taxoClassIndexCountPairVec_from_rawKmerValueVec_SE(vector<valueT>& tmpKmerValueVec,
		int& invalidKmerCount, int& repetitiveKmerCount, int& discriminativeKmerCount,
		vector< pair<int,int> >& taxoClassIndexCountPairVec)
	{
		int tmpInvalidKmerCount = 0, tmpRepetitiveKmerCount = 0, tmpDiscriminativeKmerCount = 0;
		map<int,int> taxoClassId2countMap;
		for(int tmp = 0; tmp < tmpKmerValueVec.size(); tmp ++)
		{
			valueT tmpValue = tmpKmerValueVec[tmp];
			//cout << "tmpValue: " << tmpValue << endl; 
			if((tmpValue < 1)||(tmpValue > genome_num + 1))
				tmpInvalidKmerCount ++;
			else if(tmpValue == genome_num + 1)
				tmpRepetitiveKmerCount ++;
			else
			{
				tmpDiscriminativeKmerCount ++;
				int tmpTaxoClassIndex = rawKmerClassIndex2taxoClassIndexVec[tmpValue - 1];
				//cout << "tmpTaxoClassIndex: " << tmpTaxoClassIndex << endl;
				map< int,int >::iterator tmpIter = taxoClassId2countMap.find(tmpTaxoClassIndex);
				if(tmpIter == taxoClassId2countMap.end()) // not found
					taxoClassId2countMap.insert(pair<int,int>(tmpTaxoClassIndex, 1));
				else // found
					(tmpIter->second) ++;
			}
		}
		for(map<int,int>::iterator tmpIter = taxoClassId2countMap.begin();
			tmpIter != taxoClassId2countMap.end(); tmpIter ++)
		{
			int tmpTaxoClassIndex = tmpIter->first;
			int tmpTaxoClassCount = tmpIter->second;
			taxoClassIndexCountPairVec.push_back(pair<int,int>(tmpTaxoClassIndex, tmpTaxoClassCount));
		}
		invalidKmerCount = tmpInvalidKmerCount;
		repetitiveKmerCount = tmpRepetitiveKmerCount;
		discriminativeKmerCount = tmpDiscriminativeKmerCount;
	}

	void get_taxoClassIndexCountPairVec_from_rawKmerValueVec_PE(
		vector<valueT>& tmpKmerValueVec_1, vector<valueT>& tmpKmerValueVec_2,
		int& invalidKmerCount, int& repetitiveKmerCount, int& discriminativeKmerCount,
		vector< pair<int,int> >& taxoClassIndexCountPairVec)
	{
		int tmpInvalidKmerCount = 0, tmpRepetitiveKmerCount = 0, tmpDiscriminativeKmerCount = 0;
		map<int,int> taxoClassId2countMap;
		for(int tmp = 0; tmp < tmpKmerValueVec_1.size(); tmp ++)
		{
			valueT tmpValue = tmpKmerValueVec_1[tmp];
			//cout << endl << "tmpValue: " << tmpValue << endl; 
			if((tmpValue < 1)||(tmpValue > genome_num + 1))
				tmpInvalidKmerCount ++;
			else if(tmpValue == genome_num + 1)
				tmpRepetitiveKmerCount ++;
			else
			{
				tmpDiscriminativeKmerCount ++;
				int tmpTaxoClassIndex = rawKmerClassIndex2taxoClassIndexVec[tmpValue - 1];
				//cout << "discriminative Kmer !" << endl;
				//cout << "tmpTaxoClassIndex: " << tmpTaxoClassIndex << endl; 
				map<int,int>::iterator tmpIter = taxoClassId2countMap.find(tmpTaxoClassIndex);
				if(tmpIter == taxoClassId2countMap.end()) // not found
					taxoClassId2countMap.insert(pair<int,int>(tmpTaxoClassIndex, 1));
				else // found
					(tmpIter->second) ++;
			}
		}
		for(int tmp = 0; tmp < tmpKmerValueVec_2.size(); tmp ++)
		{
			valueT tmpValue = tmpKmerValueVec_2[tmp];
			//cout << endl << "tmpValue: " << tmpValue << endl; 
			if((tmpValue < 1)||(tmpValue > genome_num + 1))
				tmpInvalidKmerCount ++;
			else if(tmpValue == genome_num + 1)
				tmpRepetitiveKmerCount ++;
			else
			{
				tmpDiscriminativeKmerCount ++;				
				int tmpTaxoClassIndex = rawKmerClassIndex2taxoClassIndexVec[tmpValue - 1];
				//cout << "discriminative Kmer !" << endl;
				//cout << "tmpTaxoClassIndex: " << tmpTaxoClassIndex << endl; 
				map<int,int>::iterator tmpIter = taxoClassId2countMap.find(tmpTaxoClassIndex);
				if(tmpIter == taxoClassId2countMap.end()) // not found
					taxoClassId2countMap.insert(pair<int,int>(tmpTaxoClassIndex, 1));
				else // found
					(tmpIter->second) ++;
			}
		}		
		for(map<int,int>::iterator tmpIter = taxoClassId2countMap.begin();
			tmpIter != taxoClassId2countMap.end(); tmpIter ++)
		{
			int tmpTaxoClassIndex = tmpIter->first;
			int tmpTaxoClassCount = tmpIter->second;
			taxoClassIndexCountPairVec.push_back(pair<int,int>(tmpTaxoClassIndex, tmpTaxoClassCount));
		}
		invalidKmerCount = tmpInvalidKmerCount;
		repetitiveKmerCount = tmpRepetitiveKmerCount;
		discriminativeKmerCount = tmpDiscriminativeKmerCount;
	}

	void get_rawGroupIdCountPairVec_from_rawKmerValueVec_multiTaxoRank_SE(
		vector<valueT>& tmpKmerValueVec,
		int& invalidKmerCount, int& repetitiveKmerCount, int& discriminativeKmerCount,
		vector< pair<int,int> >& rawGroupIdCountPairVec)
	{
		int tmpInvalidKmerCount = 0, tmpRepetitiveKmerCount = 0, tmpDiscriminativeKmerCount = 0;
		map<int,int> groupId2countMap;
		for(int tmp = 0; tmp < tmpKmerValueVec.size(); tmp ++)
		{
			valueT tmpValue = tmpKmerValueVec[tmp];
			if((tmpValue < 1)||(tmpValue > total_taxo_set_num + 1))
				tmpInvalidKmerCount ++;
			else if(tmpValue == total_taxo_set_num + 1)
				tmpRepetitiveKmerCount ++;
			else
			{
				tmpDiscriminativeKmerCount ++;
				int tmpGroupId = tmpValue;
				map< int,int >::iterator tmpIter = groupId2countMap.find(tmpGroupId);
				if(tmpIter == groupId2countMap.end()) // not found
					groupId2countMap.insert(pair<int,int>(tmpGroupId, 1));
				else // found
					(tmpIter->second) ++;
			}
		}
		for(map<int,int>::iterator tmpIter = groupId2countMap.begin();
			tmpIter != groupId2countMap.end(); tmpIter ++)
		{
			int tmpGroupId = tmpIter->first;
			int tmpGroupCount = tmpIter->second;
			rawGroupIdCountPairVec.push_back(pair<int,int>(tmpGroupId, tmpGroupCount));
		}
		invalidKmerCount = tmpInvalidKmerCount;
		repetitiveKmerCount = tmpRepetitiveKmerCount;
		discriminativeKmerCount = tmpDiscriminativeKmerCount;	
	}

	void get_rawGroupIdCountPairVec_from_rawKmerValueVec_multiTaxoRank_PE(
		vector<valueT>& tmpKmerValueVec_1, vector<valueT>& tmpKmerValueVec_2,
		int& invalidKmerCount, int& repetitiveKmerCount, int& discriminativeKmerCount,
		vector< pair<int,int> >& rawGroupIdCountPairVec)
	{
		int tmpInvalidKmerCount = 0, tmpRepetitiveKmerCount = 0, tmpDiscriminativeKmerCount = 0;
		map<int,int> groupId2countMap;
		for(int tmp = 0; tmp < tmpKmerValueVec_1.size(); tmp ++)
		{
			valueT tmpValue = tmpKmerValueVec_1[tmp];
			if((tmpValue < 1)||(tmpValue > total_taxo_set_num + 1))
				tmpInvalidKmerCount ++;
			else if(tmpValue == total_taxo_set_num + 1)
				tmpRepetitiveKmerCount ++;
			else
			{
				tmpDiscriminativeKmerCount ++;
				int tmpGroupId = tmpValue;
				map< int,int >::iterator tmpIter = groupId2countMap.find(tmpGroupId);
				if(tmpIter == groupId2countMap.end()) // not found
					groupId2countMap.insert(pair<int,int>(tmpGroupId, 1));
				else // found
					(tmpIter->second) ++;
			}
		}
		for(int tmp = 0; tmp < tmpKmerValueVec_2.size(); tmp ++)
		{
			valueT tmpValue = tmpKmerValueVec_2[tmp];
			if((tmpValue < 1)||(tmpValue > total_taxo_set_num + 1))
				tmpInvalidKmerCount ++;
			else if(tmpValue == total_taxo_set_num + 1)
				tmpRepetitiveKmerCount ++;
			else
			{
				tmpDiscriminativeKmerCount ++;
				int tmpGroupId = tmpValue;
				map< int,int >::iterator tmpIter = groupId2countMap.find(tmpGroupId);
				if(tmpIter == groupId2countMap.end()) // not found
					groupId2countMap.insert(pair<int,int>(tmpGroupId, 1));
				else // found
					(tmpIter->second) ++;
			}
		}		
		for(map<int,int>::iterator tmpIter = groupId2countMap.begin();
			tmpIter != groupId2countMap.end(); tmpIter ++)
		{
			int tmpGroupId = tmpIter->first;
			int tmpGroupCount = tmpIter->second;
			rawGroupIdCountPairVec.push_back(pair<int,int>(tmpGroupId, tmpGroupCount));
		}
		invalidKmerCount = tmpInvalidKmerCount;
		repetitiveKmerCount = tmpRepetitiveKmerCount;
		discriminativeKmerCount = tmpDiscriminativeKmerCount;	
	}

	void get_taxoReissuedIdCountPairVec_from_rawGroupIdCountPairVec_multiTaxoRank_specificRank(
		vector< pair<int,int> >& rawGroupIdCountPairVec, int specificRank,
		vector< pair<int,int> >& groupIdCountVec_lowerRank,
		vector< pair<int,int> >& groupIdCountVec_thisRank,
		vector< pair<int,int> >& groupIdCountVec_higherRank)
	{
		for(int tmp = 0; tmp < rawGroupIdCountPairVec.size(); tmp++)
		{
			int tmpGroupId = rawGroupIdCountPairVec[tmp].first;
			int tmpTaxoRank, tmpTaxoReissuedId;
			this->get_taxoReissuedIdRank_from_groupId(tmpGroupId, tmpTaxoReissuedId, tmpTaxoRank);
			if(tmpTaxoRank == specificRank)
				groupIdCountVec_thisRank.push_back(rawGroupIdCountPairVec[tmp]);
			else if(tmpTaxoRank < specificRank)
				groupIdCountVec_higherRank.push_back(rawGroupIdCountPairVec[tmp]);
			else if(tmpTaxoRank > specificRank)
				groupIdCountVec_lowerRank.push_back(rawGroupIdCountPairVec[tmp]);
			else
			{}
		}
	}

	void get_rawTaxoIdCountPairVec_from_rawKmerValueVec_PE(
		vector<valueT>& tmpKmerValueVec_1, vector<valueT>& tmpKmerValueVec_2,
		int& invalidKmerCount, int& repetitiveKmerCount, int& discriminativeKmerCount,
		vector< pair<int,int> >& rawTaxoIdCountPairVec)
	{
		int tmpInvalidKmerCount = 0, tmpRepetitiveKmerCount = 0, tmpDiscriminativeKmerCount = 0;
		map<int,int> taxoId2countMap;
		for(int tmp = 0; tmp < tmpKmerValueVec_1.size(); tmp ++)
		{
			valueT tmpValue = tmpKmerValueVec_1[tmp];
			if((tmpValue < 1)||(tmpValue > taxo_class_num + 1))
				tmpInvalidKmerCount ++;
			else if(tmpValue == taxo_class_num + 1)
				tmpRepetitiveKmerCount ++;
			else
			{
				tmpDiscriminativeKmerCount ++;
				int tmpTaxoId = tmpValue;
				map< int,int >::iterator tmpIter = taxoId2countMap.find(tmpTaxoId);
				if(tmpIter == taxoId2countMap.end()) // not found
					taxoId2countMap.insert(pair<int,int>(tmpTaxoId, 1));
				else // found
					(tmpIter->second) ++;
			}
		}
		for(int tmp = 0; tmp < tmpKmerValueVec_2.size(); tmp ++)
		{
			valueT tmpValue = tmpKmerValueVec_2[tmp];
			if((tmpValue < 1)||(tmpValue > taxo_class_num + 1))
				tmpInvalidKmerCount ++;
			else if(tmpValue == taxo_class_num + 1)
				tmpRepetitiveKmerCount ++;
			else
			{
				tmpDiscriminativeKmerCount ++;
				int tmpTaxoId = tmpValue;
				map< int,int >::iterator tmpIter = taxoId2countMap.find(tmpTaxoId);
				if(tmpIter == taxoId2countMap.end()) // not found
					taxoId2countMap.insert(pair<int,int>(tmpTaxoId, 1));
				else // found
					(tmpIter->second) ++;
			}
		}		
		for(map<int,int>::iterator tmpIter = taxoId2countMap.begin();
			tmpIter != taxoId2countMap.end(); tmpIter ++)
		{
			int tmpTaxoId = tmpIter->first;
			int tmpTaxoCount = tmpIter->second;
			rawTaxoIdCountPairVec.push_back(pair<int,int>(tmpTaxoId, tmpTaxoCount));
		}
		invalidKmerCount = tmpInvalidKmerCount;
		repetitiveKmerCount = tmpRepetitiveKmerCount;
		discriminativeKmerCount = tmpDiscriminativeKmerCount;		
	}

	void get_rawTaxoIdCountPairVec_from_rawKmerValueVec_SE(vector<valueT>& tmpKmerValueVec,
		int& invalidKmerCount, int& repetitiveKmerCount, int& discriminativeKmerCount,
		vector< pair<int,int> >& rawTaxoIdCountPairVec)
	{
		int tmpInvalidKmerCount = 0, tmpRepetitiveKmerCount = 0, tmpDiscriminativeKmerCount = 0;
		map<int,int> taxoId2countMap;
		for(int tmp = 0; tmp < tmpKmerValueVec.size(); tmp ++)
		{
			valueT tmpValue = tmpKmerValueVec[tmp];
			if((tmpValue < 1)||(tmpValue > taxo_class_num + 1))
				tmpInvalidKmerCount ++;
			else if(tmpValue == taxo_class_num + 1)
				tmpRepetitiveKmerCount ++;
			else
			{
				tmpDiscriminativeKmerCount ++;
				int tmpTaxoId = tmpValue;
				map< int,int >::iterator tmpIter = taxoId2countMap.find(tmpTaxoId);
				if(tmpIter == taxoId2countMap.end()) // not found
					taxoId2countMap.insert(pair<int,int>(tmpTaxoId, 1));
				else // found
					(tmpIter->second) ++;
			}
		}
		for(map<int,int>::iterator tmpIter = taxoId2countMap.begin();
			tmpIter != taxoId2countMap.end(); tmpIter ++)
		{
			int tmpTaxoId = tmpIter->first;
			int tmpTaxoCount = tmpIter->second;
			rawTaxoIdCountPairVec.push_back(pair<int,int>(tmpTaxoId, tmpTaxoCount));
		}
		invalidKmerCount = tmpInvalidKmerCount;
		repetitiveKmerCount = tmpRepetitiveKmerCount;
		discriminativeKmerCount = tmpDiscriminativeKmerCount;		
	}

	void get_genomeIdCountPairVec_from_rawKmerValueVec_SE(vector<valueT>& tmpKmerValueVec,
		int& invalidKmerCount, int& repetitiveKmerCount, int& discriminativeKmerCount,
		vector< pair<int,int> >& genomeIdCountPairVec)
	{
		int tmpInvalidKmerCount = 0, tmpRepetitiveKmerCount = 0, tmpDiscriminativeKmerCount = 0;
		map<int,int> genomeId2countMap;
		for(int tmp = 0; tmp < tmpKmerValueVec.size(); tmp ++)
		{
			valueT tmpValue = tmpKmerValueVec[tmp];
			if((tmpValue < 1)||(tmpValue > genome_num + 1))
				tmpInvalidKmerCount ++;
			else if(tmpValue == genome_num + 1)
				tmpRepetitiveKmerCount ++;
			else
			{
				tmpDiscriminativeKmerCount ++;
				int tmpGenomeId = tmpValue;
				map< int,int >::iterator tmpIter = genomeId2countMap.find(tmpGenomeId);
				if(tmpIter == genomeId2countMap.end()) // not found
					genomeId2countMap.insert(pair<int,int>(tmpGenomeId, 1));
				else // found
					(tmpIter->second) ++;
			}
		}
		for(map<int,int>::iterator tmpIter = genomeId2countMap.begin();
			tmpIter != genomeId2countMap.end(); tmpIter ++)
		{
			int tmpGenomeId = tmpIter->first;
			int tmpGenomeCount = tmpIter->second;
			genomeIdCountPairVec.push_back(pair<int,int>(tmpGenomeId, tmpGenomeCount));
		}
		invalidKmerCount = tmpInvalidKmerCount;
		repetitiveKmerCount = tmpRepetitiveKmerCount;
		discriminativeKmerCount = tmpDiscriminativeKmerCount;
	}

	void get_genomeIdCountPairVec_from_rawKmerValueVec_PE(
		vector<valueT>& tmpKmerValueVec_1, vector<valueT>& tmpKmerValueVec_2,
		int& invalidKmerCount, int& repetitiveKmerCount, int& discriminativeKmerCount,
		vector< pair<int,int> >& genomeIdCountPairVec)
	{
		int tmpInvalidKmerCount = 0, tmpRepetitiveKmerCount = 0, tmpDiscriminativeKmerCount = 0;
		map<int,int> genomeId2countMap;
		for(int tmp = 0; tmp < tmpKmerValueVec_1.size(); tmp ++)
		{
			valueT tmpValue = tmpKmerValueVec_1[tmp];
			if((tmpValue < 1)||(tmpValue > genome_num + 1))
				tmpInvalidKmerCount ++;
			else if(tmpValue == genome_num + 1)
				tmpRepetitiveKmerCount ++;
			else
			{
				tmpDiscriminativeKmerCount ++;
				int tmpGenomeId = tmpValue;
				map< int,int >::iterator tmpIter = genomeId2countMap.find(tmpGenomeId);
				if(tmpIter == genomeId2countMap.end()) // not found
					genomeId2countMap.insert(pair<int,int>(tmpGenomeId, 1));
				else // found
					(tmpIter->second) ++;
			}
		}
		for(int tmp = 0; tmp < tmpKmerValueVec_2.size(); tmp ++)
		{
			valueT tmpValue = tmpKmerValueVec_2[tmp];
			if((tmpValue < 1)||(tmpValue > genome_num + 1))
				tmpInvalidKmerCount ++;
			else if(tmpValue == genome_num + 1)
				tmpRepetitiveKmerCount ++;
			else
			{
				tmpDiscriminativeKmerCount ++;
				int tmpGenomeId = tmpValue;
				map< int,int >::iterator tmpIter = genomeId2countMap.find(tmpGenomeId);
				if(tmpIter == genomeId2countMap.end()) // not found
					genomeId2countMap.insert(pair<int,int>(tmpGenomeId, 1));
				else // found
					(tmpIter->second) ++;
			}
		}		
		for(map<int,int>::iterator tmpIter = genomeId2countMap.begin();
			tmpIter != genomeId2countMap.end(); tmpIter ++)
		{
			int tmpGenomeId = tmpIter->first;
			int tmpGenomeCount = tmpIter->second;
			genomeIdCountPairVec.push_back(pair<int,int>(tmpGenomeId, tmpGenomeCount));
		}
		invalidKmerCount = tmpInvalidKmerCount;
		repetitiveKmerCount = tmpRepetitiveKmerCount;
		discriminativeKmerCount = tmpDiscriminativeKmerCount;
	}

	void get_taxoIdCount_best_secondBest(vector< pair<int,int> >& taxoClassIndexCountPairVec,
		int invalidKmerCount, int repetitiveKmerCount, int discriminativeKmerCount,
		int& bestTaxoClass_index, int& bestTaxoClass_count, int& secondBestTaxoClass_index, int& secondBestTaxoClass_count)
	{
		for(int tmp = 0; tmp < taxoClassIndexCountPairVec.size(); tmp++)
		{
			int tmpTaxoClassIndex = taxoClassIndexCountPairVec[tmp].first;
			int tmpTaxoClassCount = taxoClassIndexCountPairVec[tmp].second;
			if(tmpTaxoClassCount > bestTaxoClass_count)
			{
				secondBestTaxoClass_index = bestTaxoClass_index;
				secondBestTaxoClass_count = bestTaxoClass_count;
				bestTaxoClass_index = tmpTaxoClassIndex;
				bestTaxoClass_count = tmpTaxoClassCount;
			}
			else if(tmpTaxoClassCount > secondBestTaxoClass_count)
			{
				secondBestTaxoClass_index = tmpTaxoClassIndex;
				secondBestTaxoClass_count = tmpTaxoClassCount;				
			}
		}		
	}

	void get_taxoReissuedId_specificityAndConfidenceScore_from_discriminativeAndRepetitiveKmerCount(
		vector< pair<int, pair<int,int> > >& tmpTaxoReissuedId_discriminativeAndRepetitiveKmerCountPair_vec,
		vector< pair<int, pair<double,double> > >& tmpTaxoReissuedId_specificityAndConfidenceScore_vec,
		int total_Kmer_count)
	{
		for(int tmp = 0; tmp < tmpTaxoReissuedId_discriminativeAndRepetitiveKmerCountPair_vec.size(); tmp++)
		{
			int tmpTaxoReissuedId = tmpTaxoReissuedId_discriminativeAndRepetitiveKmerCountPair_vec[tmp].first;
			int tmpTaxoCount_specific = (tmpTaxoReissuedId_discriminativeAndRepetitiveKmerCountPair_vec[tmp].second).first;
			int tmpTaxoCount_confident = (tmpTaxoReissuedId_discriminativeAndRepetitiveKmerCountPair_vec[tmp].second).second;
			double tmpTaxo_specificity_score = (double)tmpTaxoCount_specific/(double)total_Kmer_count;
			double tmpTaxo_confidence_score = (double)(tmpTaxoCount_specific+tmpTaxoCount_confident)/(double)total_Kmer_count;
			tmpTaxoReissuedId_specificityAndConfidenceScore_vec.push_back(pair<int, pair<double,double> >(tmpTaxoReissuedId,
				pair<double, double> (tmpTaxo_specificity_score, tmpTaxo_confidence_score)));
		}
	}

	void get_assignedTaxoReissuedId_best_secondBest(
		vector< pair<int, pair<double,double> > >& tmpTaxoReissuedId_specificityAndConfidenceScore_vec, 
		int& tmp_taxo_reissuedId_best, 
		double& tmp_taxo_reissuedId_best_specificity_score, double& tmp_taxo_reissuedId_best_confidence_score,
		int& tmp_taxo_reissuedId_secondBest, 
		double& tmp_taxo_reissuedId_secondBest_specificity_score, double& tmp_taxo_reissuedId_secondBest_confidence_score,
		double assignment_specificity_score_min, double assignment_confidence_score_min)	
	{
		tmp_taxo_reissuedId_best = -1;
		double tmp_taxo_reissuedId_best_score = 0.0;
		tmp_taxo_reissuedId_best_specificity_score = 0.0;
		tmp_taxo_reissuedId_best_confidence_score = 0.0;

		tmp_taxo_reissuedId_secondBest = -1;
		double tmp_taxo_reissuedId_secondBest_score = 0.0;
		tmp_taxo_reissuedId_secondBest_specificity_score = 0.0;
		tmp_taxo_reissuedId_secondBest_confidence_score = 0.0;
		
		for(int tmp = 0; tmp < tmpTaxoReissuedId_specificityAndConfidenceScore_vec.size(); tmp++)
		{
			int tmpTaxo_reissuedId = tmpTaxoReissuedId_specificityAndConfidenceScore_vec[tmp].first;
			double tmpTaxo_score = (tmpTaxoReissuedId_specificityAndConfidenceScore_vec[tmp].second).first
				+ (tmpTaxoReissuedId_specificityAndConfidenceScore_vec[tmp].second).second;
			if(tmpTaxo_score > tmp_taxo_reissuedId_best_score)
			{
				tmp_taxo_reissuedId_secondBest = tmp_taxo_reissuedId_best;
				tmp_taxo_reissuedId_secondBest_score = tmp_taxo_reissuedId_best_score;
				tmp_taxo_reissuedId_secondBest_specificity_score = tmp_taxo_reissuedId_best_specificity_score;
				tmp_taxo_reissuedId_secondBest_confidence_score = tmp_taxo_reissuedId_best_confidence_score;				
				tmp_taxo_reissuedId_best = tmpTaxo_reissuedId;
				tmp_taxo_reissuedId_best_score = tmpTaxo_score;
				tmp_taxo_reissuedId_best_specificity_score = (tmpTaxoReissuedId_specificityAndConfidenceScore_vec[tmp].second).first;
				tmp_taxo_reissuedId_best_confidence_score = (tmpTaxoReissuedId_specificityAndConfidenceScore_vec[tmp].second).second;
			}
			else if(tmpTaxo_score > tmp_taxo_reissuedId_secondBest_score)
			{
				tmp_taxo_reissuedId_secondBest = tmpTaxo_reissuedId;
				tmp_taxo_reissuedId_secondBest_score = tmpTaxo_score;
				tmp_taxo_reissuedId_secondBest_specificity_score = (tmpTaxoReissuedId_specificityAndConfidenceScore_vec[tmp].second).first;
				tmp_taxo_reissuedId_secondBest_confidence_score = (tmpTaxoReissuedId_specificityAndConfidenceScore_vec[tmp].second).second;
			}
			else
			{}
		}

		if(tmp_taxo_reissuedId_best != -1)
		{
			if((tmp_taxo_reissuedId_best_specificity_score < assignment_specificity_score_min)
				||(tmp_taxo_reissuedId_best_confidence_score < assignment_confidence_score_min))
				tmp_taxo_reissuedId_best = -2;
				//tmp_taxo_reissuedId_secondBest = -1;
		}
	}

	int determineAssignedTaxoClassIndex_from_taxoClassIndexCountPairVec(
		double discriminativeKmer_proportion_best_secondBest_diff_min, int discriminativeKmer_count_best_min, 
		int invalidKmerCount, int repetitiveKmerCount, int discriminativeKmerCount,
		vector< pair<int,int> >& taxoClassIndexCountPairVec,
		int& bestTaxoClass_index, int& bestTaxoClass_count, int& secondBestTaxoClass_index, int& secondBestTaxoClass_count)
	{		
		for(int tmp = 0; tmp < taxoClassIndexCountPairVec.size(); tmp++)
		{
			int tmpTaxoClassIndex = taxoClassIndexCountPairVec[tmp].first;
			int tmpTaxoClassCount = taxoClassIndexCountPairVec[tmp].second;
			if(tmpTaxoClassCount > bestTaxoClass_count)
			{
				secondBestTaxoClass_index = bestTaxoClass_index;
				secondBestTaxoClass_count = bestTaxoClass_count;
				bestTaxoClass_index = tmpTaxoClassIndex;
				bestTaxoClass_count = tmpTaxoClassCount;
			}
			else if(tmpTaxoClassCount > secondBestTaxoClass_count)
			{
				secondBestTaxoClass_index = tmpTaxoClassIndex;
				secondBestTaxoClass_count = tmpTaxoClassCount;				
			}
		}
		//if(bestTaxoClass_index == -1)
		//	return -1;
		//else
		if(bestTaxoClass_index != -1)
		{
			//double tmp_discriminativeKmer_best_proportion = (double)bestTaxoClass_count/(double)discriminativeKmerCount;
			//double tmp_discriminativeKmer_secondBest_proportion = (double)secondBestTaxoClass_count/(double)discriminativeKmerCount;
			//double tmp_discriminativeKmer_best_secondBest_diff_divide
			//	= tmp_discriminativeKmer_best_proportion / tmp_discriminativeKmer_secondBest_proportion;
			double tmp_discriminativeKmer_best_secondBest_diff_divide;
			if((secondBestTaxoClass_count == 0)&&(bestTaxoClass_count > 0))
				tmp_discriminativeKmer_best_secondBest_diff_divide = 99;
			else if((secondBestTaxoClass_count == 0)&&(bestTaxoClass_count == 0))
				tmp_discriminativeKmer_best_secondBest_diff_divide = 0;
			else
				tmp_discriminativeKmer_best_secondBest_diff_divide = (double)bestTaxoClass_count/(double)secondBestTaxoClass_count;
			if(!((tmp_discriminativeKmer_best_secondBest_diff_divide >= discriminativeKmer_proportion_best_secondBest_diff_min)
				&&(bestTaxoClass_count >= discriminativeKmer_count_best_min)))
				bestTaxoClass_index = -1;
			//else
			//	return -1;
		}
		return bestTaxoClass_index;
	}

	int determineAssignedTaxoClassIndex_from_genomeIdCountPairVec(double discriminativeKmer_proportion_best_secondBest_diff_min,
		int discriminativeKmer_count_best_min, int invalidKmerCount, int repetitiveKmerCount, int discriminativeKmerCount,
		vector< pair<int,int> >& genomeIdCountPairVec, ReissuedGenomeID2TaxoID_Info& inputReissuedGenomeId2TaxoIdInfo,
		int& bestGenome_Id, int& bestGenome_count, int& secondBestGenome_Id, int& secondBestGenome_count)
	{
		int tmp_bestGenome_id = this->determineAssignedTaxoClassIndex_from_taxoClassIndexCountPairVec(discriminativeKmer_proportion_best_secondBest_diff_min,
				discriminativeKmer_count_best_min, invalidKmerCount, repetitiveKmerCount, discriminativeKmerCount,
				genomeIdCountPairVec,
				bestGenome_Id, bestGenome_count, secondBestGenome_Id, secondBestGenome_count);
		if(tmp_bestGenome_id == -1)
			return -1;
		else 
			return inputReissuedGenomeId2TaxoIdInfo.return_reissuedTaxoId_from_genomeId(tmp_bestGenome_id, rank);
	}*/

	// int determineAssignedTaxoClassIndex_from_rawTaxoIdCountPairVec(
	// 	double discriminativeKmer_proportion_best_secondBest_diff_min, int discriminativeKmer_count_best_min, 
	// 	int invalidKmerCount, int repetitiveKmerCount, int discriminativeKmerCount,
	// 	vector< pair<int,int> >& rawTaxoIdCountPairVec, ReissuedGenomeID2TaxoID_Info& inputReissuedGenomeId2TaxoIdInfo,
	// 	int& bestTaxo_Id, int& bestTaxo_count, int& secondBestTaxo_Id, int& secondBestTaxo_count)
	// {
	// 	int tmp_bestGenome_id = this->determineAssignedTaxoClassIndex_from_taxoClassIndexCountPairVec(
	// 			discriminativeKmer_proportion_best_secondBest_diff_min,
	// 			discriminativeKmer_count_best_min, invalidKmerCount, repetitiveKmerCount, discriminativeKmerCount,
	// 			genomeIdCountPairVec,
	// 			bestGenome_Id, bestGenome_count, secondBestGenome_Id, secondBestGenome_count);
	// 	if(tmp_bestGenome_id == -1)
	// 		return -1;
	// 	else 
	// 		return inputReissuedGenomeId2TaxoIdInfo.return_reissuedTaxoId_from_genomeId(tmp_bestGenome_id, rank);
	// }	
};
#endif