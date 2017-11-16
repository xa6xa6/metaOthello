// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#ifndef FUSIONASSIGNEMNT_INFO_H
#define FUSIONASSIGNMENT_IFNO_H
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

class FusionAssignment_Info
{
private:
	int discriminativeClassId_min;
	int discriminativeClassId_max;
	int repetitiveClassId;
public:
	FusionAssignment_Info()
	{}

	void initiate(int tmpDiscriminativeClassId_min, int tmpDiscriminativeClassId_max, int tmpRepetitiveClassId)
	{
		discriminativeClassId_min = tmpDiscriminativeClassId_min;
		discriminativeClassId_max = tmpDiscriminativeClassId_max;
		repetitiveClassId = tmpRepetitiveClassId;
	}

	bool examineTwoKmerWindowDistance(int window_1_rightBoundary_compatible, int window_2_leftBoundary_compatible, 
		int tmpQuerySeqKmerNum_1, int Kmer_length)
	{
		if((window_1_rightBoundary_compatible > tmpQuerySeqKmerNum_1)||(window_2_leftBoundary_compatible <= tmpQuerySeqKmerNum_1))
		{
			int tmpWindowDistance = window_2_leftBoundary_compatible - window_1_rightBoundary_compatible - 1;
			if((tmpWindowDistance >= Kmer_length - FUSION_GENE_PAIR_DETECTION_BUFFER_MAX)
				&&(tmpWindowDistance <= Kmer_length + FUSION_GENE_PAIR_DETECTION_BUFFER_MAX))
				return true;
			else
				return false;
		}
		else if(((window_2_leftBoundary_compatible - window_1_rightBoundary_compatible - 1) <= (Kmer_length + FUSION_GENE_PAIR_DETECTION_BUFFER_MAX))
			&&(((window_1_rightBoundary_compatible >= tmpQuerySeqKmerNum_1 - FUSION_GENE_PAIR_DETECTION_BUFFER_MAX)
				&&(window_1_rightBoundary_compatible <= tmpQuerySeqKmerNum_1))
				||((window_2_leftBoundary_compatible >= tmpQuerySeqKmerNum_1 + 1)
					&&(window_2_leftBoundary_compatible <= tmpQuerySeqKmerNum_1 + 1 + FUSION_GENE_PAIR_DETECTION_BUFFER_MAX))))
			return true;
		else
			return false;
	}

	int getCompatibleKmerCount_1stWindow(int window_1_rightBoundary_compatible, 
		int fusionAssignment_geneReissuedId_1, vector<valueT>& tmpKmerValueVec)
	{
		int tmpCompatibleKmerCount = 0;
		for(int tmp = 0; tmp <= window_1_rightBoundary_compatible; tmp++)	
		{
			if((tmpKmerValueVec[tmp] == repetitiveClassId)||(tmpKmerValueVec[tmp] == fusionAssignment_geneReissuedId_1))
				tmpCompatibleKmerCount ++;
		}
		return tmpCompatibleKmerCount;
	}

	int getCompatibleKmerCount_2ndWindow(int window_2_leftBoundary_compatible,
		int fusionAssignment_geneReissuedId_2, vector<valueT>& tmpKmerValueVec)
	{
		int tmpCompatibleKmerCount = 0;
		for(int tmp = window_2_leftBoundary_compatible; tmp < tmpKmerValueVec.size(); tmp++)
		{
			if((tmpKmerValueVec[tmp] == repetitiveClassId)||(tmpKmerValueVec[tmp] == fusionAssignment_geneReissuedId_2))
				tmpCompatibleKmerCount ++;
		}
		return tmpCompatibleKmerCount;
	}

	int updateWindowRightBoundaryWithRepeatKmer(int windowRightBoundary, vector<valueT>& tmpKmerValueVec)
	{
		int windowRightBoundary_new = windowRightBoundary;
		for(int tmpIndex = windowRightBoundary + 1; tmpIndex < tmpKmerValueVec.size(); tmpIndex++)
		{
			if(tmpKmerValueVec[tmpIndex] == repetitiveClassId)
				windowRightBoundary_new ++;
			else
				break;
		}
		return windowRightBoundary_new;
	}

	int updateWindowLeftBoundaryWithRepeatKmer(int windowLeftBoundary, vector<valueT>& tmpKmerValueVec)
	{
		int windowLeftBoundary_new = windowLeftBoundary;
		for(int tmpIndex = windowLeftBoundary - 1; tmpIndex >= 0; tmpIndex--)
		{
			if(tmpKmerValueVec[tmpIndex] == repetitiveClassId)
				windowLeftBoundary_new --;
			else
				break;			
		}
		return windowLeftBoundary_new;
	}

	void getFusionAssignment(int& fusionAssignment_geneReissuedId_1, int& fusionAssignment_geneReissuedId_2,
		
		int& kmerCount_geneReissuedId_1, int& kmerCount_geneReissuedId_2, 
		int& fusionSite, int& kmerWindowDistance,
		
		int& kmerCount_geneReissuedId_1_compatible, int& kmerCount_geneReissuedId_2_compatible, 
		int& fusionSite_compatible, int& kmerWindowDistance_compatible,
		
		vector<valueT>& tmpKmerValueVec, int tmpQuerySeqKmerNum_1, int Kmer_length)
	{
		int invalidKmerCount = 0;
		int repetitiveKmerCount = 0;
		int discriminativeKmerCount = 0;
		map<int, int> geneId2countMap;
		for(int tmp = 0; tmp < tmpKmerValueVec.size(); tmp++)
		{
			int tmpKmerValue = (int)tmpKmerValueVec[tmp];
			//cout << "tmpKmerValue: " << tmpKmerValue << endl;
			if((tmpKmerValue >= discriminativeClassId_min)&&(tmpKmerValue <= discriminativeClassId_max))
			{
				discriminativeKmerCount ++;
				map<int, int>::iterator tmpIter = geneId2countMap.find(tmpKmerValue);
				if(tmpIter != geneId2countMap.end()) // found
					(tmpIter->second) ++;
				else // not found, new kmer value
					geneId2countMap.insert(pair<int,int>(tmpKmerValue, 1));
			}
			else if(tmpKmerValue == repetitiveClassId)
				repetitiveKmerCount ++;
			else
				invalidKmerCount ++;
		}

		int gene_reissuedId_best = -1;
		int gene_count_best = 0;
		int gene_reissuedId_secondBest = -1;
		int gene_count_secondBest = 0;
		int gene_reissuedId_thirdBest = -1;
		int gene_count_thirdBest = 0;		
		for(map<int, int>::iterator tmpIter = geneId2countMap.begin(); tmpIter != geneId2countMap.end(); tmpIter ++)
		{
			int tmpGene_reissuedId = tmpIter->first;
			int tmpGene_count = tmpIter->second;
			if(tmpGene_count > gene_count_best)
			{
				gene_reissuedId_thirdBest = gene_reissuedId_secondBest;
				gene_count_thirdBest = gene_count_secondBest;
				gene_reissuedId_secondBest = gene_reissuedId_best;
				gene_count_secondBest = gene_count_best;
				gene_reissuedId_best = tmpGene_reissuedId;
				gene_count_best = tmpGene_count;
			}	
			else if(tmpGene_count > gene_count_secondBest)
			{
				gene_reissuedId_thirdBest = gene_reissuedId_secondBest;
				gene_count_thirdBest = gene_count_secondBest;				
				gene_reissuedId_secondBest = tmpGene_reissuedId;
				gene_count_secondBest = tmpGene_count;
			}
			else if(tmpGene_count > gene_count_thirdBest)
			{
				gene_reissuedId_thirdBest = tmpGene_reissuedId;
				gene_count_thirdBest = tmpGene_count;
			}
			else
			{}
		}
		//cout << "gene_reissuedId_best: " << gene_reissuedId_best << endl;
		//cout << "gene_count_best: " << gene_count_best << endl;
		//cout << "gene_reissuedId_secondBest: " << gene_reissuedId_secondBest << endl;
		//cout << "gene_count_secondBest: " << gene_count_secondBest << endl;
		//cout << "gene_reissuedId_thirdBest: " << gene_reissuedId_thirdBest << endl;
		//cout << "gene_count_thirdBest: " << gene_count_thirdBest << endl;				
		// initiate fusionAssignment_geneReissuedId
		fusionAssignment_geneReissuedId_1 = -1;
		fusionAssignment_geneReissuedId_2 = -1;

		if((gene_reissuedId_best >= 0)&&(gene_reissuedId_secondBest >= 0)&&(gene_count_secondBest != gene_count_thirdBest))
		{
			int best_gene_index_min = -1;
			int best_gene_index_max = -1;
			int secondBest_gene_index_min = -1;
			int secondBest_gene_index_max = -1;			
			for(int tmp = 0; tmp < tmpKmerValueVec.size(); tmp++)
			{
				int tmpKmerValue = tmpKmerValueVec[tmp];
				if(tmpKmerValue == gene_reissuedId_best)
				{
					if(best_gene_index_min == -1)
						best_gene_index_min = tmp;
					best_gene_index_max = tmp;
				}
				else if(tmpKmerValue == gene_reissuedId_secondBest)
				{
					if(secondBest_gene_index_min == -1)
						secondBest_gene_index_min = tmp;
					secondBest_gene_index_max = tmp;
				}
				else
				{}
			}

			// updates window boundaries
			if(best_gene_index_min > secondBest_gene_index_max)
			{
				// updates window right boundary
				int window_1_rightBoundary_compatible
					= this->updateWindowRightBoundaryWithRepeatKmer(secondBest_gene_index_max, tmpKmerValueVec);
				// updates window left boundary
				int window_2_leftBoundary_compatible
					= this->updateWindowLeftBoundaryWithRepeatKmer(best_gene_index_min, tmpKmerValueVec);
				if((window_1_rightBoundary_compatible < window_2_leftBoundary_compatible)
					&&(gene_count_best >= FUSION_GENE_PAIR_DETECTION_ANCHOR_KMER_NUM_MIN)
					&&(gene_count_secondBest >= FUSION_GENE_PAIR_DETECTION_ANCHOR_KMER_NUM_MIN))
				{
					bool tmpTwoKmerWindowDistance_valid_bool = this->examineTwoKmerWindowDistance(
						window_1_rightBoundary_compatible, window_2_leftBoundary_compatible, tmpQuerySeqKmerNum_1, Kmer_length);
					if(!tmpTwoKmerWindowDistance_valid_bool)
						return;
					fusionAssignment_geneReissuedId_1 = gene_reissuedId_secondBest;
					fusionAssignment_geneReissuedId_2 = gene_reissuedId_best;

					kmerCount_geneReissuedId_1 = gene_count_secondBest;
					kmerCount_geneReissuedId_2 = gene_count_best;
					fusionSite = this->getFusionSite(secondBest_gene_index_max, best_gene_index_min, 
						tmpQuerySeqKmerNum_1, Kmer_length); //(best_gene_index_min + secondBest_gene_index_max + Kmer_length + 2)/2;
					kmerWindowDistance = best_gene_index_min - secondBest_gene_index_max;

					kmerCount_geneReissuedId_1_compatible = this->getCompatibleKmerCount_1stWindow(window_1_rightBoundary_compatible, 
						fusionAssignment_geneReissuedId_1, tmpKmerValueVec);
					kmerCount_geneReissuedId_2_compatible = this->getCompatibleKmerCount_2ndWindow(window_2_leftBoundary_compatible,
						fusionAssignment_geneReissuedId_2, tmpKmerValueVec);
					fusionSite_compatible = this->getFusionSite(window_1_rightBoundary_compatible, window_2_leftBoundary_compatible,
						tmpQuerySeqKmerNum_1, Kmer_length); //(window_1_rightBoundary_compatible + window_2_leftBoundary_compatible + Kmer_length + 2)/2;
					kmerWindowDistance_compatible = window_2_leftBoundary_compatible - window_1_rightBoundary_compatible;
				}
			}
			else if(best_gene_index_max < secondBest_gene_index_min)
			{
				// updates window right boundary
				int window_1_rightBoundary_compatible = 
					this->updateWindowRightBoundaryWithRepeatKmer(best_gene_index_max, tmpKmerValueVec);
				// updates window left boundary
				int window_2_leftBoundary_compatible = 
					this->updateWindowLeftBoundaryWithRepeatKmer(secondBest_gene_index_min, tmpKmerValueVec);
				if(//(window_1_rightBoundary_compatible + 21 >= window_2_leftBoundary_compatible)&&
					(window_1_rightBoundary_compatible < window_2_leftBoundary_compatible)
					&&(gene_count_best >= FUSION_GENE_PAIR_DETECTION_ANCHOR_KMER_NUM_MIN)
					&&(gene_count_secondBest >= FUSION_GENE_PAIR_DETECTION_ANCHOR_KMER_NUM_MIN))
				{
					bool tmpTwoKmerWindowDistance_valid_bool = this->examineTwoKmerWindowDistance(
						window_1_rightBoundary_compatible, window_2_leftBoundary_compatible, tmpQuerySeqKmerNum_1, Kmer_length);
					if(!tmpTwoKmerWindowDistance_valid_bool)
						return;
					fusionAssignment_geneReissuedId_1 = gene_reissuedId_best;
					fusionAssignment_geneReissuedId_2 = gene_reissuedId_secondBest;
				
					kmerCount_geneReissuedId_1 = gene_count_best;
					kmerCount_geneReissuedId_2 = gene_count_secondBest;
					fusionSite = this->getFusionSite(best_gene_index_max, secondBest_gene_index_min,
						tmpQuerySeqKmerNum_1, Kmer_length); //(best_gene_index_max + secondBest_gene_index_min + Kmer_length + 2)/2;
					kmerWindowDistance = secondBest_gene_index_min - best_gene_index_max;

					kmerCount_geneReissuedId_1_compatible = this->getCompatibleKmerCount_1stWindow(window_1_rightBoundary_compatible, 
						fusionAssignment_geneReissuedId_1, tmpKmerValueVec);
					kmerCount_geneReissuedId_2_compatible = this->getCompatibleKmerCount_2ndWindow(window_2_leftBoundary_compatible,
						fusionAssignment_geneReissuedId_2, tmpKmerValueVec);
					fusionSite_compatible = this->getFusionSite(window_1_rightBoundary_compatible, window_2_leftBoundary_compatible,
						tmpQuerySeqKmerNum_1, Kmer_length); //(window_1_rightBoundary_compatible + window_2_leftBoundary_compatible + Kmer_length + 2)/2;
					kmerWindowDistance_compatible = window_2_leftBoundary_compatible - window_1_rightBoundary_compatible;
				}				
			}
			else
			{}
		}
	}

	int getFusionSite(int window_1_rightBoundary, int window_2_leftBoundary, int tmpQuerySeqKmerNum_1, int Kmer_length)
	{
		return (window_1_rightBoundary + Kmer_length - 1 + window_2_leftBoundary)/2 + 1;
		// if((window_1_rightBoundary > tmpQuerySeqKmerNum_1)||(window_2_leftBoundary <= tmpQuerySeqKmerNum_1))
		// 	return (window_1_rightBoundary + Kmer_length - 1 + window_2_leftBoundary)/2;
		// else
		// {
		// 	int distanceToRead1End_window1 = tmpQuerySeqKmerNum_1 - window_1_rightBoundary;
		// 	int distanceToRead2Start_window2 = window_2_leftBoundary - (tmpQuerySeqKmerNum_1 + 1);
		// 	if((distanceToRead1End_window1 <= FUSION_GENE_PAIR_DETECTION_BUFFER_MAX)
		// 		&&(distanceToRead2Start_window2 <= FUSION_GENE_PAIR_DETECTION_BUFFER_MAX))
		// 		return tmpQuerySeqKmerNum_1 + Kmer_length - 1;
		// 	else if(distanceToRead1End_window1 >= distanceToRead2Start_window2)

		// }
		// if((window_1_rightBoundary <= tmpQuerySeqKmerNum_1)&&(window_2_leftBoundary > tmpQuerySeqKmerNum_1))
		// 	return (window_1_rightBoundary + window_2_leftBoundary)/2;
		// else if(window_1_rightBoundary > tmpQuerySeqKmerNum_1)
		// {}
		// else
		// {

		// }
	}
};

#endif