// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#ifndef REGIONASSIGNMENT_INFO_H
#define REGIONASSIGNMENT_INFO_H
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
#include "../../query/general/queryConstantDef.h"

using namespace std;

class RegionAssignment_Info
{
private:
	int validRepetitiveKmerClassIndex;
	int invalidKmerClassIndex;

	int validRepetitiveClassKmerCount;
	int validDiscrimitiveClassKmerCount_total;
   	int invalidClassKmerCount;
   	int totalKmerCount;

	int validDiscrimitiveClassIndex_best;
	int validDiscrimitiveClassKmerCount_best;

	int validDiscrimitiveClassIndex_secondBest;
	int validDiscrimitiveClassKmerCount_secondBest;

	double confidence_score_best;
	double uniqueness_score_best;
public:
	RegionAssignment_Info()
	{
		validRepetitiveKmerClassIndex = -1;
		invalidKmerClassIndex = -1;

		validRepetitiveClassKmerCount = 0;
		validDiscrimitiveClassKmerCount_total = 0;
		invalidClassKmerCount = 0;
		totalKmerCount = 0;

		validDiscrimitiveClassIndex_best = -1;
		validDiscrimitiveClassKmerCount_best = 0;
	
		validDiscrimitiveClassIndex_secondBest = -1;
		validDiscrimitiveClassKmerCount_secondBest = 0;
	}

	int return_validDiscrimitiveClassIndex_best()
	{
		return validDiscrimitiveClassIndex_best;
	}

	int return_validDiscrimitiveClassIndex_secondBest()
	{
		return validDiscrimitiveClassIndex_secondBest;
	}

	int return_validDiscrimitiveClassKmerCount_best()
	{
		return validDiscrimitiveClassKmerCount_best;
	}

	int return_validDiscrimitiveClassKmerCount_secondBest()
	{
		return validDiscrimitiveClassKmerCount_secondBest;
	}

	int return_validDiscrimitiveClassKmerCount_total()
	{
		return validDiscrimitiveClassKmerCount_total;
	}

	int return_validClassKmerCount_repetitive()
	{
		return validRepetitiveClassKmerCount;
	}	

	int return_invalidClassKmerCount()
	{
		return invalidClassKmerCount;
	}

	double return_confidence_score_best()
	{
		return confidence_score_best;
	}

	double return_uniqueness_score_best()
	{
		return uniqueness_score_best;
	}

	void initiate(int tmpValidRepetitiveKmerClassIndex, int tmpInvalidKmerClassIndex) 
	// though Kmers assigned to other classes can also be alien Kmers
	{
		validRepetitiveKmerClassIndex = tmpValidRepetitiveKmerClassIndex;
		invalidKmerClassIndex = tmpInvalidKmerClassIndex;
	}

	void getMostKmerCountRegion_SE_1st_2nd(vector< pair<valueT,int> >& tmpQuerySeq_classIndexCountPairVec_1)
	{
		int tmpValidRepetitiveClassKmerCount = 0;
		int tmpValidDiscrimitiveClassKmerCount_total = 0;
    	int tmpInvalidClassKmerCount = 0;
    	int tmpTotalKmerCount = 0;

    	int tmpValidDiscrimitiveClassKmerCount_best = 0;
    	int tmpValidDiscrimitiveClassIndex_best = -1;
    	
    	//get the best validDiscrimitiveClassKmerCount and validDiscrimitiveClassIndex
    	vector< pair<valueT,int> > tmpValidDiscrimitiveClassIndexCountPairVec;
    	for(int tmp = 0; tmp < tmpQuerySeq_classIndexCountPairVec_1.size(); tmp++)
    	{
    		int tmpClassCount = tmpQuerySeq_classIndexCountPairVec_1[tmp].second;
    		int tmpClassIndex = (int)tmpQuerySeq_classIndexCountPairVec_1[tmp].first;
    		tmpTotalKmerCount += tmpClassCount;
    		if(tmpClassIndex == invalidKmerClassIndex) // invalid class
    			tmpInvalidClassKmerCount += tmpClassCount;
    		else if(tmpClassIndex == validRepetitiveKmerClassIndex) // valid repetitive class
    			tmpValidRepetitiveClassKmerCount += tmpClassCount;
    		else // valid discrimitive class
    		{
    			tmpValidDiscrimitiveClassIndexCountPairVec.push_back(tmpQuerySeq_classIndexCountPairVec_1[tmp]);
    			tmpValidDiscrimitiveClassKmerCount_total += tmpClassCount;
	    		if(tmpClassCount > tmpValidDiscrimitiveClassKmerCount_best)
	    		{
	    			tmpValidDiscrimitiveClassKmerCount_best = tmpClassCount;
	    			tmpValidDiscrimitiveClassIndex_best = tmpClassIndex;
	    		}
    		}
    	}

    	int tmpValidDiscrimitiveClassKmerCount_secondBest = 0;
    	int tmpValidDiscrimitiveClassIndex_secondBest = -1;
    	for(int tmp = 0; tmp < tmpValidDiscrimitiveClassIndexCountPairVec.size(); tmp++)
    	{
     		int tmpClassCount = tmpValidDiscrimitiveClassIndexCountPairVec[tmp].second;
    		int tmpClassIndex = (int)tmpValidDiscrimitiveClassIndexCountPairVec[tmp].first;
    		if(tmpClassIndex != tmpValidDiscrimitiveClassIndex_best) // not the best valid discrimitive class
    		{
	    		if(tmpClassCount > tmpValidDiscrimitiveClassKmerCount_secondBest)
	    		{
	    			tmpValidDiscrimitiveClassKmerCount_secondBest = tmpClassCount;
	    			tmpValidDiscrimitiveClassIndex_secondBest = tmpClassIndex;
	    		}
    		}   		
    	}

		validRepetitiveClassKmerCount = tmpValidRepetitiveClassKmerCount;
		validDiscrimitiveClassKmerCount_total = tmpValidDiscrimitiveClassKmerCount_total;
    	invalidClassKmerCount = tmpInvalidClassKmerCount;
    	totalKmerCount = tmpTotalKmerCount;
		validDiscrimitiveClassIndex_best = tmpValidDiscrimitiveClassIndex_best;
		validDiscrimitiveClassIndex_secondBest = tmpValidDiscrimitiveClassIndex_secondBest;
		validDiscrimitiveClassKmerCount_best = tmpValidDiscrimitiveClassKmerCount_best;
		validDiscrimitiveClassKmerCount_secondBest = tmpValidDiscrimitiveClassKmerCount_secondBest;
	}

	void getMostKmerCountRegion_PE_1st_2nd(vector< pair<valueT,int> >& tmpQuerySeq_classIndexCountPairVec_1,
		vector< pair<valueT,int> >& tmpQuerySeq_classIndexCountPairVec_2)
	{
		// merge tmpQuerySeq_classIdCountPairVec_1 and tmpQuerySeq_classIdCountPairVec_2
		// to tmpQuerySeq_classIdCountPairVec_PE
		vector< pair<valueT,int> > tmpQuerySeq_classIndexCountPairVec_PE;
		for(int tmp = 0; tmp < tmpQuerySeq_classIndexCountPairVec_1.size(); tmp++)
			tmpQuerySeq_classIndexCountPairVec_PE.push_back(tmpQuerySeq_classIndexCountPairVec_1[tmp]);
		for(int tmp = 0; tmp < tmpQuerySeq_classIndexCountPairVec_2.size(); tmp++)
		{
			valueT tmpValueT_new = tmpQuerySeq_classIndexCountPairVec_2[tmp].first;
			int tmpCount_new = tmpQuerySeq_classIndexCountPairVec_2[tmp].second;
			int current_tmpQuerySeq_classIndexCountPairVec_PE_size = tmpQuerySeq_classIndexCountPairVec_PE.size();
			bool tmpValueT_exist_bool = false;
			for(int tmp2 = 0; tmp2 < current_tmpQuerySeq_classIndexCountPairVec_PE_size; tmp2++)
			{
				valueT tmpValueT_existing = tmpQuerySeq_classIndexCountPairVec_PE[tmp2].first;
				if(tmpValueT_existing == tmpValueT_new)
				{
					tmpValueT_exist_bool = true;
					(tmpQuerySeq_classIndexCountPairVec_PE[tmp2].second) += tmpCount_new;
					break;
				}
			}
			if(!tmpValueT_exist_bool)
				tmpQuerySeq_classIndexCountPairVec_PE.push_back(pair<valueT,int>(tmpValueT_new, tmpCount_new));
		}
		this->getMostKmerCountRegion_SE_1st_2nd(tmpQuerySeq_classIndexCountPairVec_PE);
	}

	void get_confidence_score_best()
	{
		int validTotalCount_best = validDiscrimitiveClassKmerCount_best + validRepetitiveClassKmerCount;
		confidence_score_best = (double)validTotalCount_best/(double)totalKmerCount; 
	}

	void get_uniqueness_score_best()
	{
		int validTotalCount_best = validDiscrimitiveClassKmerCount_best + validRepetitiveClassKmerCount;
		int validTotalCount_valid = validDiscrimitiveClassKmerCount_total + validRepetitiveClassKmerCount;
		uniqueness_score_best = (double)validTotalCount_best/(double)validTotalCount_valid;
	}	
};
#endif