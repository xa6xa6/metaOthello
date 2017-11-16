// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#ifndef MSKW_INFO_H
#define MSKW_INFO_H
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
//#include "taxoClassAssignment_info.h"
using namespace std;
// <taxoReissuedId, vector< pair< 1stKmerIndexInVec, KmerCount > > >
typedef map<int, vector< pair<int,int> > > TaxoReissuedId_2_Window1stKmerVecIndexKmerCountPairVec_Map;
typedef TaxoReissuedId_2_Window1stKmerVecIndexKmerCountPairVec_Map::iterator Taxo2WindowVecMapIter;

class MaximalSpecificKmerWindow_Info
{
private:
	TaxoReissuedId_2_Window1stKmerVecIndexKmerCountPairVec_Map taxoId2windwoVecMap;

public:
	MaximalSpecificKmerWindow_Info()
	{

	}

	void get_taxoReissuedId_maximalWindowSquareSum_best(int& tmpMaximalWindowKmerCountSquareSum_taxoReissuedId, 
		int& tmpMaximalWindowKmerCountSquareSum, int assignment_KmerCountSquareSum_min)
	{
		tmpMaximalWindowKmerCountSquareSum_taxoReissuedId = -2;
		tmpMaximalWindowKmerCountSquareSum = 0;
		for(Taxo2WindowVecMapIter tmpIter = taxoId2windwoVecMap.begin(); tmpIter != taxoId2windwoVecMap.end(); tmpIter ++)
		{
			int tmpReissuedId = tmpIter->first;
			int tmpReissuedIdWindowVecSize = (tmpIter->second).size();
			int tmpSquareSum = 0;
			for(int tmp = 0; tmp < tmpReissuedIdWindowVecSize; tmp++)
			{
				int tmpKmerCount = (tmpIter->second)[tmp].second;
				tmpSquareSum += (tmpKmerCount * tmpKmerCount);
			}	
			if(tmpSquareSum > tmpMaximalWindowKmerCountSquareSum)
			{
				tmpMaximalWindowKmerCountSquareSum_taxoReissuedId = tmpReissuedId;
				tmpMaximalWindowKmerCountSquareSum = tmpSquareSum;
			}
		}
	}

	void get_taxoReissuedId_maximalWindow_best(int& tmpMaximalWindowKmerCount_taxoReissuedId,
		int& tmpMaximalWindowKmerCount, int assignment_KmerCount_min)
	{
		tmpMaximalWindowKmerCount_taxoReissuedId = -2;
		tmpMaximalWindowKmerCount = 0;
		for(Taxo2WindowVecMapIter tmpIter = taxoId2windwoVecMap.begin(); tmpIter != taxoId2windwoVecMap.end(); tmpIter ++)
		{
			int tmpReissuedId = tmpIter->first;
			int tmpReissuedIdWindowVecSize = (tmpIter->second).size();
			for(int tmp = 0; tmp < tmpReissuedIdWindowVecSize; tmp++)
			{
				int tmpKmerCount = (tmpIter->second)[tmp].second;
				if(tmpKmerCount > tmpMaximalWindowKmerCount)
				{	
					tmpMaximalWindowKmerCount_taxoReissuedId = tmpReissuedId;
					tmpMaximalWindowKmerCount = tmpKmerCount;
				}
			}	
		}
	}

	void get_taxoReissuedId_maximalWindowSquareSum_best_secondBest(
		int& tmpMaximalWindowKmerCountSquareSum_taxoReissuedId_best, 
		int& tmpMaximalWindowKmerCountSquareSum_best, 
		int& tmpMaximalWindowKmerCountSquareSum_taxoReissuedId_secondBest, 
		int& tmpMaximalWindowKmerCountSquareSum_secondBest, 
		double assignment_windowSizeSquareSum_min)
	{
		tmpMaximalWindowKmerCountSquareSum_taxoReissuedId_best = -2;
		tmpMaximalWindowKmerCountSquareSum_best = 0;
		tmpMaximalWindowKmerCountSquareSum_taxoReissuedId_secondBest = -2;
		tmpMaximalWindowKmerCountSquareSum_secondBest = 0;		
		for(Taxo2WindowVecMapIter tmpIter = taxoId2windwoVecMap.begin(); tmpIter != taxoId2windwoVecMap.end(); tmpIter ++)
		{
			int tmpReissuedId = tmpIter->first;
			int tmpReissuedIdWindowVecSize = (tmpIter->second).size();
			int tmpSquareSum = 0;
			for(int tmp = 0; tmp < tmpReissuedIdWindowVecSize; tmp++)
			{
				int tmpKmerCount = (tmpIter->second)[tmp].second;
				tmpSquareSum += (tmpKmerCount * tmpKmerCount);
			}	
			if((tmpSquareSum > tmpMaximalWindowKmerCountSquareSum_best)&&(tmpSquareSum >= assignment_windowSizeSquareSum_min))
			{
				tmpMaximalWindowKmerCountSquareSum_taxoReissuedId_secondBest = tmpMaximalWindowKmerCountSquareSum_taxoReissuedId_best;
				tmpMaximalWindowKmerCountSquareSum_secondBest = tmpMaximalWindowKmerCountSquareSum_best;
				tmpMaximalWindowKmerCountSquareSum_taxoReissuedId_best = tmpReissuedId;
				tmpMaximalWindowKmerCountSquareSum_best = tmpSquareSum;
			}
			else if((tmpSquareSum > tmpMaximalWindowKmerCountSquareSum_secondBest)&&(tmpSquareSum >= assignment_windowSizeSquareSum_min))
			{
				tmpMaximalWindowKmerCountSquareSum_taxoReissuedId_secondBest = tmpReissuedId;
				tmpMaximalWindowKmerCountSquareSum_secondBest = tmpSquareSum;
			}
			else
			{}
		}
	}

	void get_taxoReissuedId_maximalWindow_best_secondBest(
		int& tmpMaximalWindowKmerCount_taxoReissuedId_best,
		int& tmpMaximalWindowKmerCount_best, 		
		int& tmpMaximalWindowKmerCount_taxoReissuedId_secondBest,
		int& tmpMaximalWindowKmerCount_secondBest, 
		int assignment_windowSize_min)
	{
		tmpMaximalWindowKmerCount_taxoReissuedId_best = -2;
		tmpMaximalWindowKmerCount_best = 0;
		tmpMaximalWindowKmerCount_taxoReissuedId_secondBest = -2;
		tmpMaximalWindowKmerCount_secondBest = 0;		
		for(Taxo2WindowVecMapIter tmpIter = taxoId2windwoVecMap.begin(); tmpIter != taxoId2windwoVecMap.end(); tmpIter ++)
		{
			int tmpReissuedId = tmpIter->first;
			//cout << "tmpReissuedId: " << tmpReissuedId << endl;
			int tmpReissuedIdWindowVecSize = (tmpIter->second).size();
			//cout << "tmpReissuedIdWindowVecSize: " << tmpReissuedIdWindowVecSize << endl;
			for(int tmp = 0; tmp < tmpReissuedIdWindowVecSize; tmp++)
			{
				//cout << "tmp1stKmerIndex: " << (tmpIter->second)[tmp].first << endl;
				int tmpKmerCount = (tmpIter->second)[tmp].second;
				//cout << "tmpKmerCount: " << tmpKmerCount << endl;
				if((tmpKmerCount > tmpMaximalWindowKmerCount_best)&&(tmpKmerCount >= assignment_windowSize_min))
				{		
					tmpMaximalWindowKmerCount_taxoReissuedId_secondBest = tmpMaximalWindowKmerCount_taxoReissuedId_best;
					tmpMaximalWindowKmerCount_secondBest = tmpMaximalWindowKmerCount_best;							
					tmpMaximalWindowKmerCount_taxoReissuedId_best = tmpReissuedId;
					tmpMaximalWindowKmerCount_best = tmpKmerCount;
				}
				else if((tmpKmerCount > tmpMaximalWindowKmerCount_secondBest)&&(tmpKmerCount >= assignment_windowSize_min))
				{
					tmpMaximalWindowKmerCount_taxoReissuedId_secondBest = tmpReissuedId;
					tmpMaximalWindowKmerCount_secondBest = tmpKmerCount;						
				}
				else
				{}
			}	
		}
		//cout << "tmpMaximalWindowKmerCount_taxoReissuedId_best: " << tmpMaximalWindowKmerCount_taxoReissuedId_best << endl;
		//cout << "tmpMaximalWindowKmerCount_best: " << tmpMaximalWindowKmerCount_best << endl;
		//cout << "tmpMaximalWindowKmerCount_taxoReissuedId_secondBest: " << tmpMaximalWindowKmerCount_taxoReissuedId_secondBest << endl;
		//cout << "tmpMaximalWindowKmerCount_secondBest: " << tmpMaximalWindowKmerCount_secondBest << endl;
	}


	string print_taxoId2windowMap_2_string()
	{
		string tmpStr;
		for(Taxo2WindowVecMapIter tmpIter = taxoId2windwoVecMap.begin(); tmpIter != taxoId2windwoVecMap.end(); tmpIter ++)
		{
			int tmpReissuedId = tmpIter->first;
			tmpStr += int_to_str(tmpReissuedId);
			tmpStr += ":";
			int tmpReissuedIdWindowVecSize = (tmpIter->second).size();
			for(int tmp = 0; tmp < tmpReissuedIdWindowVecSize; tmp++)
			{
				int tmp1stKmerIndex = (tmpIter->second)[tmp].first;
				int tmpKmerCount = (tmpIter->second)[tmp].second;
				tmpStr += int_to_str(tmp1stKmerIndex);
				tmpStr += "-";
				tmpStr += int_to_str(tmpKmerCount);
				tmpStr += ";";
			}	
			tmpStr += "\n";
		}		
		return tmpStr;
	}

	void insert_window_taxoReissuedId_1stKmerIndexInVec_KmerCount(int currentWindow_taxoReissuedId, 
		int currentWindow_1stKmerIndexInVec, int currentWindow_KmerCount)
	{
		Taxo2WindowVecMapIter tmpIter = taxoId2windwoVecMap.find(currentWindow_taxoReissuedId);
		if(tmpIter == taxoId2windwoVecMap.end()) // new, insert
		{
			vector< pair<int,int> > tmpPairVec;
			tmpPairVec.push_back(pair<int,int>(currentWindow_1stKmerIndexInVec, currentWindow_KmerCount));
			taxoId2windwoVecMap.insert(pair<int, vector< pair<int,int> > >(currentWindow_taxoReissuedId, tmpPairVec));
		}
		else // existing, update vec
			(tmpIter->second).push_back(pair<int,int>(currentWindow_1stKmerIndexInVec, currentWindow_KmerCount));
	}


};
#endif