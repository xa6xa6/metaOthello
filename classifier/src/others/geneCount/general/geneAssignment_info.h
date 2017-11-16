// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#ifndef GENEASSIGNEMNT_INFO_H
#define GENEASSIGNMENT_IFNO_H
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

class GeneAssignment_Info
{
private:
	int discriminativeClassId_min;
	int discriminativeClassId_max;
	int repetitiveClassId;
public:
	GeneAssignment_Info()
	{}

	void initiate(int tmpDiscriminativeClassId_min, int tmpDiscriminativeClassId_max, int tmpRepetitiveClassId)
	{
		discriminativeClassId_min = tmpDiscriminativeClassId_min;
		discriminativeClassId_max = tmpDiscriminativeClassId_max;
		repetitiveClassId = tmpRepetitiveClassId;
	}

	int getGeneAssignment(vector<valueT>& tmpKmerValueVec)
	{
		int invalidKmerCount = 0;
		int repetitiveKmerCount = 0;
		int discriminativeKmerCount = 0;
		map<int, int> geneId2countMap;
		for(int tmp = 0; tmp < tmpKmerValueVec.size(); tmp++)
		{
			int tmpKmerValue = (int)tmpKmerValueVec[tmp];
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
		for(map<int, int>::iterator tmpIter = geneId2countMap.begin(); tmpIter != geneId2countMap.end(); tmpIter ++)
		{
			int tmpGene_reissuedId = tmpIter->first;
			int tmpGene_count = tmpIter->second;
			if(tmpGene_count > gene_count_best)
			{
				gene_reissuedId_secondBest = gene_reissuedId_best;
				gene_count_secondBest = gene_count_best;
				gene_reissuedId_best = tmpGene_reissuedId;
				gene_count_best = tmpGene_count;
			}	
			else if(tmpGene_count > gene_count_secondBest)
			{
				gene_reissuedId_secondBest = tmpGene_reissuedId;
				gene_count_secondBest = tmpGene_count;
			}
			else
			{}
		}
		if(gene_count_best == gene_count_secondBest)
			return -1;
		else
			return gene_reissuedId_best;
	}
};
#endif