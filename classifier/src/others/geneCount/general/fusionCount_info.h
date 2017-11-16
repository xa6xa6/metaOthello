// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#ifndef FUSIONCOUNT_INFO_H
#define FUSIONCOUNT_IFNO_H
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

class FusionCount_Info
{
private:
	vector<string> geneIdVec;
	int fusionReadCount;
	map< int, map<int,int> > fusionGenePairCountMap;
public:
	FusionCount_Info()
	{
		fusionReadCount = 0;
	}

	int return_geneTotalNum()
	{
		return geneIdVec.size();
	}

	string return_geneIdStr(int tmp)
	{
		return geneIdVec[tmp-1];
	}

	void initaite_gene2reissuedFile(string& geneIdListFile)
	{
		ifstream geneInfo_ifs(geneIdListFile.c_str());
		while(!geneInfo_ifs.eof())
		{
			string tmpStr;
			getline(geneInfo_ifs, tmpStr);
			if(tmpStr == "")
				break;
			geneIdVec.push_back(tmpStr);
		}
		geneInfo_ifs.close();
	}

	void addFusionCount(int tmpGeneReissuedId_raw_1, int tmpGeneReissuedId_raw_2)
	{
		int tmpGeneReissuedId_1, tmpGeneReissuedId_2;
		if(tmpGeneReissuedId_raw_1 < tmpGeneReissuedId_raw_2)
		{
			tmpGeneReissuedId_1 = tmpGeneReissuedId_raw_1;
			tmpGeneReissuedId_2 = tmpGeneReissuedId_raw_2;
		}
		else
		{
			tmpGeneReissuedId_1 = tmpGeneReissuedId_raw_2;
			tmpGeneReissuedId_2 = tmpGeneReissuedId_raw_1;
		}

		map< int, map<int,int> >::iterator tmpIter_1 = fusionGenePairCountMap.find(tmpGeneReissuedId_1);
		if(tmpIter_1 != fusionGenePairCountMap.end()) // gene 1 found 
		{
			map<int,int>::iterator tmpIter_2 = (tmpIter_1->second).find(tmpGeneReissuedId_2);
			if(tmpIter_2 != (tmpIter_1->second).end()) // gene 2 found, gene pair found
				(tmpIter_2->second)++;
			else // gene 1 found , gene 2 not found
				(tmpIter_1->second).insert(pair<int,int>(tmpGeneReissuedId_2, 1));
		}
		else // gene 1 not found
		{
			map<int,int> tmpGene2countMap;
			tmpGene2countMap.insert(pair<int,int>(tmpGeneReissuedId_2, 1));
			fusionGenePairCountMap.insert(pair< int, map<int,int> >(tmpGeneReissuedId_1, tmpGene2countMap));
		}
		fusionReadCount ++;
	}

	void printFusionCountInfo(string& tmpOutputFile)
	{
		ofstream fusionCount_ofs(tmpOutputFile.c_str());
		fusionCount_ofs << "fusion read #: " << fusionReadCount << endl << endl;

		for(map< int, map<int,int> >::iterator tmpIter_1 = fusionGenePairCountMap.begin();
			tmpIter_1 != fusionGenePairCountMap.end(); tmpIter_1 ++)
		{
			int tmpGeneReissuedId_1 = tmpIter_1->first;
			for(map<int,int>::iterator tmpIter_2 = (tmpIter_1->second).begin(); 
				tmpIter_2 != (tmpIter_1->second).end(); tmpIter_2++)
			{
				int tmpGeneReissuedId_2 = tmpIter_2->first;
				int tmpFusionGenePairCount = tmpIter_2->second;
				string tmpGeneId_1 = geneIdVec[tmpGeneReissuedId_1 - 1];
				string tmpGeneId_2 = geneIdVec[tmpGeneReissuedId_2 - 1];
				fusionCount_ofs << tmpGeneReissuedId_1 << "\t" << tmpGeneReissuedId_2 << "\t" 
					<< tmpGeneId_1 << "\t" << tmpGeneId_2 << "\t" << tmpFusionGenePairCount << endl;
			}
		}
		fusionCount_ofs.close();
	}
};
#endif