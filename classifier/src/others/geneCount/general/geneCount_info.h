// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#ifndef GENECOUNT_INFO_H
#define GENECOUNT_IFNO_H
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

class GeneCount_Info
{
private:
	vector<string> geneIdVec;
	vector<int> geneCountVec;

	int unassignedReadCount;
public:
	GeneCount_Info()
	{
		unassignedReadCount = 0;
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
			geneCountVec.push_back(0);
		}
		geneInfo_ifs.close();
	}

	void addGeneCount(int tmpGeneReissuedId)
	{
		int tmpGeneIndex = tmpGeneReissuedId - 1;
		geneCountVec[tmpGeneIndex] ++;
	}

	void addUnassignedReadCount()
	{
		unassignedReadCount ++;
	}

	void printGeneCountInfo(string& tmpOutputFile)
	{
		ofstream geneCount_ofs(tmpOutputFile.c_str());
		for(int tmp = 0; tmp < geneCountVec.size(); tmp++)
			geneCount_ofs << geneIdVec[tmp] << "\t" << geneCountVec[tmp] << endl;
		//geneCount_ofs << endl << "Unassigned read #:\t" << unassignedReadCount << endl;
		geneCount_ofs.close();
	}

	void printAssignStatsInfo(string& tmpOutputFile)
	{
		ofstream assignStats_ofs(tmpOutputFile.c_str());
		int assignedReadCout_total = 0;
		for(int tmp = 0; tmp < geneCountVec.size(); tmp++)
			assignedReadCout_total += geneCountVec[tmp];
		assignStats_ofs << "Assigned read #:\t" << assignedReadCout_total << endl;
		assignStats_ofs << "Unassigned read #:\t" << unassignedReadCount << endl;
		assignStats_ofs.close();
	}
};
#endif