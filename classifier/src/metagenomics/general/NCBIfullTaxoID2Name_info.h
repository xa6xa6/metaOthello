// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#ifndef NCBIFULLTAXOID2NAME_INFO_H
#define NCBIFULLTAXOID2NAME_IFNO_H
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

class NCBIfullTaxoID2Name_Info
{
private:
	vector<string> taxoNameVec;
public:
	NCBIfullTaxoID2Name_Info()
	{}

	int return_NCBIfullTaxoIdMax()
	{
		return taxoNameVec.size() - 1;
	}

	string return_taxoName(int taxoId)
	{
		return taxoNameVec[taxoId];
	}

	void initiate_taxoID2NameFile(string& taxoId2NameFile)
	{
		vector<int> interIdVec;
		vector<string> interNameVec;
		ifstream taxoId2Name_ifs(taxoId2NameFile.c_str());
		int currentMaxId = 0;
		while(!taxoId2Name_ifs.eof())
		{
			string tmpStr;
			getline(taxoId2Name_ifs, tmpStr);
			if(tmpStr == "")
				break;
			vector<string> tmpFieldVec;
			this->parseStr2fieldVec(tmpFieldVec, tmpStr);
			string tmpIdStr = tmpFieldVec[0];
			//cout << "tmpIdStr: " << tmpIdStr << endl;
			int tmpIdInt = atoi(tmpIdStr.c_str());
			if(tmpIdInt > currentMaxId)
				currentMaxId = tmpIdInt;
			interIdVec.push_back(tmpIdInt);
			interNameVec.push_back(tmpFieldVec[2]);
		}
		taxoId2Name_ifs.close();
		cout << "currentMaxId: " << currentMaxId << endl;
		cout << "interIdVec.size(): " << interIdVec.size() << endl;
		for(int tmp = 0; tmp <= currentMaxId; tmp++)
			taxoNameVec.push_back("NULL");
		for(int tmp = 0; tmp < interIdVec.size(); tmp++)
		{
			int tmpId = interIdVec[tmp];
			string tmpName = interNameVec[tmp];
			taxoNameVec[tmpId] = tmpName;
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