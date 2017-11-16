// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
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
#include "general/fusionDetectionEvaluation_info.h"

//typedef map<string, int> GeneId2indexMap;

using namespace std;

void parseStr2fieldVec_tab(vector<string>& tmpFieldVec, string& tmpStr)
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

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 geneId_1" << endl;
		cout << "#2 geneId_2" << endl;
		cout << "#2 fusionReadDetectionResults" << endl;
		cout << "#3 specificFusionReadFile" << endl;
		exit(1);
	}
	int pairReadLength = 200;

	string geneId_1 = argv[1];
	string geneId_2 = argv[2];
	string fusionReadDetectionResults = argv[3];
	string specificFusionReadFile = argv[4];

	ifstream fusionRead_ifs(fusionReadDetectionResults.c_str());
	ofstream specificFusionRead_ofs(specificFusionReadFile.c_str());
	while(!fusionRead_ifs.eof())
	{
		string tmpStr;
		getline(fusionRead_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpFieldVec;
		parseStr2fieldVec_tab(tmpFieldVec, tmpStr);
		string tmp_gene_id_1 = tmpFieldVec[1];
		string tmp_gene_id_2 = tmpFieldVec[2];
		if((tmp_gene_id_1 == geneId_1)&&(tmp_gene_id_2 == geneId_2))
			specificFusionRead_ofs << tmpStr << endl;
		else if((tmp_gene_id_1 == geneId_2)&&(tmp_gene_id_2 == geneId_1))
			specificFusionRead_ofs << tmpFieldVec[0] << "\t" 
				<< tmpFieldVec[2] << "\t" << tmpFieldVec[1] << "\t"
				<< tmpFieldVec[4] << "\t" << tmpFieldVec[3] << "\t"
				<< (pairReadLength - atoi(tmpFieldVec[5].c_str())) << "\t" << tmpFieldVec[6] << "\t"
				<< tmpFieldVec[8] << "\t" << tmpFieldVec[7] << "\t"
				<< (pairReadLength - atoi(tmpFieldVec[9].c_str())) << "\t" << tmpFieldVec[10] << endl;
		else
		{}
	}
	fusionRead_ifs.close();
	specificFusionRead_ofs.close();
	return 0;
}