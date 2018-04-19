#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

using namespace std;

vector<string> parseStr(string& tmpStr)
{
	vector<string> tmpFieldVec;
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
	return tmpFieldVec;
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 input_path_to_taxonAssignment_file" << endl;
		cout << "#2 input_path_to_taxonList_file" << endl;		
		cout << "#3 output_path_to_taxonReadCount_file" << endl;
		exit(1);
	}
	string input_path_to_taxonAssignment_file = argv[1];
	string input_path_to_taxonList_file = argv[2];
	string output_path_to_taxonReadCount_file = argv[3];

	unordered_map<int, pair<string, int> > taxon2countMap;

	ifstream taxonList_ifs(input_path_to_taxonList_file.c_str());
	string tmpHeader;
	getline(taxonList_ifs, tmpHeader);
	vector<string> headerFieldVec = parseStr(tmpHeader);
	while(!taxonList_ifs.eof())
	{
		string tmpStr;
		getline(taxonList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpFieldVec = parseStr(tmpStr);
		int tmpId = atoi(tmpFieldVec[1].c_str());
		if(tmpId < 0)
			continue;
		string tmpName = tmpFieldVec[2];
		taxon2countMap.insert(pair<int, pair<string,int> >(
			tmpId, pair<string, int>(tmpName, 0)));
	}
	taxonList_ifs.close();
	cout << "loading taxonAssignment " << endl;
	ifstream taxonAssignment_ifs(input_path_to_taxonAssignment_file.c_str());
	while(!taxonAssignment_ifs.eof())
	{
		string tmpStr;
		getline(taxonAssignment_ifs, tmpStr);
		if(tmpStr == "")
			break;
		//cout << "tmpStr: " << endl << tmpStr << endl;
		vector<string> tmpFieldVec = parseStr(tmpStr);
		for(int tmp = 1; tmp < tmpFieldVec.size(); tmp++)
		{
			string tmpIdStr = tmpFieldVec[tmp];
			int tmpId = atoi(tmpIdStr.c_str());
			//cout << "tmpId: " << tmpId << endl;
			if(tmpId < 0)
				continue;
			unordered_map<int, pair<string, int> >::iterator tmpIter 
				= taxon2countMap.find(tmpId);
			if(tmpIter != taxon2countMap.end()) // found
				((tmpIter->second).second)++;
		}
	}
	taxonAssignment_ifs.close();
	cout << "printing results ..." << endl;
	ofstream taxonReadCount_ofs(output_path_to_taxonReadCount_file.c_str());
	taxonReadCount_ofs << headerFieldVec[1] << "\t" << headerFieldVec[2] 
		<< "\tRead_count" << endl;
	for(unordered_map<int, pair<string, int> >::iterator tmpIter 
			= taxon2countMap.begin(); tmpIter != taxon2countMap.end(); tmpIter++)
	{
		taxonReadCount_ofs << tmpIter->first << "\t" << (tmpIter->second).first
			<< "\t" << (tmpIter->second).second << endl;
	}
	taxonReadCount_ofs.close();
	return 0;
}