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
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <errno.h>

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputSpeciesId2taxoIdFile" << endl;
		cout << "#2 inputFqFile" << endl;
		cout << "#3 outputSpeciesReadCountFile" << endl;
		exit(1);
	}
	cout << "starts ..." << endl;
	string inputSpeciesId2taxoIdFile = argv[1];
	string inputFqFile = argv[2];
	string outputSpeciesReadCountFile = argv[3];
	unsigned long long SpeciesId_max = 0;
	ifstream speciesInfo_ifs(inputSpeciesId2taxoIdFile.c_str());
	vector<unsigned long long> speciesIdVec;
	while(!speciesInfo_ifs.eof())
	{
		string tmpStr;
		getline(speciesInfo_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpSpeciesIdStr = tmpStr.substr(0, tabLoc);
		unsigned long long tmpSpeciesId = atoll(tmpSpeciesIdStr.c_str());
		speciesIdVec.push_back(tmpSpeciesId);
		if(tmpSpeciesId > SpeciesId_max)
			SpeciesId_max = tmpSpeciesId;
	}
	speciesInfo_ifs.close();
	cout << "SpeciesId_max: " << SpeciesId_max << endl;
	vector<int> speciesCountVec;
	for(int tmp = 0; tmp < SpeciesId_max; tmp++)
		speciesCountVec.push_back(0);

	ifstream fq_ifs(inputFqFile.c_str());
	while(!fq_ifs.eof())
	{
		string tmpStr;
		getline(fq_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int lineLoc_1 = tmpStr.find("_");
		int lineLoc_2 = tmpStr.find("_", lineLoc_1 + 1);
		string tmpSpeciesIdStr = tmpStr.substr(lineLoc_1 + 1, lineLoc_2 - lineLoc_1 - 1);
		//cout << "tmpSpeciesIdStr: " << tmpSpeciesIdStr << endl;
		unsigned long long tmpSpeciesId = atoll(tmpSpeciesIdStr.c_str());
		speciesCountVec[tmpSpeciesId - 1] ++;
		string tmpStr_2, tmpStr_3, tmpStr_4;
		getline(fq_ifs, tmpStr_2);
		getline(fq_ifs, tmpStr_3);
		getline(fq_ifs, tmpStr_4);
	}
	fq_ifs.close();
	cout << "start to do count" << endl;
	ofstream count_ofs(outputSpeciesReadCountFile.c_str());
	for(int tmp = 0; tmp < speciesIdVec.size(); tmp++)
	{
		unsigned long long tmpSpeciesId = speciesIdVec[tmp];
		int tmpCount = speciesCountVec[tmpSpeciesId - 1];
		count_ofs << tmpSpeciesId << "\t" << tmpCount << endl;
	}
	count_ofs.close();
	return 0;
}