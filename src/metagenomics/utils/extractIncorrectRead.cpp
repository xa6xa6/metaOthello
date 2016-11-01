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

void parseStr2fieldVec_line(vector<string>& tmpFieldVec, string& tmpStr)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find("_", startLoc);
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
	if(argc != 3)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputReadAssignmentFile" << endl;
		cout << "#2 incorrectlyAssignedReadFile" << endl;
		exit(1);
	}
	string inputReadAssignmentFile = argv[1];
	string incorrectlyAssignedReadFile = argv[2];
	ifstream readAssig_ifs(inputReadAssignmentFile.c_str());
	ofstream incorrect_ofs(incorrectlyAssignedReadFile.c_str());
	while(!readAssig_ifs.eof())
	{
		string tmpStr;
		getline(readAssig_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpFieldVec_tab;
		vector<string> tmpFieldVec_line;
		parseStr2fieldVec_tab(tmpFieldVec_tab, tmpStr);
		string tmpName = tmpFieldVec_tab[0];
		string tmpAssignedIdStr = tmpFieldVec_tab[1];
		//int tmpAssignedId = atoi(tmpAssignedIdStr.c_str());
		parseStr2fieldVec_line(tmpFieldVec_line, tmpName);
		string tmpTruth = (tmpFieldVec_line[0]).substr(1);
		if(tmpAssignedIdStr != tmpTruth)
			incorrect_ofs << tmpStr << endl;
	}
	incorrect_ofs.close();
	readAssig_ifs.close();
	return 0;
}