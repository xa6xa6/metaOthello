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
#include "../../mps3Lib/read_block_test.h"
#include "../../mps3Lib/otherFunc.h"
#include "../../mps3Lib/index_info.h"
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
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 AssignmentResuls_1" << endl;
		cout << "#2 AssignmentResuls_2" << endl;
		cout << "#3 outputResults" << endl;
		exit(1);
	}
	string AssignmentResuls_1 = argv[1];
	string AssignmentResuls_2 = argv[2];
	string outputResults = argv[3];
	ifstream assign_ifs_1(AssignmentResuls_1.c_str());
	ifstream assign_ifs_2(AssignmentResuls_2.c_str());
	ofstream assign_ofs(outputResults.c_str());
	while((!assign_ifs_1.eof())&&(!assign_ifs_2.eof()))
	{
		string tmpStr_1;
		string tmpStr_2;
		getline(assign_ifs_1, tmpStr_1);
		getline(assign_ifs_2, tmpStr_2);
		if((tmpStr_1 == "")||(tmpStr_2 == ""))
			break;
		vector<string> tmpFieldVec_1;
		vector<string> tmpFieldVec_2;
		parseStr2fieldVec_tab(tmpFieldVec_1, tmpStr_1);
		parseStr2fieldVec_tab(tmpFieldVec_2, tmpStr_2);
		string tmpAssignmentStr_1 = tmpFieldVec_1[1];
		string tmpAssignmentStr_2 = tmpFieldVec_2[1];
		int tmpAssignment_1 = atoi(tmpAssignmentStr_1.c_str());
		int tmpAssignment_2 = atoi(tmpAssignmentStr_2.c_str());
		if(tmpAssignment_1 < 0)
			assign_ofs << tmpStr_2 << endl;
		else
			assign_ofs << tmpStr_1 << endl;
	}
	assign_ifs_1.close();
	assign_ifs_2.close();
	assign_ofs.close();
	return 0;
}