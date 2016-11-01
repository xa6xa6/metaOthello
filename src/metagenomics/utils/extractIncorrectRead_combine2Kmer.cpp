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
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputReadAssignmentFile_Kmer1" << endl;
		cout << "#2 inputReadAssignmentFile_Kmer2" << endl;
		cout << "#3 incorrectlyAssignedReadFile" << endl;
		exit(1);
	}
	string inputReadAssignmentFile_1 = argv[1];
	string inputReadAssignmentFile_2 = argv[2];
	string incorrectlyAssignedReadFile = argv[3];
	ifstream readAssig_ifs_1(inputReadAssignmentFile_1.c_str());
	ifstream readAssig_ifs_2(inputReadAssignmentFile_2.c_str());
	ofstream incorrect_ofs(incorrectlyAssignedReadFile.c_str());
	while(!readAssig_ifs_1.eof())
	{
		string tmpStr_Kmer1;
		getline(readAssig_ifs_1, tmpStr_Kmer1);
		string tmpStr_Kmer2;
		getline(readAssig_ifs_2, tmpStr_Kmer2);		
		if(tmpStr_Kmer1 == "")
			break;
		vector<string> tmpFieldVec_tab_Kmer1;
		vector<string> tmpFieldVec_line_Kmer1;
		parseStr2fieldVec_tab(tmpFieldVec_tab_Kmer1, tmpStr_Kmer1);
		string tmpName_Kmer1 = tmpFieldVec_tab_Kmer1[0];
		string tmpAssignedPhylumId_Kmer1 = tmpFieldVec_tab_Kmer1[1];
		parseStr2fieldVec_line(tmpFieldVec_line_Kmer1, tmpName_Kmer1);
		string tmpTruth_Kmer1 = tmpFieldVec_line_Kmer1[6];
		vector<string> tmpFieldVec_tab_Kmer2;
		vector<string> tmpFieldVec_line_Kmer2;
		parseStr2fieldVec_tab(tmpFieldVec_tab_Kmer2, tmpStr_Kmer2);
		string tmpName_Kmer2 = tmpFieldVec_tab_Kmer2[0];
		string tmpAssignedPhylumId_Kmer2 = tmpFieldVec_tab_Kmer2[1];
		parseStr2fieldVec_line(tmpFieldVec_line_Kmer2, tmpName_Kmer2);
		string tmpTruth_Kmer2 = tmpFieldVec_line_Kmer2[6];		
		if((tmpAssignedPhylumId_Kmer1 != tmpTruth_Kmer1)&&(tmpAssignedPhylumId_Kmer2 != tmpTruth_Kmer2)&&(tmpAssignedPhylumId_Kmer1 == tmpAssignedPhylumId_Kmer2))
			incorrect_ofs << tmpStr_Kmer1 << "\t" << tmpStr_Kmer2 << endl;
	}
	incorrect_ofs.close();
	readAssig_ifs_1.close();
	readAssig_ifs_2.close();
	return 0;
}