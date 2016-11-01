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

typedef map<string, string> NCid2taxoInfoMap; 

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

string get_NCid_from_rawReadId(string& tmpReadId_raw)
{
	int ref_loc = tmpReadId_raw.find("|ref|NC_");
	int dot_loc = tmpReadId_raw.find(".", ref_loc + 1);
	return tmpReadId_raw.substr(ref_loc + 5, dot_loc - ref_loc - 4 - 1);
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 NCid2taxoInfo_file" << endl;
		cout << "#2 inputFa" << endl;
		cout << "#3 outputFile_prefix" << endl;
 		exit(1);
 	}
 	string NCid2taxoInfo_file = argv[1];
 	string inputFa = argv[2];
 	string outputFile_prefix = argv[3];

 	string outputFile_invalid = outputFile_prefix + ".invalid.fa";
 	string outputFile_fa = outputFile_prefix + ".fa";
 	cout << "start to load NCid2taxoInfo_file" << endl;
 	NCid2taxoInfoMap ncId2taxoMap;
 	ifstream NCid2taxoInfo_ifs(NCid2taxoInfo_file.c_str());
 	while(!NCid2taxoInfo_ifs.eof())
 	{
 		string tmpStr;
 		getline(NCid2taxoInfo_ifs, tmpStr);
 		if(tmpStr == "")
 			break;
 		vector<string> tmpFieldVec;
 		parseStr2fieldVec_tab(tmpFieldVec, tmpStr);
 		string tmpNCid = tmpFieldVec[0];
 		string tmpSpeciesId = tmpFieldVec[1];
 		string tmpGenusId = tmpFieldVec[2];
 		string tmpFamilyId = tmpFieldVec[3];
 		string tmpOrderId = tmpFieldVec[4];
 		string tmpClassId = tmpFieldVec[5];
 		string tmpPhylumId = tmpFieldVec[6];
  		if(tmpSpeciesId == "UNKNOWN")
 			tmpSpeciesId = "-1";	
  		if(tmpGenusId == "UNKNOWN")
 			tmpGenusId = "-1";			
  		if(tmpFamilyId == "UNKNOWN")
 			tmpFamilyId = "-1";	
  		if(tmpOrderId == "UNKNOWN")
 			tmpOrderId = "-1";	
  		if(tmpClassId == "UNKNOWN")
 			tmpClassId = "-1";	
  		if(tmpPhylumId == "UNKNOWN")
 			tmpPhylumId = "-1";
 		string tmpMapKeyStr = tmpNCid;
 		string tmpMapValueStr = tmpNCid + "_" + tmpSpeciesId + "_" + tmpGenusId + "_" 
 			+ tmpFamilyId + "_" + tmpOrderId + "_" + tmpClassId + "_" + tmpPhylumId;
 		ncId2taxoMap.insert(pair<string,string>(tmpMapKeyStr, tmpMapValueStr));
 	}
 	NCid2taxoInfo_ifs.close();

 	ifstream fa_ifs(inputFa.c_str());
 	ofstream valid_ofs(outputFile_fa.c_str());
 	ofstream invalid_ofs(outputFile_invalid.c_str());
 	unsigned long long readNum = 0;
 	while(!fa_ifs.eof())
 	{
 		string tmpReadId_raw, tmpSeq;
 		getline(fa_ifs, tmpReadId_raw);
 		if(tmpReadId_raw == "")
 			break;
 		readNum ++;
 		getline(fa_ifs, tmpSeq);
 		string tmpNCid = get_NCid_from_rawReadId(tmpReadId_raw);
 		NCid2taxoInfoMap::iterator tmpIter = ncId2taxoMap.find(tmpNCid);
 		if(tmpIter == ncId2taxoMap.end())
 			invalid_ofs << tmpReadId_raw << endl << tmpSeq << endl;
 		else
 		{
 			string tmpMapValueStr = tmpIter->second;
 			vector<string> tmpLineFieldVec;
 			parseStr2fieldVec_line(tmpLineFieldVec, tmpMapValueStr);
 			string tmpNCid_inMapValue = tmpLineFieldVec[0] + "_" + tmpLineFieldVec[1];
 			string tmpTaxoInfo_inMapValue = tmpLineFieldVec[2] + "_" + tmpLineFieldVec[3] + "_" 
 				+ tmpLineFieldVec[4] + "_" + tmpLineFieldVec[5] + "_" + tmpLineFieldVec[6] + "_"
 				+ tmpLineFieldVec[7];
 			if(tmpNCid_inMapValue != tmpNCid)
 			{
 				cout << "error in searching for NCid in Map!\n query Ncid: " << tmpNCid << endl
 					<< "Ncid in valud: " << tmpNCid_inMapValue << endl;
 				exit(1);
 			}
 			else
 				valid_ofs << ">" << tmpTaxoInfo_inMapValue << "_" << readNum << endl << tmpSeq << endl;
 		}
 	}
 	fa_ifs.close();
 	valid_ofs.close();
 	invalid_ofs.close();
	return 0;
}