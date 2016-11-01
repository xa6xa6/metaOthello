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
		cout << "#1 inputSpeciesIdCountInfoFile" << endl;
		cout << "#2 inputReadAssignmentFile" << endl;
		cout << "#3 outputFile" << endl;
		exit(1);
	}

	string inputSpeciesIdCountInfoFile = argv[1];
	string inputReadAssignmentFile = argv[2];
	string outputFile = argv[3];

	ifstream speciesIdCount_ifs(inputSpeciesIdCountInfoFile.c_str());

	vector<int> speciesIdVec;
	vector<int> speciesCountVec;
	string tmp1stLine;
	getline(speciesIdCount_ifs, tmp1stLine);
	while(!speciesIdCount_ifs.eof())
	{
		string tmpStr;
		getline(speciesIdCount_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpFieldVec;
		parseStr2fieldVec_tab(tmpFieldVec, tmpStr);
		string speciesIdStr = tmpFieldVec[0];
		int speciesId = atoi(speciesIdStr.c_str());
		string speciesCountStr = tmpFieldVec[6];
		int speciesCount = atoi(speciesCountStr.c_str());
		speciesIdVec.push_back(speciesId);
		speciesCountVec.push_back(speciesCount);
	}

	int speciesNum = speciesIdVec.size();
	cout << "speciesNum: " << speciesNum << endl;
	map<int, int> speciesId2CorrectCountMap;
	map<int, int> speciesId2IncorrectCountMap;
	map<int, int> speciesId2UnmapCountMap;
	for(int tmp = 0; tmp < speciesNum; tmp++)
	{
		speciesId2CorrectCountMap.insert(pair<int,int>(speciesIdVec[tmp], 0));
		speciesId2IncorrectCountMap.insert(pair<int,int>(speciesIdVec[tmp], 0));
		speciesId2UnmapCountMap.insert(pair<int,int>(speciesIdVec[tmp], 0));
	}

	int total_num = 0;
	int unmapped_num = 0;
	unsigned long long tmpLineNO = 0;
	ifstream readAssignment_ifs(inputReadAssignmentFile.c_str());
	while(!readAssignment_ifs.eof())
	{
		string tmpStr;
		getline(readAssignment_ifs, tmpStr);
		if(tmpStr == "")
			break;
        tmpLineNO ++;
        unsigned long long tmpThousandIndex = tmpLineNO / 5000000;
        if(tmpLineNO == tmpThousandIndex * 5000000)
            cout << "Processed Line #: " << tmpLineNO << endl;

        total_num ++;
        vector<string> tmpFieldVec_tab;
        parseStr2fieldVec_tab(tmpFieldVec_tab, tmpStr);

        string tmpAssignmentStr = tmpFieldVec_tab[1];
        int tmpAssignmentId = atoi(tmpAssignmentStr.c_str());

        string tmpReadIdStr = tmpFieldVec_tab[0];
        string tmpReadIdStr_trimmed;
       	if((tmpReadIdStr.at(0) == '>')||(tmpReadIdStr.at(0) == '@'))
        	tmpReadIdStr_trimmed = tmpReadIdStr.substr(1);
        else
        	tmpReadIdStr_trimmed = tmpReadIdStr;
        vector<string> tmpFieldVec_line;
        parseStr2fieldVec_line(tmpFieldVec_line, tmpReadIdStr_trimmed);
        string tmpReadTrueTaxoIdStr = tmpFieldVec_line[0];
        int tmpReadTrueTaxoId = atoi(tmpReadTrueTaxoIdStr.c_str());

        if(tmpReadTrueTaxoId < 0)
        {
        	cout << "invalid tmpReadTrueTaxoIdStr: " << tmpReadTrueTaxoIdStr << endl;
        	exit(1);
        }
        else
        {
        	if(tmpAssignmentId < 0)
        	{
        		map<int, int>::iterator tmpIter = speciesId2UnmapCountMap.find(tmpReadTrueTaxoId);
        		if(tmpIter == speciesId2UnmapCountMap.end())
        		{
        		 	cout << "invalid tmpReadTrueTaxoIdStr in speciesId2UnmapCountMap: " << tmpReadTrueTaxoIdStr << endl;
        			exit(1);       			
        		}
        		else
        			(tmpIter->second) ++;
        	}
        	else if(tmpReadTrueTaxoId == tmpAssignmentId)
        	{
        		map<int, int>::iterator tmpIter = speciesId2CorrectCountMap.find(tmpReadTrueTaxoId);
        		if(tmpIter == speciesId2CorrectCountMap.end())
        		{
        		 	cout << "invalid tmpReadTrueTaxoIdStr in speciesId2CorrectCountMap: " << tmpReadTrueTaxoIdStr << endl;
        			exit(1);       			
        		}
        		else
        			(tmpIter->second) ++;
        	}
        	else
        	{
        		map<int, int>::iterator tmpIter = speciesId2IncorrectCountMap.find(tmpReadTrueTaxoId);
        		if(tmpIter == speciesId2IncorrectCountMap.end())
        		{
        		 	cout << "invalid tmpReadTrueTaxoIdStr in speciesId2IncorrectCountMap: " << tmpReadTrueTaxoIdStr << endl;
        			exit(1);       			
        		}
        		else
        			(tmpIter->second) ++;        		
        	}
        }
	}

	ofstream summary_ofs(outputFile.c_str());	
	for(int tmp = 0; tmp < speciesIdVec.size(); tmp++)
	{
		int tmpSpeciesId = speciesIdVec[tmp];
		int tmpSpeciesCount = speciesCountVec[tmp];
		int tmpSpeciesCount_correct, tmpSpeciesCount_incorrect, tmpSpeciesCount_unmapped;
		map<int, int>::iterator tmpIter_correct = speciesId2CorrectCountMap.find(tmpSpeciesId);
		if(tmpIter_correct == speciesId2CorrectCountMap.end())
		{
			cout << "error!, tmpIter_correct == speciesId2CorrectCountMap.end()" << endl;
			exit(1);
		}
		tmpSpeciesCount_correct = tmpIter_correct->second;
		map<int, int>::iterator tmpIter_incorrect = speciesId2IncorrectCountMap.find(tmpSpeciesId);
		if(tmpIter_incorrect == speciesId2IncorrectCountMap.end())
		{
			cout << "error!, tmpIter_incorrect == speciesId2IncorrectCountMap.end()" << endl;
			exit(1);
		}
		tmpSpeciesCount_incorrect = tmpIter_incorrect->second;		
		map<int, int>::iterator tmpIter_unmapped = speciesId2UnmapCountMap.find(tmpSpeciesId);
		if(tmpIter_unmapped == speciesId2UnmapCountMap.end())
		{
			cout << "error!, tmpIter_unmapped == speciesId2UnmapCountMap.end()" << endl;
			exit(1);
		}
		tmpSpeciesCount_unmapped = tmpIter_unmapped->second;
		summary_ofs << tmpSpeciesId << "\t" << tmpSpeciesCount << "\t" 
			<< tmpSpeciesCount_correct << "\t" << (double)tmpSpeciesCount_correct/(double)tmpSpeciesCount << "\t"
			<< tmpSpeciesCount_incorrect << "\t" << (double)tmpSpeciesCount_incorrect/(double)tmpSpeciesCount << "\t"
			<< tmpSpeciesCount_unmapped << "\t" << (double)tmpSpeciesCount_unmapped/(double)tmpSpeciesCount << endl;
	}
	summary_ofs.close();
	return 0;
}	