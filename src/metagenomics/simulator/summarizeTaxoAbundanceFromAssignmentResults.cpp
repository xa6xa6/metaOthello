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

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputTaxoRankInfoFile" << endl;
		cout << "#2 inputReadAssignmentFile" << endl;
		cout << "#3 outputFile" << endl;
		cout << "#4 taxo_name" << endl;
		exit(1);
	}
	string inputTaxoRankInfoFile = argv[1];
	ifstream taxoRankInfo_ifs(inputTaxoRankInfoFile.c_str());
	string tmp1stLine;
	getline(taxoRankInfo_ifs, tmp1stLine);
	vector<int> taxoIdVec;
	vector<string> taxoNameVec;
	while(!taxoRankInfo_ifs.eof())
	{
		string tmpStr;
		getline(taxoRankInfo_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpFieldVec_tab;
		parseStr2fieldVec_tab(tmpFieldVec_tab, tmpStr);
		string tmpTaxoIdStr = tmpFieldVec_tab[1];
		int tmpTaxoId = atoi(tmpTaxoIdStr.c_str());// << endl;
		string tmpTaxoName = tmpFieldVec_tab[2];
		if(tmpTaxoId >= 0)
		{
			taxoIdVec.push_back(tmpTaxoId);
			taxoNameVec.push_back(tmpTaxoName);
		}
	}

	taxoRankInfo_ifs.close();
	string inputReadAssignmentFile = argv[2];
	string outputDir = argv[3];
	int taxo_rank;
	string taxo_level_name = argv[4];
	if(taxo_level_name == "Phylum")
		taxo_rank = 3;
	else if(taxo_level_name == "Class")
		taxo_rank = 4;
	else if(taxo_level_name == "Order")
		taxo_rank = 5;
	else if(taxo_level_name == "Family")
		taxo_rank = 6;
	else if(taxo_level_name == "Genus")
		taxo_rank = 7;
	else if(taxo_level_name == "Species")
		taxo_rank = 8;
	else
	{
		cout << "invalid taxo_level_name: " << taxo_level_name << endl;
		exit(1);
	}

	string summaryFile = argv[3];
	map<int, int> taxoReadCountMap;

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
        
        string tmpAssignmentStr = tmpFieldVec_tab[8 - taxo_rank + 1];
        int tmpAssignmentId = atoi(tmpAssignmentStr.c_str());

        if(tmpAssignmentId < 0)
        	unmapped_num ++;
        else
        {
        	map<int, int>::iterator tmpIter = taxoReadCountMap.find(tmpAssignmentId);
        	if(tmpIter == taxoReadCountMap.end())
        		taxoReadCountMap.insert(pair<int,int>(tmpAssignmentId, 1));
        	else
        		(tmpIter->second) ++;
        }
	}

	//vector<int> taxoIdVec;
	//vector<string> taxoNameVec;
	vector<int> taxoReadCountVec;
	int tmpTaxoVecSize = taxoIdVec.size();
	for(int tmp = 0; tmp < tmpTaxoVecSize; tmp++)
		taxoReadCountVec.push_back(0);
	double unmapped_perc = (double)unmapped_num/(double)total_num;
	for(map<int, int>::iterator tmpIter = taxoReadCountMap.begin(); 
		tmpIter != taxoReadCountMap.end(); tmpIter++)
	{
		int tmpTaxoId = tmpIter->first;
		int tmpTaxoReadCount = tmpIter->second;
		for(int tmp = 0; tmp < taxoIdVec.size(); tmp++)
		{
			if(taxoIdVec[tmp] == tmpTaxoId)
			{	
				taxoReadCountVec[tmp] = tmpTaxoReadCount;
				break;
			}
		}
	}
	readAssignment_ifs.close();
	ofstream summary_ofs(summaryFile.c_str());
	//summary_ofs << taxo_level_name << "_index\t" << taxo_level_name << "_ID\t" 
	//	<< taxo_level_name << "_name\t" << taxo_level_name << "_readCount\t" 
	//	<< taxo_level_name << "_percInMappedRead\t"
	//	<< taxo_level_name << "_percInTotalRead" << endl;
	for(int tmp = 0; tmp < taxoIdVec.size(); tmp++)
		summary_ofs << tmp << "\t" << taxoIdVec[tmp] << "\t"
			<< taxoNameVec[tmp] << "\t" << taxoReadCountVec[tmp] << "\t"
			<< (double)taxoReadCountVec[tmp]/(double)(total_num - unmapped_num) << "\t"
			<< (double)taxoReadCountVec[tmp]/(double)(total_num) << endl;
	
	//summary_ofs << endl << "Total_#:\t" << total_num << endl 
	//	<< "Mapped_#:\t" << total_num - unmapped_num << endl
	//	<< "Unmapped_#:\t" << unmapped_num << endl; 
	summary_ofs.close();
	return 0;
}