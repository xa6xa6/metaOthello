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

int return_read_taxo_id(string& tmpReadId, int taxo_rank)
{
	string tmpTrimmedReadId;
	if((tmpReadId.at(0) == '>')||(tmpReadId.at(0) == '@'))
		tmpTrimmedReadId = tmpReadId.substr(1);
	else
		tmpTrimmedReadId = tmpReadId;
	vector<string> tmpFieldVec_line;	
	parseStr2fieldVec_line(tmpFieldVec_line, tmpTrimmedReadId);
	string tmpTaxoIdStr = tmpFieldVec_line[8-taxo_rank];
	int tmpId = atoi(tmpTaxoIdStr.c_str());
	if(tmpId < 0)
		return -1;
	else
		return tmpId;
}

void return_assigned_taxo_id(string& tmpResultsStr, int taxo_rank, int& tmpReadId, int& tmpAssignmentId)
{
	vector<string> tmpFieldVec_tab;
	parseStr2fieldVec_tab(tmpFieldVec_tab, tmpResultsStr);
	string tmpReadIdStr = tmpFieldVec_tab[0];
	tmpReadId = return_read_taxo_id(tmpReadIdStr, taxo_rank);
	string tmpTaxoIdStr = tmpFieldVec_tab[8-taxo_rank + 1];
	int tmpId = atoi(tmpTaxoIdStr.c_str());
	if(tmpId < 0)
		tmpAssignmentId = -1;
	else
		tmpAssignmentId = tmpId;
}


int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "#0 Executable " << endl;
		cout << "#1 inputResults_1" << endl;
		cout << "#2 inputResults_2" << endl;
		cout << "#3 Taxo_level_name" << endl;
		cout << "#4 outputDir" << endl;
		exit(1);
	}
	string inputResults_1 = argv[1];
	string inputResults_2 = argv[2];
	string Taxo_name = argv[3];
	string outputDir = argv[4];
	int taxo_rank;
	if((Taxo_name == "Species")||(Taxo_name == "species"))
		taxo_rank = 8;
	else if((Taxo_name == "Genus")||(Taxo_name == "genus"))
		taxo_rank = 7;
	else if((Taxo_name == "Phylum")||(Taxo_name == "phylum"))
		taxo_rank = 3;
	else
	{
		cout << "tmp Taxo_name is not supported!" << endl << "Taxo_name: " << Taxo_name << endl;
		exit(1);
	}

	outputDir += "/";
	string cmd_mkdir = "mkdir " + outputDir;
	system(cmd_mkdir.c_str());
	string correct_incorrect_file = outputDir + "1_correct_2_incorrect.txt";
	string incorrect_correct_file = outputDir + "1_incorrect_2_correct.txt";
	string unmap_incorrect_file = outputDir + "1_unmap_2_incorrect.txt";
	string unmap_correct_file = outputDir + "1_unmap_2_correct.txt";
	string correct_unmap_file = outputDir + "1_correct_2_unmap.txt";
	string incorrect_unmap_file = outputDir + "1_incorrect_2_unmap.txt";
	//string both_correct_file = outputDir + "both_correct.txt";
	//string both_incorrect_file = outputDir + "both_incorrect.txt";
	ofstream correct_incorrect_ofs(correct_incorrect_file.c_str());
	ofstream incorrect_correct_ofs(incorrect_correct_file.c_str());
	ofstream unmap_incorrect_ofs(unmap_incorrect_file.c_str());
	ofstream unmap_correct_ofs(unmap_correct_file.c_str());
	ofstream correct_unmap_ofs(correct_unmap_file.c_str());
	ofstream incorrect_unmap_ofs(incorrect_unmap_file.c_str());	
	//ofstream correct_incorrect_ofs(correct_incorrect_file.c_str());
	//ofstream correct_incorrect_ofs(correct_incorrect_file.c_str());
	ifstream results_1_ifs(inputResults_1.c_str());
	ifstream results_2_ifs(inputResults_2.c_str());
	while((!results_1_ifs.eof())&&(!results_2_ifs.eof()))
	{
		string tmpStr_1, tmpStr_2;
		getline(results_1_ifs, tmpStr_1);
		getline(results_2_ifs, tmpStr_2);
		if((tmpStr_1 == "")||(tmpStr_2 == ""))
			break;

		int tmpReadId_1, tmpReadId_2, tmpAssignmentId_1, tmpAssignmentId_2;
		return_assigned_taxo_id(tmpStr_1, taxo_rank, tmpReadId_1, tmpAssignmentId_1);
		return_assigned_taxo_id(tmpStr_2, taxo_rank, tmpReadId_2, tmpAssignmentId_2);
		if(tmpReadId_1 != tmpReadId_2)
		{
			cout << "error! tmpReadId_1 != tmpReadId_2" << endl;
			exit(1);
		}
		if(tmpReadId_1 < 0)
			continue;
		if((tmpAssignmentId_1 < 0)&&(tmpAssignmentId_2 < 0))
			continue;
		if((tmpAssignmentId_1 == tmpReadId_1)&&(tmpAssignmentId_2 != tmpReadId_1)&&(tmpAssignmentId_2 >= 0))
			correct_incorrect_ofs << tmpStr_1 << "\t" << tmpStr_2 << endl;
		else if((tmpAssignmentId_1 != tmpReadId_1)&&(tmpAssignmentId_2 == tmpReadId_1)&&(tmpAssignmentId_1 >= 0))
			incorrect_correct_ofs << tmpStr_1 << "\t" << tmpStr_2 << endl;
		else if((tmpAssignmentId_1 < 0)&&(tmpAssignmentId_2 != tmpReadId_1))
			unmap_incorrect_ofs << tmpStr_1 << "\t" << tmpStr_2 << endl;
		else if((tmpAssignmentId_1 < 0)&&(tmpAssignmentId_2 == tmpReadId_1))
			unmap_correct_ofs << tmpStr_1 << "\t" << tmpStr_2 << endl;
		else if((tmpAssignmentId_1 != tmpReadId_1)&&(tmpAssignmentId_2 < 0))
			incorrect_unmap_ofs << tmpStr_1 << "\t" << tmpStr_2 << endl;
		else if((tmpAssignmentId_1 == tmpReadId_1)&&(tmpAssignmentId_2 < 0))
			correct_unmap_ofs << tmpStr_1 << "\t" << tmpStr_2 << endl;
		else
		{}
	}
	results_1_ifs.close();
	results_2_ifs.close();
	correct_incorrect_ofs.close();
	incorrect_correct_ofs.close();
	unmap_incorrect_ofs.close();
	unmap_correct_ofs.close();
	correct_unmap_ofs.close();
	incorrect_unmap_ofs.close();		
	return 0;
}