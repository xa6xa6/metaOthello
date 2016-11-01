#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <bitset>
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

typedef map<string, unsigned long long> EMBLid2taxoIdMap;

void parseStr2fieldVec(vector<string>& tmpFieldVec, string& tmpStr)
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
	if(argc != 10)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputGb2NCBItaxoIdFile" << endl;
		cout << "#2 inputGenomeId2taxoIdFile" << endl;
		cout << "#3 GenomeId_max" << endl;
		cout << "#4 inputFqFile_1_withEMBLid" << endl;
		cout << "#5 inputFqFile_2_withEMBLid" << endl;
		cout << "#6 outputFormattedFqFile_1_withTaxoId" << endl;
		cout << "#7 outputFormattedFqFile_2_withTaxoId" << endl;
		cout << "#8 outputInvalidFile" << endl;
		cout << "#9 outputInvalidFile" << endl;
		exit(1);
	}
	cout << "start to initaite genomeId2taxoId vecs" << endl;
	string inputGb2NCBItaxoIdFile = argv[1];
	string inputGenomeId2taxoIdFile = argv[2];
	string GenomeId_max_str = argv[3];
	string inputFqFile_1 = argv[4];
	string inputFqFile_2 = argv[5];
	string outputFormattedFqFile_1 = argv[6];
	string outputFormattedFqFile_2 = argv[7];
	string outputInvalidFqFile_1 = argv[8];
	string outputInvalidFqFile_2 = argv[9];

	int GenomeId_max = atoi(GenomeId_max_str.c_str());
	vector<int> speciesIdVec;
	vector<int> genusIdVec;
	vector<int> familyIdVec;
	vector<int> orderIdVec;
	vector<int> classIdVec;
	vector<int> phylumIdVec;
	for(int tmp = 0; tmp < GenomeId_max; tmp++)
	{
		speciesIdVec.push_back(0);
		genusIdVec.push_back(0);
		familyIdVec.push_back(0);
		orderIdVec.push_back(0);
		classIdVec.push_back(0);
		phylumIdVec.push_back(0);
	}
	ifstream genomeId2taxoId_ifs(inputGenomeId2taxoIdFile.c_str());
	while(!genomeId2taxoId_ifs.eof())
	{
		string tmpStr;
		getline(genomeId2taxoId_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpFieldVec;
		parseStr2fieldVec(tmpFieldVec, tmpStr);
		int tmpGenomeId = atoi(tmpFieldVec[1].c_str());
		int tmpSpeciesId = atoi(tmpFieldVec[2].c_str());
		int tmpGenusId = atoi(tmpFieldVec[3].c_str());
		int tmpFamilyId = atoi(tmpFieldVec[4].c_str());
		int tmpOrderId = atoi(tmpFieldVec[5].c_str());
		int tmpClassId = atoi(tmpFieldVec[6].c_str());
		int tmpPhylumId = atoi(tmpFieldVec[7].c_str());
		speciesIdVec[tmpGenomeId - 1] = tmpSpeciesId;
		genusIdVec[tmpGenomeId - 1] = tmpGenusId;
		familyIdVec[tmpGenomeId - 1] = tmpFamilyId;
		orderIdVec[tmpGenomeId - 1] = tmpOrderId;
		classIdVec[tmpGenomeId - 1] = tmpClassId;
		phylumIdVec[tmpGenomeId - 1] = tmpPhylumId;
	}
	genomeId2taxoId_ifs.close();
	cout << "end of initiating genomeId2taxoId vecs" << endl;
	cout << "start to initiate Gb2NCBItaxoId map" << endl;
	ifstream gb2NCBItaxoId_ifs(inputGb2NCBItaxoIdFile.c_str());
	string firstLine;
	getline(gb2NCBItaxoId_ifs, firstLine);
	EMBLid2taxoIdMap tmpEMBLid2taxoIdMap;
	unsigned long long tmpLineNO = 0; 
	while(!gb2NCBItaxoId_ifs.eof())
	{
		string tmpStr;
		getline(gb2NCBItaxoId_ifs, tmpStr);
		if(tmpStr == "")
			break;
        tmpLineNO ++;
        unsigned long long tmpThousandIndex = tmpLineNO / 5000000;
        if(tmpLineNO == tmpThousandIndex * 5000000)
            cout << "Processed Line #: " << tmpLineNO << endl; 		
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		string tmpEMBLid = tmpStr.substr(0, tabLoc_1);
		string tmpTaxoIdStr = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		unsigned long long tmpTaxoId = atoll(tmpTaxoIdStr.c_str());
		tmpEMBLid2taxoIdMap.insert(pair<string, unsigned long long>(tmpEMBLid, tmpTaxoId));
	}
	gb2NCBItaxoId_ifs.close();
	cout << "end of initiating Gb2NCBItaxoId map" << endl;

	cout << "start to covert raw fq with embl id files 2 fq with taxo id files" << endl;
	ifstream oriFq_ifs_1(inputFqFile_1.c_str());
	ifstream oriFq_ifs_2(inputFqFile_2.c_str());
	ofstream validFq_ofs_1(outputFormattedFqFile_1.c_str());
	ofstream validFq_ofs_2(outputFormattedFqFile_2.c_str());
	ofstream invalidFq_ofs_1(outputInvalidFqFile_1.c_str());
	ofstream invalidFq_ofs_2(outputInvalidFqFile_2.c_str());
	unsigned long long tmpReadId = 0;
	while(!oriFq_ifs_1.eof())
	{
		string tmpStr_1, tmpStr_2, tmpStr_3, tmpStr_4;
		getline(oriFq_ifs_1, tmpStr_1);
		if(tmpStr_1 == "")
			break;
		tmpReadId ++;
		getline(oriFq_ifs_1, tmpStr_2);
		getline(oriFq_ifs_1, tmpStr_3);
		getline(oriFq_ifs_1, tmpStr_4);
		string tmpStr_5, tmpStr_6, tmpStr_7, tmpStr_8;
		getline(oriFq_ifs_2, tmpStr_5);
		getline(oriFq_ifs_2, tmpStr_6);
		getline(oriFq_ifs_2, tmpStr_7);
		getline(oriFq_ifs_2, tmpStr_8);
		int lineLoc = tmpStr_1.find("-");
		string tmpEMBLid = tmpStr_1.substr(1, lineLoc - 1);
		EMBLid2taxoIdMap::iterator tmpIter = tmpEMBLid2taxoIdMap.find(tmpEMBLid);
		if(tmpIter == tmpEMBLid2taxoIdMap.end())//not found
		{
			invalidFq_ofs_1 << tmpStr_1 << endl;
			invalidFq_ofs_1 << tmpStr_2 << endl;
			invalidFq_ofs_1 << tmpStr_3 << endl;
			invalidFq_ofs_1 << tmpStr_4 << endl;
			invalidFq_ofs_2 << tmpStr_5 << endl;
			invalidFq_ofs_2 << tmpStr_6 << endl;
			invalidFq_ofs_2 << tmpStr_7 << endl;
			invalidFq_ofs_2 << tmpStr_8 << endl;
		}
		else
		{
			int tmpGenomeId = tmpIter->second;
			int tmpSpeciesId;
			if(tmpGenomeId > GenomeId_max)
				tmpSpeciesId = 0;
			else
				tmpSpeciesId = speciesIdVec[tmpGenomeId - 1];
			if(tmpSpeciesId == 0)
			{
				invalidFq_ofs_1 << tmpStr_1 << endl;
				invalidFq_ofs_1 << tmpStr_2 << endl;
				invalidFq_ofs_1 << tmpStr_3 << endl;
				invalidFq_ofs_1 << tmpStr_4 << endl;
				invalidFq_ofs_2 << tmpStr_5 << endl;
				invalidFq_ofs_2 << tmpStr_6 << endl;
				invalidFq_ofs_2 << tmpStr_7 << endl;
				invalidFq_ofs_2 << tmpStr_8 << endl;				
			}
			else
			{
				int tmpGenusId = genusIdVec[tmpGenomeId - 1];
				int tmpFamilyId = familyIdVec[tmpGenomeId - 1];
				int tmpOrderId = orderIdVec[tmpGenomeId - 1];
				int tmpClassId = classIdVec[tmpGenomeId - 1];
				int tmpPhylumId = phylumIdVec[tmpGenomeId - 1];
				string tmpTaxoPathId = int_to_str(tmpGenomeId) + "_" 
					+ int_to_str(tmpSpeciesId) + "_"
					+ int_to_str(tmpGenusId) + "_"
					+ int_to_str(tmpFamilyId) + "_"
					+ int_to_str(tmpOrderId) + "_"
					+ int_to_str(tmpClassId) + "_"
					+ int_to_str(tmpPhylumId);
				validFq_ofs_1 << "@" << tmpTaxoPathId << "_" << tmpReadId << "/1" << endl;
				validFq_ofs_1 << tmpStr_2 << endl;
				validFq_ofs_1 << tmpStr_3 << endl;
				validFq_ofs_1 << tmpStr_4 << endl;		
				validFq_ofs_2 << "@" << tmpTaxoPathId << "_" << tmpReadId << "/2" << endl;
				validFq_ofs_2 << tmpStr_6 << endl;
				validFq_ofs_2 << tmpStr_7 << endl;
				validFq_ofs_2 << tmpStr_8 << endl;						
			}
		}
	}
	cout << "All jobs done!" << endl;

	oriFq_ifs_1.close();
	oriFq_ifs_2.close();
	validFq_ofs_1.close();
	validFq_ofs_2.close();
	invalidFq_ofs_1.close();
	invalidFq_ofs_2.close();
	return 0;
}