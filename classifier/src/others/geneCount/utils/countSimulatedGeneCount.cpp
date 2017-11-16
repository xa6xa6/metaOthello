// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
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

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputGeneIdList" << endl;
		cout << "#2 inputReadFile_singleEnd" << endl;
		cout << "#3 outputGeneCountFile" << endl;
		exit(1);
	}
	string inputGeneIdList = argv[1];
	string inputReadFile_singleEnd = argv[2];
	string outputGeneCountFile = argv[3];
	
	vector<string> geneIdVec;
	ifstream geneId_ifs(inputGeneIdList.c_str());
	while(!geneId_ifs.eof())
	{
		string tmpStr;
		getline(geneId_ifs, tmpStr);
		if(tmpStr == "")
			break;
		geneIdVec.push_back(tmpStr);
	}
	geneId_ifs.close();
	cout << "geneIdVec.size(): " << geneIdVec.size();

	vector<int> geneCountVec;
	for(int tmp = 0; tmp < geneIdVec.size(); tmp++)
		geneCountVec.push_back(0);

	int tmpLineNO = 0;
	ifstream fq_ifs(inputReadFile_singleEnd.c_str());
	while(!fq_ifs.eof())
	{
		string tmpId;
		getline(fq_ifs, tmpId);
		if(tmpId == "")
			break;
		tmpLineNO ++;
		int tmpThousandIndex = tmpLineNO / 1000000;
		if(tmpLineNO == tmpThousandIndex * 1000000)
			cout << "Processed Line #: " << tmpLineNO << endl;		
		string tmpSeq, tmpComm, tmpQual;
		getline(fq_ifs, tmpSeq);
		getline(fq_ifs, tmpComm);
		getline(fq_ifs, tmpQual);
		int underLineLoc = tmpId.find("_");
		int tmpGeneIndex = atoi((tmpId.substr(1, underLineLoc - 1)).c_str()) - 1;
		geneCountVec[tmpGeneIndex] ++;
	}
	fq_ifs.close();

	ofstream geneCountFile_ofs(outputGeneCountFile.c_str());
	for(int tmp = 0; tmp < geneIdVec.size(); tmp++)
		geneCountFile_ofs << tmp + 1 << "\t" << geneIdVec[tmp] << "\t" << geneCountVec[tmp] << endl; 
	geneCountFile_ofs.close();
	return 0;
}