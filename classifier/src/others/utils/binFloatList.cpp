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
		cout << "#1 inputFloatListFile" << endl;
		cout << "#2 outputBinCount" << endl;
		cout << "#3 binSize" << endl;
		exit(1);
	}
	string inputFloatListFile = argv[1];
	string outputBinCountFile = argv[2];
	string binSizeStr = argv[3];
	double binSize = atof(binSizeStr.c_str());
	int binNum = 1/binSize;
	if(binNum * binSize != 1)
	{
		cout << "invalid binSize!" << endl;
		cout << "binNum * binSize != 1" << endl;
		exit(1);
	}
	vector<int> binCountVec;
	for(int tmp = 0; tmp < binNum; tmp++)
		binCountVec.push_back(0);
	ifstream floatList_ifs(inputFloatListFile.c_str());
	while(!floatList_ifs.eof())
	{
		string tmpStr;
		getline(floatList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		double tmpFloat = atof(tmpStr.c_str());
		int tmpBinIndex = tmpFloat/binSize;
		if(tmpBinIndex == binNum)
			tmpBinIndex = tmpBinIndex - 1;
		binCountVec[tmpBinIndex] ++;
	}
	floatList_ifs.close();

	ofstream binCount_ofs(outputBinCountFile.c_str());
	for(int tmp = binNum - 1; tmp >= 0; tmp--)
		binCount_ofs << tmp+1 << "\t" << binCountVec[tmp] << endl;
	binCount_ofs.close();
	return 0;
}