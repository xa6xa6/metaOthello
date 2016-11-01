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
		cout << "#1 inputKmerSetIdFile" << endl;
		cout << "#2 total_set_num" << endl;
		cout << "#3 KmerSetIdDistributionFile" << endl;
		exit(1);
	}
	string inputKmerSetIdFile = argv[1];
	string total_set_num_str = argv[2];
	string KmerSetIdDistributionFile = argv[3];
	int total_set_num = atoi(total_set_num_str.c_str());
	unsigned long long* freqArray = new unsigned long long[total_set_num + 1];
	ifstream KmerSetId_ifs(inputKmerSetIdFile.c_str());
	unsigned long long tmpLineNO = 0;
	while(!KmerSetId_ifs.eof())
	{
		string tmpStr;
		getline(KmerSetId_ifs, tmpStr);
		if(tmpStr == "")
			break;
        
        tmpLineNO ++;
        unsigned long long tmpThousandIndex = tmpLineNO / 5000000;
        if(tmpLineNO == tmpThousandIndex * 5000000)
            cout << "Processed Line #: " << tmpLineNO << endl;

		int tmpSetId = atoi(tmpStr.c_str());
		if((tmpSetId < 0)||(tmpSetId > total_set_num))
		{
			cout << "invalid tmpSetId: " << tmpSetId << endl;
			exit(1);
		}
		freqArray[tmpSetId] ++;
	}
	KmerSetId_ifs.close();
	ofstream KmerSetDistribution_ofs(KmerSetIdDistributionFile.c_str());
	for(int tmp = 0; tmp <= total_set_num; tmp++)
		KmerSetDistribution_ofs << "Set\t" << tmp << "\t" << freqArray[tmp] << endl;
	KmerSetDistribution_ofs.close();

	delete freqArray;
	return 0;
}