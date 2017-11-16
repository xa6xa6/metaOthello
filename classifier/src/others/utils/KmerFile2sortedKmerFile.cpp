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
#include <chrono>
#include <inttypes.h>
#include "../mps3Lib/read_block_test.h"
#include "../mps3Lib/otherFunc.h"
#include "../mps3Lib/index_info.h"
#include "../general/sortKmerBin_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputRawKmerFile outputSortedKmerFile tmpDir" << endl;
		exit(1);
	}
	string inputRawKmerFile = argv[1];
	string outputSortedKmerFile = argv[2];
	string tmpDir = argv[3];
	tmpDir += "/";
	SortKmerBin_Info tmpSortKmerBinInfo;
	tmpSortKmerBinInfo.initiate_mkTmpDir(inputRawKmerFile, outputSortedKmerFile, tmpDir);
	tmpSortKmerBinInfo.sort_Kmer();
	return 0;
}