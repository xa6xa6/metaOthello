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
#include "../../Ye_implementation/othello.h"
#include "../../Ye_implementation/mulothindex.h"
#include <chrono>
#include <inttypes.h>
#include "../../Ye_implementation/io_helper.h"
#include "../../mps3Lib/read_block_test.h"
#include "../../mps3Lib/otherFunc.h"
#include "../../mps3Lib/index_info.h"

using namespace std;
typedef unsigned long long keyT;
typedef uint64_t valueT;

vector<keyT> keys;
vector<valueT> values;

IOHelper<keyT,valueT> *helper;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable splitbit Kmer_length inputKmerTwoColumnFile outputKmerIndexBitArray" << endl;
		exit(1);
	}
	cout << "Kmer 2 index starts ..." << endl;
	int splitbit;
	string splitbitStr = argv[1];
	splitbit = atoi(splitbitStr.c_str());

	string Kmer_length_str = argv[2];
	int Kmer_length = atoi(Kmer_length_str.c_str());

	cout << "splitbit: " << splitbit << endl;
	cout << "Kmer_length: " << Kmer_length << endl;

	helper = new ConstantLengthKmerHelper<keyT,valueT>(Kmer_length, splitbit);

	cout << "start to index Kmers" << endl;
	MulOthIndex<keyT> * moth;
	moth = new MulOthIndex<keyT>(argv[3], splitbit, helper, true);
	if (!moth->buildsucc) return 1;

	cout << "start to write indexes to file" << endl;
    moth->writeToFile(argv[4]);

    delete moth;
	delete helper;
	return 0;
}