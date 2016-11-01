#ifndef JF2KMER_INFO_H
#define JF2KMER_INFO_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>

using namespace std;

class Jf2Kmer_Info
{
private:
	string inputJfFile;
	string outputKmerFile;
	string jellyFishBinPath;
public:
	Jf2Kmer_Info()
	{}

	void initiate(string& tmpInputJfFile, string& tmpOutputKmerFile, string& tmpJellyFishBinPath)
	{
		inputJfFile = tmpInputJfFile;
		outputKmerFile = tmpOutputKmerFile;
		jellyFishBinPath = tmpJellyFishBinPath;
	}

	void jf2Kmer()
	{
		string cmd_dump = jellyFishBinPath + "/jellyfish dump -t -c -o " + outputKmerFile + " " + inputJfFile;
		system(cmd_dump.c_str());
	}

	string return_outputKmerFile()
	{
		return outputKmerFile;
	}
};
#endif