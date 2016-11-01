#ifndef MERGEJF_INFO_H
#define MERGEJF_INFO_H
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

class MergeJf_Info
{
private:
	vector<string> inputJfFileVec;
	string outputMergedJfFile;
	string jellyFishBinPath;
public:
	MergeJf_Info()
	{}

	void initaite(vector<string>& tmpInputJfFileVec, string& tmpOutputMergedJfFile, string& tmpJellyFishBinPath)
	{
		//string inputJfFile = tmpInputJfFile;
		for(int tmp = 0; tmp < tmpInputJfFileVec.size(); tmp++)
		{
			string tmpInputJfFile = tmpInputJfFileVec[tmp];
			inputJfFileVec.push_back(tmpInputJfFile);
		}
		string outputMergedJfFile = tmpOutputMergedJfFile;
		string jellyFishBinPath = tmpJellyFishBinPath;
	}

	void mergeJf()
	{
		string cmd_merge = jellyFishBinPath + "/jellyfish merge -o " + outputMergedJfFile;
		for(int tmp = 0; tmp < inputJfFileVec.size(); tmp++)
		{
			cmd_merge += " ";
			cmd_merge += inputJfFileVec[tmp];
		}
		system(cmd_merge.c_str());
	}

	string return_outputMergedJfFile()
	{
		return outputMergedJfFile;
	}
};
#endif