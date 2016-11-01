#ifndef FA2JF_INFO_H
#define FA2JF_INFO_H
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

class Fa2jf_Info
{
private:
	string inputFaFile;
	string outputJfFile;
	string jellyFishBinPath;
	string threads_num_str;
	string count_min_str;
	string Kmer_length_str;
	string bf_size_str;
public:
	Fa2jf_Info()
	{}

	void initaite_withMinCount(string& tmpInputFaFile, string& tmpOutputJfFile, string& tmpJellyFishBinPath,
		string& tmp_threads_num_str, string& tmp_count_min_str, string& tmp_Kmer_length_str, string& tmp_bf_size_str)
	{
		inputFaFile = tmpInputFaFile;
		outputJfFile = tmpOutputJfFile;
		jellyFishBinPath = tmpJellyFishBinPath;
		threads_num_str = tmp_threads_num_str;
		count_min_str = tmp_count_min_str;
		Kmer_length_str = tmp_Kmer_length_str;
		bf_size_str = tmp_bf_size_str;
	}

	void initaite_noMinCount(string& tmpInputFaFile, string& tmpOutputJfFile, string& tmpJellyFishBinPath,
		string& tmp_threads_num_str, string& tmp_Kmer_length_str, string& tmp_bf_size_str)
	{
		inputFaFile = tmpInputFaFile;
		outputJfFile = tmpOutputJfFile;
		jellyFishBinPath = tmpJellyFishBinPath;
		threads_num_str = tmp_threads_num_str;
		count_min_str = "1";
		Kmer_length_str = tmp_Kmer_length_str;
		bf_size_str = tmp_bf_size_str;
	}

	void fa2jf()
	{
		//string threads_num_str = int_to_str();
		if(count_min_str == "1")
		{
			string cmd_dump_noMinCount = jellyFishBinPath + "/jellyfish count -o " 
				+ outputJfFile + " -m " + Kmer_length_str + " -t " + threads_num_str 
				+ " -s " + bf_size_str + " -C " + inputFaFile;
			system(cmd_dump_noMinCount.c_str());	
		}
		else
		{	
			string cmd_dump_withMinCount = jellyFishBinPath + "/jellyfish count -o " + outputJfFile 
				+ " -m " + Kmer_length_str + " -t " + threads_num_str + " -s " + bf_size_str 
				+ " -C -L " + count_min_str + " " + inputFaFile;
			system(cmd_dump_withMinCount.c_str());
		}
	}

	string return_outputJfFile()
	{
		return outputJfFile;
	}
};
#endif