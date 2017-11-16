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
#include <hash_map>
#include <map>
#include <set>
#include "../mps3Lib/read_block_test.h"
#include "../mps3Lib/otherFunc.h"
#include "../mps3Lib/index_info.h"
#include "general/lothello_node_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 8)
	{
		cout << "Executable inputJellyfishBin nodeType outputKmerClassFile outputMergedJfFile num_of_lowerNodeOrSample"; 
		cout << " tmpDir jf_1 (jf_2 jf_3 ...) (readableKmer_1 (readableKmer_2 readableKmer_3 ...))" << endl;
		exit(1);
	}
	string inputJellyfishBin = argv[1];
	inputJellyfishBin += "/";
	string nodeTypeStr = argv[2];
	bool leaf_or_nonLeaf_bool;
	if((nodeTypeStr == "Leaf")||(nodeTypeStr == "LEAF")||(nodeTypeStr == "leaf"))
		leaf_or_nonLeaf_bool = true;
	else if((nodeTypeStr == "Root")||(nodeTypeStr == "ROOT")||(nodeTypeStr == "root")
		||(nodeTypeStr == "Internal")||(nodeTypeStr == "INTERNAL")||(nodeTypeStr == "internal"))
		leaf_or_nonLeaf_bool = false;
	else
	{
		cout << "invalid nodeType: " << nodeTypeStr << endl;
		cout << "Please specify: Leaf or Root or Internal !" << endl;
		exit(1);
	}

	string outputKmerClassFile = argv[3];
	string outputMergedJfFile = argv[4];
	string lowerNodeOrSampleNumStr = argv[5];
	int lowerNodeOrSampleNum = atoi(lowerNodeOrSampleNumStr.c_str());
	string tmpDir = argv[6];
	tmpDir += "/";
	if(leaf_or_nonLeaf_bool) // leaf node
	{
		int jf_file_num = argc - 7;
		if(jf_file_num != lowerNodeOrSampleNum)
		{
			cout << "error ! jf_file_num != lowerNodeOrSampleNum" << endl;
			exit(1);
		}
		vector<string> jfFileVec;
		for(int tmp = 0; tmp < jf_file_num; tmp++)
		{
			string tmpJfFile = argv[7+tmp];
			jfFileVec.push_back(tmpJfFile);
		}

		Lothello_Node_Info tmpNodeInfo; // leaf node
		cout << "start to initiate lothello node" << endl;
		tmpNodeInfo.initaite_lothello_node_leaf(jfFileVec, 
			outputKmerClassFile, outputMergedJfFile, inputJellyfishBin);
		cout << "start to generateAndPrintThisNodeKmerClassInfoFile" << endl;
		tmpNodeInfo.generateAndPrintThisNodeKmerClassInfoFile_leaf(tmpDir);
	}
	else // nonLeaf node, root or internal
	{
		int jf_file_num = (argc - 7)/2;
		if((jf_file_num * 2) != argc)
		{
			cout << "error ! (jf_file_num * 2) != argc" << endl;
			exit(1);
		}
		if(jf_file_num != lowerNodeOrSampleNum)
		{
			cout << "error ! jf_file_num != lowerNodeOrSampleNum" << endl;
			exit(1);
		}
		vector<string> jfFileVec;
		vector<string> readableKmerFileVec;
		for(int tmp = 0; tmp < jf_file_num; tmp++)
		{
			string tmpJfFile = argv[7+tmp];
			jfFileVec.push_back(tmpJfFile);
			string tmpReadableKmerFile = argv[7 + jf_file_num + tmp];
			readableKmerFileVec.push_back(tmpReadableKmerFile);
		}
		Lothello_Node_Info tmpNodeInfo; // nonLeaf node, root or internal
		cout << "start to initiate lothello node" << endl;
		tmpNodeInfo.initaite_lothello_node_nonLeaf(jfFileVec, readableKmerFileVec,
			outputKmerClassFile, outputMergedJfFile, inputJellyfishBin);
		cout << "start to generateAndPrintThisNodeKmerClassInfoFile" << endl;
		tmpNodeInfo.generateAndPrintThisNodeKmerClassInfoFile_nonLeaf(tmpDir);
	}
	return 0;
}