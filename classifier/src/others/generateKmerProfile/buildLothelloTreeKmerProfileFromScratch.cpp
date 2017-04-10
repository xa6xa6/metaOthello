// 2, // level 0 -- root -- source: 0~2 at level 1
// 2,5,8, // level 1 -- source: NODE 0: 0~2 at level 2, NODE 1: 3~5 at level 2; NODE 2: 6~8 at level 2
// 2,5,8,11,14,17,20,23,26,
// level 2 -- source: NODE 0: 0~2 at level 3, NODE 1: 3~5 at level 3; NODE 2: 6~8 at level 3
// level 2 -- source: NODE 3: 9~11 at level 3, NODE 4: 12~14 at level 3; NODE 5: 15~17 at level 3
// level 2 -- source: NODE 6: 18~20 at level 3, NODE 7: 21~23 at level 3; NODE 8: 24~26 at level 3
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
#include "general/lothello_tree_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable jellyFishBinPath inputTreeInfo inputJfPathListFile outputFolder" << endl;
		exit(1);
	}
	string outputFolder = argv[4];
	outputFolder += "/";
	string cmd_mkdir_outputFolder = "mkdir " + outputFolder;
	system(cmd_mkdir_outputFolder.c_str());
	string log_file = outputFolder + "/log.txt";
	ofstream log_ofs(log_file.c_str());	

	log_ofs << "start to initiate jfPathVec from jfPathListFile" << endl;
	string inputJfPathListFile = argv[3];
	vector<string> jfPathVec;
	ifstream jfPathList_ifs(inputJfPathListFile.c_str());
	while(!jfPathList_ifs.eof())
	{
		string tmpStr;
		getline(jfPathList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		jfPathVec.push_back(tmpStr);
	}
	jfPathList_ifs.close();

	log_ofs << "start to initiate lothello Tree Info" << endl;
	string inputTreeInfo = argv[2];
	string outputFolder_jf_readableKmer = outputFolder + "jf_readableKmer";	
	Lothello_Tree_Info lothelloTreeInfo;
	lothelloTreeInfo.initiate(inputTreeInfo, outputFolder_jf_readableKmer, jfPathVec);	

	log_ofs << "start to print out lothello Tree Info" << endl;
	string outputFolder_treeInfo = outputFolder + "tree_info";
	lothelloTreeInfo.print_lothelloTreeStructure(outputFolder_treeInfo);

	string jellyFishBinPath = argv[1];
	string outputFolder_tmpIntermediateDir = outputFolder + "tmpIntermediateDir";
	string cmd_mkdir_outputFolder_tmpIntermediateDir = "mkdir " + outputFolder_tmpIntermediateDir;
	system(cmd_mkdir_outputFolder_tmpIntermediateDir.c_str());
	lothelloTreeInfo.generate_lothello_allNodes(jellyFishBinPath, outputFolder_tmpIntermediateDir);

	log_ofs.close();
	return 0;
}