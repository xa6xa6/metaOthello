// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
// 2, // level 0 -- root -- source: 0~2 at level 1
// 2,5,8, // level 1 -- source: NODE 0: 0~2 at level 2, NODE 1: 3~5 at level 2; NODE 2: 6~8 at level 2
// 2,5,8,11,14,17,20,23,26,
// level 2 -- source: NODE 0: 0~2 at level 3, NODE 1: 3~5 at level 3; NODE 2: 6~8 at level 3
// level 2 -- source: NODE 3: 9~11 at level 3, NODE 4: 12~14 at level 3; NODE 5: 15~17 at level 3
// level 2 -- source: NODE 6: 18~20 at level 3, NODE 7: 21~23 at level 3; NODE 8: 24~26 at level 3
#ifndef LOTHELLO_TREE_INFO_H
#define LOTHELLO_TREE_INFO_H
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
#include "../../mps3Lib/read_block_test.h"
#include "../../mps3Lib/otherFunc.h"
#include "../../mps3Lib/index_info.h"
#include "lothello_node_info.h"
using namespace std;

class Lothello_Tree_Info
{
private:
	int levelNum; // start from level_0 -- root
	vector<int> nodeNumVecInEachLevel;
	vector< vector< pair<int,int> > > lowerLevelSourceIndexRangeVecVec; 

	vector< vector<string> > nodeFileVecVec_jf;
	vector< vector<string> > nodeFileVecVec_readableKmer;
	
	//source could be nodes or raw samples
	vector< vector< vector<string> > > sourceFileVecVecVec_jf;
	vector< vector< vector<string> > > sourceFileVecVecVec_readableKmer;

public:
	Lothello_Tree_Info()
	{}

	void initiate(string& treeInfoFile, string& outputFolder_jf_readableKmer, vector<string>& providedRawSampleJfFilePathVec)
	{
		string cmd_mkdir_outputFolder_jf_readableKmer = "mkdir " + outputFolder_jf_readableKmer;
		system(cmd_mkdir_outputFolder_jf_readableKmer.c_str());
		this->initiate_levelNum_nodeNumVecInEachLevel_lowerLevelSourceIndexRangeVecVec(treeInfoFile);
		this->initiate_nodeFileVecVec_jf_readableKmer(outputFolder_jf_readableKmer);
		this->initiate_sourceFileVecVecVec_jf_readableKmer(outputFolder_jf_readableKmer, providedRawSampleJfFilePathVec);
	}

	void initiate_levelNum_nodeNumVecInEachLevel_lowerLevelSourceIndexRangeVecVec(string& treeInfoFile)
	{
		int tmpLevelNum = 0;
		ifstream treeInfo_ifs(treeInfoFile.c_str());
		while(!treeInfo_ifs.eof())
		{
			string tmpStr;
			getline(treeInfo_ifs, tmpStr);
			if(tmpStr == "")
				break;
			vector<int> tmpLastLowerLevelSourceIndexVec;
			this->parseStr2intVec_comma(tmpLastLowerLevelSourceIndexVec, tmpStr);
			
			int tmpLevelNodeNum = tmpLastLowerLevelSourceIndexVec.size();
			nodeNumVecInEachLevel.push_back(tmpLevelNodeNum); // initiate nodeNumVecInEachLevel
			
			vector< pair<int,int> > tmpLowerLevelSourceIndexRangeVec;
			int tmpSourceIndexRange_start = 0;
			for(int tmpLevelNodeIndex = 0; tmpLevelNodeIndex < tmpLevelNodeNum; tmpLevelNodeIndex ++)
			{
				int tmpSourceIndexRange_end = tmpLastLowerLevelSourceIndexVec[tmpLevelNodeIndex];
				tmpLowerLevelSourceIndexRangeVec.push_back(pair<int,int>
					(tmpSourceIndexRange_start, tmpSourceIndexRange_end));
				tmpSourceIndexRange_start = tmpSourceIndexRange_end + 1;
			}
			lowerLevelSourceIndexRangeVecVec.push_back(tmpLowerLevelSourceIndexRangeVec);
			
			tmpLevelNum ++;
		}
		levelNum = tmpLevelNum; // levelNum
		treeInfo_ifs.close();
	}

	void initiate_nodeFileVecVec_jf_readableKmer(string& outputFolder_jf_readableKmer)
	{
		for(int tmpLevelIndex = 0; tmpLevelIndex < levelNum; tmpLevelIndex ++)
		{
			int tmpLevelNodeNum = nodeNumVecInEachLevel[tmpLevelIndex];
			vector<string> tmpLevelJfFileVec;
			vector<string> tmpLevelReadableKmerFileVec;
			for(int tmpNodeIndex = 0; tmpNodeIndex < tmpLevelNodeNum; tmpNodeIndex ++)
			{
				string tmp_jf_file_path = outputFolder_jf_readableKmer + "level_" 
					+ int_to_str(tmpLevelIndex) + "_" + int_to_str(tmpNodeIndex) + ".jf";
				string tmp_readableKmer_path = outputFolder_jf_readableKmer + "level_"
					+ int_to_str(tmpLevelIndex) + "_" + int_to_str(tmpNodeIndex) + ".readableKmer";
				tmpLevelJfFileVec.push_back(tmp_jf_file_path);
				tmpLevelReadableKmerFileVec.push_back(tmp_readableKmer_path);
			}
			nodeFileVecVec_jf.push_back(tmpLevelJfFileVec);
			nodeFileVecVec_readableKmer.push_back(tmpLevelReadableKmerFileVec);
		}
	}

	void initiate_sourceFileVecVecVec_jf_readableKmer(string& outputFolder_jf_readableKmer, vector<string>& providedRawSampleJfFilePathVec)
	{
		// nonLeaf nodes -- root & internal nodes
		for(int tmpLevelIndex = 0; tmpLevelIndex < levelNum - 1; tmpLevelIndex ++)
		{
			vector< vector<string> > tmpLevelSourceFileVecVec_jf;
			vector< vector<string> > tmpLevelSourceFileVecVec_readableKmer;
			int tmpLevelNodeNum = nodeNumVecInEachLevel[tmpLevelIndex];
			for(int tmpNodeIndex = 0; tmpNodeIndex < tmpLevelNodeNum; tmpNodeIndex ++)
			{			
				int tmpSourceIndexRange_start = (lowerLevelSourceIndexRangeVecVec[tmpLevelIndex])[tmpNodeIndex].first;
				int tmpSourceIndexRange_end = (lowerLevelSourceIndexRangeVecVec[tmpLevelIndex])[tmpNodeIndex].second;
				vector<string> tmpNodeSourceFileVec_jf;
				vector<string> tmpNodeSourceFileVec_readableKmer;
				for(int tmpSourceFileIndex = tmpSourceIndexRange_start; tmpSourceFileIndex <= tmpSourceIndexRange_end; tmpSourceFileIndex ++)
				{
					string tmpSourceFile_jf = outputFolder_jf_readableKmer + "level_" + int_to_str(tmpLevelIndex + 1) 
						+ "_" + int_to_str(tmpSourceFileIndex) + ".jf";
					string tmpSourceFile_readableKmer = outputFolder_jf_readableKmer + "level_" + int_to_str(tmpLevelIndex + 1) 
						+ "_" + int_to_str(tmpSourceFileIndex) + ".readableKmer";						
					tmpNodeSourceFileVec_jf.push_back(tmpSourceFile_jf);
					tmpNodeSourceFileVec_readableKmer.push_back(tmpSourceFile_readableKmer);
				}	
				tmpLevelSourceFileVecVec_jf.push_back(tmpNodeSourceFileVec_jf);
				tmpLevelSourceFileVecVec_readableKmer.push_back(tmpNodeSourceFileVec_readableKmer);
			}
			sourceFileVecVecVec_jf.push_back(tmpLevelSourceFileVecVec_jf);
			sourceFileVecVecVec_readableKmer.push_back(tmpLevelSourceFileVecVec_readableKmer);
		}
		// leaf nodes
		int leafLevelIndex = levelNum - 1;
		int leafNodeNum = nodeNumVecInEachLevel[levelNum - 1];
		vector< vector<string> > leafLevelSourceFileVecVec_jf;
		vector< vector<string> > leafLevelSourceFileVecVec_readableKmer;
		for(int tmpNodeIndex = 0; tmpNodeIndex < leafNodeNum; tmpNodeIndex ++)
		{
			int tmpRawSampleSourceIndexRange_start = (lowerLevelSourceIndexRangeVecVec[leafLevelIndex])[tmpNodeIndex].first;
			int tmpRawSampleSourceIndexRange_end = (lowerLevelSourceIndexRangeVecVec[leafLevelIndex])[tmpNodeIndex].second;			
			vector<string> leafNodeSourceFileVec_jf;
			vector<string> leafNodeSourceFileVec_readableKmer;
			for(int tmpSourceFileIndex = tmpRawSampleSourceIndexRange_start; tmpSourceFileIndex <= tmpRawSampleSourceIndexRange_end; tmpSourceFileIndex ++)
			{
				leafNodeSourceFileVec_jf.push_back(providedRawSampleJfFilePathVec[tmpSourceFileIndex]);
				leafNodeSourceFileVec_readableKmer.push_back("NULL");
			}
			leafLevelSourceFileVecVec_jf.push_back(leafNodeSourceFileVec_jf);
			leafLevelSourceFileVecVec_readableKmer.push_back(leafNodeSourceFileVec_readableKmer);
		}
		sourceFileVecVecVec_jf.push_back(leafLevelSourceFileVecVec_jf);
		sourceFileVecVecVec_readableKmer.push_back(leafLevelSourceFileVecVec_readableKmer);		
	}

	void generate_lothello_allNodes(string& tmpJellyFishBinDirPath, string& tmpIntermediateDir)
	{
		// leaf node
		int leafLevelIndex = levelNum - 1;
		int leafNodeNum = nodeNumVecInEachLevel[leafLevelIndex];
		for(int tmpNodeIndex = 0; tmpNodeIndex < leafNodeNum; tmpNodeIndex++)
		{
			string tmpNode_intermediateDir = tmpIntermediateDir + "/tmpIntermediateDir_L_" + int_to_str(leafLevelIndex) + "_" + int_to_str(tmpNodeIndex); 
			this->generate_lothello_node_readableKmer_leaf(tmpNodeIndex, tmpJellyFishBinDirPath, tmpNode_intermediateDir);
		}
		// nonLeaf node
		for(int tmpLevelIndex = leafLevelIndex - 1; tmpLevelIndex >= 0; tmpLevelIndex --)
		{
			int tmpLevelNodeNum = nodeNumVecInEachLevel[tmpLevelIndex];
			for(int tmpNodeIndex = 0; tmpNodeIndex < tmpLevelNodeNum; tmpNodeIndex ++)
			{
				string tmpNode_intermediateDir = tmpIntermediateDir + "/tmpIntermediateDir_L_" + int_to_str(tmpLevelIndex) + "_" + int_to_str(tmpNodeIndex);
				this->generate_lothello_node_readableKmer_nonLeaf(tmpLevelIndex, tmpNodeIndex, tmpJellyFishBinDirPath, tmpNode_intermediateDir);
			}
		}
	}

	void generate_lothello_node_readableKmer_leaf(int leafNodeIndex, string& tmpJellyFishBinDirPath, string& tmpIntermediateDir)
	{
		int leafLevelIndex = levelNum - 1;
		Lothello_Node_Info* tmpNodeInfo = new Lothello_Node_Info();
		tmpNodeInfo->initaite_lothello_node_leaf((sourceFileVecVecVec_jf[leafLevelIndex])[leafNodeIndex],
			(nodeFileVecVec_readableKmer[leafLevelIndex])[leafNodeIndex], (nodeFileVecVec_jf[leafLevelIndex])[leafNodeIndex], tmpJellyFishBinDirPath);
		tmpNodeInfo->generateAndPrintThisNodeKmerClassInfoFile_leaf(tmpIntermediateDir);
		delete tmpNodeInfo;
	}

	void generate_lothello_node_readableKmer_nonLeaf(int nonLeafLevelIndex, int nonLeafNodeIndex, string& tmpJellyFishBinDirPath, string& tmpIntermediateDir)
	{
		Lothello_Node_Info* tmpNodeInfo = new Lothello_Node_Info();
		tmpNodeInfo->initaite_lothello_node_nonLeaf((sourceFileVecVecVec_jf[nonLeafLevelIndex])[nonLeafNodeIndex],
			(sourceFileVecVecVec_readableKmer[nonLeafLevelIndex])[nonLeafNodeIndex], (nodeFileVecVec_readableKmer[nonLeafLevelIndex])[nonLeafNodeIndex], 
			(nodeFileVecVec_jf[nonLeafLevelIndex])[nonLeafNodeIndex], tmpJellyFishBinDirPath);
		tmpNodeInfo->generateAndPrintThisNodeKmerClassInfoFile_nonLeaf(tmpIntermediateDir);
		delete tmpNodeInfo;
	}

	void print_lothelloTreeStructure(string& outputDir)
	{
		outputDir += "/";
		string cmd_mkdir = "mkdir " + outputDir;
		system(cmd_mkdir.c_str());
		string output_treeStats = outputDir + "tree.stats";
		string output_treeStructure = outputDir + "tree.structure";
		string output_nodeJfReadableKmerFilePath = outputDir + "node_jf_readableKmer.path";
		string output_sourceJfReadableKmerFilePath = outputDir + "source_jf_readableKmer.path";

		// print tree stats
		ofstream treeStats_ofs(output_treeStats.c_str());
		treeStats_ofs << "Level #: " << levelNum << endl << endl;
		treeStats_ofs << "Total Node #: " << this->sum_intVec(nodeNumVecInEachLevel) << endl;
		for(int tmpLevelIndex = 0; tmpLevelIndex < levelNum; tmpLevelIndex ++)
			treeStats_ofs << "     Level " << tmpLevelIndex << ": " << nodeNumVecInEachLevel[tmpLevelIndex] << endl;  
		treeStats_ofs.close();

		// print tree structure
		ofstream treeStructure_ofs(output_treeStructure.c_str());
		treeStructure_ofs << "Lothello Tree Root (L_0_0)" << endl;
		// nonLeaf node
		string tmpIndent = "";
		for(int tmpLevelIndex = 0; tmpLevelIndex < levelNum - 1; tmpLevelIndex ++)
		{
			int tmpLevelNodeNum = nodeNumVecInEachLevel[tmpLevelIndex];
			int tmpNextLevelIndex = tmpLevelIndex + 1;
			for(int tmpNodeIndex = 0; tmpNodeIndex < tmpLevelNodeNum; tmpNodeIndex ++)
			{
				string tmpNodeId = "L_" + int_to_str(tmpLevelIndex) + "_" + int_to_str(tmpNodeIndex);
				int tmpSourceIndexRange_start = (lowerLevelSourceIndexRangeVecVec[tmpLevelIndex])[tmpNodeIndex].first;
 				int tmpSourceIndexRange_end = (lowerLevelSourceIndexRangeVecVec[tmpLevelIndex])[tmpNodeIndex].second;
				for(int tmpSourceIndex = tmpSourceIndexRange_start; tmpSourceIndex <= tmpSourceIndexRange_end; tmpSourceIndex ++)
				{
					string tmpSourceId = "L_" + int_to_str(tmpNextLevelIndex) + "_" + int_to_str(tmpSourceIndex);
					treeStructure_ofs << tmpIndent << tmpNodeId << "<--" << tmpSourceId << endl;
				}
			}
			tmpIndent += "\t";
		}
		// leaf node
		int leafLevelIndex = levelNum - 1;
		int leafNodeNum = nodeNumVecInEachLevel[leafLevelIndex];
		for(int tmpNodeIndex = 0; tmpNodeIndex < leafNodeNum; tmpNodeIndex ++)
		{
			int tmpSourceIndexRange_start = (lowerLevelSourceIndexRangeVecVec[leafLevelIndex])[tmpNodeIndex].first;
 			int tmpSourceIndexRange_end = (lowerLevelSourceIndexRangeVecVec[leafLevelIndex])[tmpNodeIndex].second;
			for(int tmpSourceIndex = tmpSourceIndexRange_start; tmpSourceIndex <= tmpSourceIndexRange_end; tmpSourceIndex ++)
				treeStructure_ofs << tmpIndent << "L_" << leafLevelIndex << "_" << tmpNodeIndex << "<--Sample_" << tmpSourceIndex << endl;
		}
		treeStructure_ofs.close();

		// print node_jf_readableKmer
		ofstream node_jf_readableKmer_ofs(output_nodeJfReadableKmerFilePath.c_str());
		for(int tmpLevelIndex = 0; tmpLevelIndex < levelNum; tmpLevelIndex ++)
		{
			int tmpLevelNodeNum = nodeNumVecInEachLevel[tmpLevelIndex];
			for(int tmpNodeIndex = 0; tmpNodeIndex < tmpLevelNodeNum; tmpNodeIndex ++)
			{
				node_jf_readableKmer_ofs << "L_" << tmpLevelIndex << "_" << tmpNodeIndex << "\t"
					<< (nodeFileVecVec_jf[tmpLevelIndex])[tmpNodeIndex] << "\t" 
					<< (nodeFileVecVec_readableKmer[tmpLevelIndex])[tmpNodeIndex] << endl;
			}
		}
		node_jf_readableKmer_ofs.close();
		// print source_jf_readableKmer

		ofstream source_jf_readableKmer_ofs(output_sourceJfReadableKmerFilePath.c_str());
		for(int tmpLevelIndex = 0; tmpLevelIndex < levelNum; tmpLevelIndex ++)
		{
			int tmpLevelNodeNum = nodeNumVecInEachLevel[tmpLevelIndex];
			for(int tmpNodeIndex = 0; tmpNodeIndex < tmpLevelNodeNum; tmpNodeIndex ++)
			{
				source_jf_readableKmer_ofs << "L_" << tmpLevelIndex << "_" << tmpNodeIndex << ":" << endl; 
				int tmpSourceIndexRange_start = (lowerLevelSourceIndexRangeVecVec[leafLevelIndex])[tmpNodeIndex].first;
 				int tmpSourceIndexRange_end = (lowerLevelSourceIndexRangeVecVec[leafLevelIndex])[tmpNodeIndex].second;
				int tmpSourceNum = tmpSourceIndexRange_end - tmpSourceIndexRange_start + 1;
				for(int tmpSourceIndex = tmpSourceIndexRange_start; tmpSourceIndex <= tmpSourceIndexRange_end; tmpSourceIndex ++)
					source_jf_readableKmer_ofs << ((sourceFileVecVecVec_jf[tmpLevelIndex])[tmpNodeIndex])[tmpSourceIndex] 
						<< "\t" << ((sourceFileVecVecVec_readableKmer[tmpLevelIndex])[tmpNodeIndex])[tmpSourceIndex] << endl;
			}
		}
		source_jf_readableKmer_ofs.close();
	}

	void parseStr2intVec_comma(vector<int>& tmpIntVec, string& tmpStr)
	{
		int tmpStartLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tmpCommaLoc = tmpStr.find(",", tmpStartLoc);
			if(tmpCommaLoc == string::npos)
				break;
			string tmpIntStr = tmpStr.substr(tmpStartLoc, tmpCommaLoc - tmpStartLoc);
			int tmpInt = atoi(tmpIntStr.c_str());
			tmpIntVec.push_back(tmpInt);
			tmpStartLoc = tmpCommaLoc + 1;
		}		
	}	

	int sum_intVec(vector<int>& tmpIntVec)
	{
		int tmpSum = 0;
		for(int tmp = 0; tmp < tmpIntVec.size(); tmp++)
			tmpSum += tmpIntVec[tmp];
		return tmpSum;
	}
};
#endif