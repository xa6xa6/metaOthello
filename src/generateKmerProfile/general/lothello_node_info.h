#ifndef LOTHELLO_NODE_INFO_H
#define LOTHELLO_NODE_INFO_H
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

class Lothello_Node_Info
{
private:
	int branchNumFromThisNode; // leaf & nonLeaf (root & internal)
	//input
	vector<string> lowerLevelNodeJfFileVec; // leaf & nonLeaf (root & internal)
	vector<string> lowerLevelNodeSortedReadableKmerFileVec; // nonLeaf: root & internal, could be kmerCountFile or kmerClassFile
	// output
	string thisNodeKmerClassInfoFile; // leaf & nonLeaf (root & internal)
	string thisNodeMergedJfFile; // leaf & nonLeaf (root & internal)
	// jellyfish tool bin path
	string jellyFishBinDirPath;
public:
	Lothello_Node_Info()
	{}

	void initaite_lothello_node_leaf(vector<string>& tmpInputJfFileVec, 
		string& tmpKmerClassInfoOutputFile, string& tmpMergedJfOutputFile, string& tmpJellyFishBinDirPath)
	{
		branchNumFromThisNode = tmpInputJfFileVec.size();
		for(int tmp = 0; tmp < branchNumFromThisNode; tmp++)
			lowerLevelNodeJfFileVec.push_back(tmpInputJfFileVec[tmp]);
		thisNodeKmerClassInfoFile = tmpKmerClassInfoOutputFile;
		thisNodeMergedJfFile = tmpMergedJfOutputFile;
		jellyFishBinDirPath = tmpJellyFishBinDirPath;
	}

	void initaite_lothello_node_nonLeaf(vector<string>& tmpJfFileVec, vector<string>& tmpSortedReadableKmerFileVec, 
		string& tmpKmerClassInfoOutputFile, string& tmpMergedJfOutputFile, string& tmpJellyFishBinDirPath)
	{
		branchNumFromThisNode = tmpJfFileVec.size();
		for(int tmp = 0; tmp < branchNumFromThisNode; tmp++)
		{
			lowerLevelNodeJfFileVec.push_back(tmpJfFileVec[tmp]);
			lowerLevelNodeSortedReadableKmerFileVec.push_back(tmpSortedReadableKmerFileVec[tmp]);
		}
		thisNodeKmerClassInfoFile = tmpKmerClassInfoOutputFile;
		thisNodeMergedJfFile = tmpMergedJfOutputFile;
		jellyFishBinDirPath = tmpJellyFishBinDirPath;
	}

	void generateAndPrintThisNodeKmerClassInfoFile_leaf(string& tmpDir)
	{
		// mkdir
		string cmd_mkdir_tmpDir = "mkdir " + tmpDir;
		system(cmd_mkdir_tmpDir.c_str());
		string tmpDir_sort = tmpDir + "/tmpSortDir/";
		string cmd_mkdir_tmpDir_sort = "mkdir " + tmpDir_sort;
		system(cmd_mkdir_tmpDir_sort.c_str());
		string tmpDir_readableKmer = tmpDir + "/tmpReadableKmer/";
		string cmd_mkdir_tmpDir_readableKmer = "mkdir " + tmpDir_readableKmer;
		system(cmd_mkdir_tmpDir_readableKmer.c_str());
		string tmpDir_readableKmer_sorted = tmpDir + "/tmpReadableKmer_sorted/";
		string cmd_mkdir_tmpDir_readableKmer_sorted =  "mkdir " + tmpDir_readableKmer_sorted;
		system(cmd_mkdir_tmpDir_readableKmer_sorted.c_str());

		// start to merge jf files to an union one
		string cmd_jellyfish_merge = jellyFishBinDirPath + "/jellyfish merge -o " + thisNodeMergedJfFile;
		for(int tmp = 0; tmp < branchNumFromThisNode; tmp++)
		{
			cmd_jellyfish_merge += " ";
			cmd_jellyfish_merge += lowerLevelNodeJfFileVec[tmp];
		}
		system(cmd_jellyfish_merge.c_str());

		// start to convert jf files in readable format: <kmer_seq> + "/t" + <kmerCount>
		string cmd_jellyfish_dump_unionKmerFile = jellyFishBinDirPath + "/jellyfish dump -t -c -o " 
			+ tmpDir + "/merged.jf.readable " + thisNodeMergedJfFile;
		system(cmd_jellyfish_dump_unionKmerFile.c_str());
		for(int tmp = 0; tmp < branchNumFromThisNode; tmp++)
		{
			string tmp_cmd_jf_dump = jellyFishBinDirPath + "/jellyfish dump -t -c -o "
				+ tmpDir_readableKmer + int_to_str(tmp) + ".readable " + lowerLevelNodeJfFileVec[tmp];
			system(tmp_cmd_jf_dump.c_str());
		}

		// start to sort each kmer count file in lexicographical order
		string cmd_sort_unionKmerFile = "sort -k1 -T=" + tmpDir_sort + " " + tmpDir 
			+ "/merged.jf.readable > " + tmpDir + "/merged.jf.readable.sorted";
		system(cmd_sort_unionKmerFile.c_str());
		for(int tmp = 0; tmp < branchNumFromThisNode; tmp++)
		{
			string tmp_cmd_sort_single_KmerFile = "sort -k1 -T=" + tmpDir_sort + " " 
				+ tmpDir_readableKmer + int_to_str(tmp) + ".readable > "
				+ tmpDir_readableKmer_sorted + int_to_str(tmp) + ".readable.sorted";
			system(tmp_cmd_sort_single_KmerFile.c_str());
		}

		// generate Kmer class info file
		string inputMergedKmerReadableFile = tmpDir + "/merged.jf.readable.sorted";
		vector<string> KmerReadableFilePathVec;
		for(int tmp = 0; tmp < branchNumFromThisNode; tmp++)
		{
			string tmpReadableKmer_sorted_file = tmpDir_readableKmer_sorted + int_to_str(tmp) + ".readable.sorted";
			KmerReadableFilePathVec.push_back(tmpReadableKmer_sorted_file);
		}

		// assign class info
		int providedKmerReadableFileNum = branchNumFromThisNode;
		vector<ifstream*> KmerReadableIfsVec;
		for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp++)
		{
			ifstream *tmpKmerReadableIfs = new ifstream(KmerReadableFilePathVec[tmp].c_str());
			KmerReadableIfsVec.push_back(tmpKmerReadableIfs);
		}		
		ofstream kmerClass_ofs(thisNodeKmerClassInfoFile.c_str());
		ifstream mergedKmer_ifs(inputMergedKmerReadableFile.c_str());

		vector<bool> kmer_file_end_bool_vec;
		for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp ++)
			kmer_file_end_bool_vec.push_back(false);

		vector<string> currentKmerStrVec;
		for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp ++)
			currentKmerStrVec.push_back("");
		//cout << "head lines" << endl;
		for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp ++)
		{
			string tmpKmerFileStr;
			getline((*KmerReadableIfsVec[tmp]), tmpKmerFileStr);
			//cout << "tmpKmerFileStr: " << tmpKmerFileStr << endl;
			int tmpTabLoc = tmpKmerFileStr.find("\t");
			currentKmerStrVec[tmp] = tmpKmerFileStr.substr(0, tmpTabLoc);
			//cout << "tmpCurrentKmerStr: " << tmp << endl << currentKmerStrVec[tmp] << endl;
		}

		while(!mergedKmer_ifs.eof())
		{
			string tmpStr;
			getline(mergedKmer_ifs, tmpStr);
			if(tmpStr == "")
				break;
			int tabLoc = tmpStr.find("\t");
			string tmpMergedStr = tmpStr.substr(0, tabLoc);
			vector<bool> tmp_kmer_exist_in_sample_bool_vec;
			for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp ++)
				tmp_kmer_exist_in_sample_bool_vec.push_back(false);

			int tmpKmerExistingClassFlag = 0;
			for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp ++)
			{
				if(tmpMergedStr == currentKmerStrVec[tmp])
				{
					tmp_kmer_exist_in_sample_bool_vec[tmp] = true;
					tmpKmerExistingClassFlag += pow(2, tmp);
					if(!kmer_file_end_bool_vec[tmp])
					{
						if((*KmerReadableIfsVec[tmp]).eof())
							kmer_file_end_bool_vec[tmp] = true;
						else
						{	
							string tmpKmerFileStr;
							getline((*KmerReadableIfsVec[tmp]), tmpKmerFileStr);
							if(tmpKmerFileStr == "")
								kmer_file_end_bool_vec[tmp] = true;
							else
							{
								int tmpTabLoc = tmpKmerFileStr.find("\t");
								currentKmerStrVec[tmp] = tmpKmerFileStr.substr(0, tmpTabLoc);
							}
						}
					}
				}
				else
					tmp_kmer_exist_in_sample_bool_vec[tmp] = false;
			}

			if(tmpKmerExistingClassFlag > 0)
				kmerClass_ofs << tmpMergedStr << "\t" << tmpKmerExistingClassFlag << endl;
			else
			{
				cout << "error ! tmpKmerExistingClassFlag <= 0" << endl;
				cout << "tmpMergedStr: " << tmpMergedStr << endl;
				exit(1);
			}
		}

		// remove intermediate dirs and files
		string cmd_rm_dir = "rm -r " + tmpDir;
		//system(cmd_rm_dir.c_str());

		for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp++)
		{
			(*KmerReadableIfsVec[tmp]).close();
			delete KmerReadableIfsVec[tmp];
		}
		mergedKmer_ifs.close();
		kmerClass_ofs.close();
	}

	void generateAndPrintThisNodeKmerClassInfoFile_nonLeaf(string& tmpDir)
	{
		// mkdir
		string cmd_mkdir_tmpDir = "mkdir " + tmpDir;
		system(cmd_mkdir_tmpDir.c_str());
		string tmpDir_sort = tmpDir + "/tmpSortDir/";
		string cmd_mkdir_tmpDir_sort = "mkdir " + tmpDir_sort;
		system(cmd_mkdir_tmpDir_sort.c_str());

		// start to merge jf files to an union one
		string cmd_jellyfish_merge = jellyFishBinDirPath + "/jellyfish merge -o " + thisNodeMergedJfFile;
		for(int tmp = 0; tmp < branchNumFromThisNode; tmp++)
		{
			cmd_jellyfish_merge += " ";
			cmd_jellyfish_merge += lowerLevelNodeJfFileVec[tmp];
		}
		system(cmd_jellyfish_merge.c_str());

		// start to convert the merged jf file to readable format: <kmer_seq> + "/t" + <kmerCount>
		string cmd_jellyfish_dump_unionKmerFile = jellyFishBinDirPath + "/jellyfish dump -t -c -o " 
			+ tmpDir + "/merged.jf.readable " + thisNodeMergedJfFile;
		system(cmd_jellyfish_dump_unionKmerFile.c_str());

		// start to sort merged readable format Kmer file in lexicographical order (the 2nd field each line is the classInfo, not the count)
		string cmd_sort_unionKmerFile = "sort -k1 -T=" + tmpDir_sort + " " + tmpDir 
			+ "/merged.jf.readable > " + tmpDir + "/merged.jf.readable.sorted";
		system(cmd_sort_unionKmerFile.c_str());

		// generate Kmer class info file
		string inputMergedKmerReadableFile = tmpDir + "/merged.jf.readable.sorted";
		vector<string> KmerReadableFilePathVec;
		for(int tmp = 0; tmp < branchNumFromThisNode; tmp++)
			KmerReadableFilePathVec.push_back(lowerLevelNodeSortedReadableKmerFileVec[tmp]);

		// assign class info
		int providedKmerReadableFileNum = branchNumFromThisNode;
		vector<ifstream*> KmerReadableIfsVec;
		for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp++)
		{
			ifstream *tmpKmerReadableIfs = new ifstream(KmerReadableFilePathVec[tmp].c_str());
			KmerReadableIfsVec.push_back(tmpKmerReadableIfs);
		}		
		ofstream kmerClass_ofs(thisNodeKmerClassInfoFile.c_str());
		ifstream mergedKmer_ifs(inputMergedKmerReadableFile.c_str());

		vector<bool> kmer_file_end_bool_vec;
		for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp ++)
			kmer_file_end_bool_vec.push_back(false);

		vector<string> currentKmerStrVec;
		for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp ++)
			currentKmerStrVec.push_back("");

		for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp ++)
		{
			string tmpKmerFileStr;
			getline((*KmerReadableIfsVec[tmp]), tmpKmerFileStr);
			int tmpTabLoc = tmpKmerFileStr.find("\t");
			currentKmerStrVec[tmp] = tmpKmerFileStr.substr(0, tmpTabLoc);
		}

		while(!mergedKmer_ifs.eof())
		{
			string tmpStr;
			getline(mergedKmer_ifs, tmpStr);
			if(tmpStr == "")
				break;
			int tabLoc = tmpStr.find("\t");
			string tmpMergedStr = tmpStr.substr(0, tabLoc);
			vector<bool> tmp_kmer_exist_in_sample_bool_vec;
			for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp ++)
				tmp_kmer_exist_in_sample_bool_vec.push_back(false);

			int tmpKmerExistingClassFlag = 0;
			for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp ++)
			{
				if(tmpMergedStr == currentKmerStrVec[tmp])
				{
					tmp_kmer_exist_in_sample_bool_vec[tmp] = true;
					tmpKmerExistingClassFlag += pow(2, tmp);
					if(!kmer_file_end_bool_vec[tmp])
					{
						if((*KmerReadableIfsVec[tmp]).eof())
							kmer_file_end_bool_vec[tmp] = true;
						else
						{	
							string tmpKmerFileStr;
							getline((*KmerReadableIfsVec[tmp]), tmpKmerFileStr);
							if(tmpKmerFileStr == "")
								kmer_file_end_bool_vec[tmp] = true;
							else
							{
								int tmpTabLoc = tmpKmerFileStr.find("\t");
								currentKmerStrVec[tmp] = tmpKmerFileStr.substr(0, tmpTabLoc);
							}
						}
					}
				}
				else
					tmp_kmer_exist_in_sample_bool_vec[tmp] = false;
			}

			if(tmpKmerExistingClassFlag > 0)
				kmerClass_ofs << tmpMergedStr << "\t" << tmpKmerExistingClassFlag << endl;
			else
			{
				cout << "error ! tmpKmerExistingClassFlag <= 0" << endl;
				cout << "tmpMergedStr: " << tmpMergedStr << endl;
				exit(1);
			}
		}

		// remove intermediate dirs and files
		string cmd_rm_dir = "rm -r " + tmpDir;
		//system(cmd_rm_dir.c_str());

		for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp++)
		{
			(*KmerReadableIfsVec[tmp]).close();
			delete KmerReadableIfsVec[tmp];
		}
		mergedKmer_ifs.close();
		kmerClass_ofs.close();	
	}
};
#endif