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

using namespace std;

void parseStr2fieldVec_tab(vector<string>& tmpFieldVec, string& tmpStr)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find("\t", startLoc);
		if(tabLoc == string::npos)
			break;
		string tmpField = tmpStr.substr(startLoc, tabLoc-startLoc);
		tmpFieldVec.push_back(tmpField);
		startLoc = tabLoc + 1;
	}
	tmpFieldVec.push_back(tmpStr.substr(startLoc));
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 fusionGroudTruth" << endl;
		cout << "#2 fusionDetectionFile" << endl;
		cout << "#3 outputDir" << endl;
		exit(1);
	}
	int gene_1_tab_index = 2;
	int gene_2_tab_index = 3;

	string fusionGroudTruth = argv[1];
	string fusionDetectionFile = argv[2];
	string outputDir = argv[3]; outputDir += "/";
	string cmd_mkdir = "mkdir " + outputDir;
	system(cmd_mkdir.c_str());
	string output_fusion_true_detected = outputDir + "detected_true.txt";
	string output_fusion_false_detected = outputDir + "detected_false.txt";
	string output_fusion_true_missed = outputDir + "missed_true.txt"; 
	cout << "start to load ground truth fusions " << endl;
	set<string> groundTruthFusionGenePairSet;
	ifstream gt_ifs(fusionGroudTruth.c_str());
	while(!gt_ifs.eof())
	{
		string tmpStr;
		getline(gt_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpGene_1 = tmpStr.substr(0, tabLoc);
		string tmpGene_2 = tmpStr.substr(tabLoc + 1);
		string tmpGenePair = tmpGene_1 + tmpGene_2;
		groundTruthFusionGenePairSet.insert(tmpGenePair);
	}
	gt_ifs.close();
	cout << "start to load detected fusion file" << endl;
	set<string> detected_true_genePair_set;
	ofstream dt_tr_ofs(output_fusion_true_detected.c_str());
	ofstream dt_fl_ofs(output_fusion_false_detected.c_str());
	ifstream dt_ifs(fusionDetectionFile.c_str());
	while(!dt_ifs.eof())
	{
		string tmpStr;
		getline(dt_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpFieldVec;
		parseStr2fieldVec_tab(tmpFieldVec, tmpStr);
		string tmpGene_1 = tmpFieldVec[gene_1_tab_index];
		string tmpGene_2 = tmpFieldVec[gene_2_tab_index];
		string tmpGenePair_1 = tmpGene_1 + tmpGene_2;
		string tmpGenePair_2 = tmpGene_2 + tmpGene_1;
		//cout << "tmpGenePair_1: " << tmpGenePair_1 << endl;
		//cout << "tmpGenePair_2: " << tmpGenePair_2 << endl;
		if(groundTruthFusionGenePairSet.find(tmpGenePair_1) != groundTruthFusionGenePairSet.end()) // tmpGenePair_1 found
		{
			detected_true_genePair_set.insert(tmpGenePair_1);
			dt_tr_ofs << tmpStr << endl;
		}
		else if(groundTruthFusionGenePairSet.find(tmpGenePair_2) != groundTruthFusionGenePairSet.end())
		{
			detected_true_genePair_set.insert(tmpGenePair_2);
			dt_tr_ofs << tmpStr << endl;
		}
		else
			dt_fl_ofs << tmpStr << endl;
	}
	dt_ifs.close();
	dt_tr_ofs.close();
	dt_fl_ofs.close();


	ofstream ms_tr_ofs(output_fusion_true_missed.c_str());
	for(set<string>::iterator tmpIter = groundTruthFusionGenePairSet.begin();
		tmpIter != groundTruthFusionGenePairSet.end(); tmpIter ++)
	{
		string tmpGtFusionGenePair = (*tmpIter);
		if(detected_true_genePair_set.find(tmpGtFusionGenePair) == detected_true_genePair_set.end()) // not detected
			ms_tr_ofs << tmpGtFusionGenePair.substr(0,15) << "\t" << tmpGtFusionGenePair.substr(15) << endl;
	}
	ms_tr_ofs.close();
	return 0;
}