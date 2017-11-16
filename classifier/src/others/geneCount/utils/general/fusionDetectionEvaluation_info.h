// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#ifndef FUSIONDETECTIONEVALUATION_INFO_H
#define FUSIONDETECTIONEVALUATION_INFO_H
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

//typedef map<string, int> GeneId2indexMap;

using namespace std;

class FusionDetectionEvaluation_Info
{
private:
	int gene_index_1;
	int gene_index_2;
	string gene_id_1;
	string gene_id_2;

	vector<string> readIdVec;
	vector<int> KmerCountVec_discriminative_1;
	vector<int> KmerCountVec_discriminative_2;
	vector<int> fusionSiteVec_discriminative;
	vector<int> windowDistanceVec_discriminative;

	vector<int> KmerCountVec_compatible_1;
	vector<int> KmerCountVec_compatible_2;
	vector<int> fusionSiteVec_compatible;
	vector<int> windowDistanceVec_compatible;
public:
	FusionDetectionEvaluation_Info()
	{
		//readCount = 0;
	}

	int return_gene_index_1()
	{
		return gene_index_1;
	}

	int return_gene_index_2()
	{
		return gene_index_2;
	}	

	string return_gene_id_1()
	{
		return gene_id_1;
	}

	string return_gene_id_2()
	{
		return gene_id_2;
	}	

	void initiate_1stRead(int tmp_gene_index_1, int tmp_gene_index_2, string& tmp_gene_id_1, string& tmp_gene_id_2, string& readId, 
		int KmerCount_discriminative_1, int KmerCount_discriminative_2, int fusionSite_discriminative, int windowDistance_discriminative,
		int KmerCount_compatible_1, int KmerCount_compatible_2, int fusionSite_compatible, int windowDistance_compatible)
	{
		gene_index_1 = tmp_gene_index_1;
		gene_index_2 = tmp_gene_index_2;
		gene_id_1 = tmp_gene_id_1;
		gene_id_2 = tmp_gene_id_2;

		this->addNewRead(readId, KmerCount_discriminative_1, KmerCount_discriminative_2, 
			fusionSite_discriminative, windowDistance_discriminative,
			KmerCount_compatible_1, KmerCount_compatible_2, 
			fusionSite_compatible, windowDistance_compatible);	
	}

	void addNewRead(string& readId, 
		int KmerCount_discriminative_1, int KmerCount_discriminative_2, int fusionSite_discriminative, int windowDistance_discriminative,
		int KmerCount_compatible_1, int KmerCount_compatible_2, int fusionSite_compatible, int windowDistance_compatible)
	{
		readIdVec.push_back(readId);

		KmerCountVec_discriminative_1.push_back(KmerCount_discriminative_1);
		KmerCountVec_discriminative_2.push_back(KmerCount_discriminative_2);
		fusionSiteVec_discriminative.push_back(fusionSite_discriminative);
		windowDistanceVec_discriminative.push_back(windowDistance_discriminative);

		KmerCountVec_compatible_1.push_back(KmerCount_compatible_1);
		KmerCountVec_compatible_2.push_back(KmerCount_compatible_2);
		fusionSiteVec_compatible.push_back(fusionSite_compatible);
		windowDistanceVec_compatible.push_back(windowDistance_compatible);
	}

	int get_readCount()
	{
		return readIdVec.size();
	}

	void get_fusionSite_discriminative_mean_sd(int& tmpMean, double& tmpSd)
	{
		int tmpSum = 0;
		for(int tmp = 0; tmp < fusionSiteVec_discriminative.size(); tmp++)
			tmpSum += fusionSiteVec_discriminative[tmp];
		tmpMean = (double)tmpSum/(double)fusionSiteVec_discriminative.size();
		int tmpSquareSum = 0;
		for(int tmp = 0; tmp < fusionSiteVec_discriminative.size(); tmp++)
			tmpSquareSum += pow((fusionSiteVec_discriminative[tmp] - tmpMean), 2);
		tmpSd = (double)sqrt((double)tmpSquareSum)/(double)fusionSiteVec_discriminative.size();		
	}

	void get_windowDistanceVec_discriminative_mean_sd(int& tmpMean, double& tmpSd)
	{
		int tmpSum = 0;
		for(int tmp = 0; tmp < windowDistanceVec_discriminative.size(); tmp++)
			tmpSum += windowDistanceVec_discriminative[tmp];
		tmpMean = (double)tmpSum/(double)windowDistanceVec_discriminative.size();
		int tmpSquareSum = 0;
		for(int tmp = 0; tmp < windowDistanceVec_discriminative.size(); tmp++)
			tmpSquareSum += pow((windowDistanceVec_discriminative[tmp] - tmpMean), 2);
		tmpSd = (double)sqrt((double)tmpSquareSum)/(double)windowDistanceVec_discriminative.size();		
	}	

	void get_KmerCount_discriminative_1_mean_sd(int& tmpMean, double& tmpSd)
	{
		int tmpSum = 0;
		for(int tmp = 0; tmp < KmerCountVec_discriminative_1.size(); tmp++)
			tmpSum += KmerCountVec_discriminative_1[tmp];
		tmpMean = (double)tmpSum/(double)KmerCountVec_discriminative_1.size();
		int tmpSquareSum = 0;
		for(int tmp = 0; tmp < KmerCountVec_discriminative_1.size(); tmp++)
			tmpSquareSum += pow((KmerCountVec_discriminative_1[tmp] - tmpMean), 2);
		tmpSd = (double)sqrt((double)tmpSquareSum)/(double)KmerCountVec_discriminative_1.size();
	}

	void get_KmerCount_discriminative_2_mean_sd(int& tmpMean, double& tmpSd)
	{
		int tmpSum = 0;
		for(int tmp = 0; tmp < KmerCountVec_discriminative_2.size(); tmp++)
			tmpSum += KmerCountVec_discriminative_2[tmp];
		tmpMean = tmpSum/KmerCountVec_discriminative_2.size();
		
		int tmpSquareSum = 0;
		for(int tmp = 0; tmp < KmerCountVec_discriminative_2.size(); tmp++)
			tmpSquareSum += pow((KmerCountVec_discriminative_2[tmp] - tmpMean), 2);
		tmpSd = (double)sqrt((double)tmpSquareSum)/(double)KmerCountVec_discriminative_2.size();	
	}

	void get_fusionSite_compatible_mean_sd(int& tmpMean, double& tmpSd)
	{
		int tmpSum = 0;
		for(int tmp = 0; tmp < fusionSiteVec_compatible.size(); tmp++)
			tmpSum += fusionSiteVec_compatible[tmp];
		tmpMean = (double)tmpSum/(double)fusionSiteVec_compatible.size();
		int tmpSquareSum = 0;
		for(int tmp = 0; tmp < fusionSiteVec_compatible.size(); tmp++)
			tmpSquareSum += pow((fusionSiteVec_compatible[tmp] - tmpMean), 2);
		tmpSd = (double)sqrt((double)tmpSquareSum)/(double)fusionSiteVec_compatible.size();		
	}

	void get_windowDistanceVec_compatible_mean_sd(int& tmpMean, double& tmpSd)
	{
		int tmpSum = 0;
		for(int tmp = 0; tmp < windowDistanceVec_compatible.size(); tmp++)
			tmpSum += windowDistanceVec_compatible[tmp];
		tmpMean = (double)tmpSum/(double)windowDistanceVec_compatible.size();
		int tmpSquareSum = 0;
		for(int tmp = 0; tmp < windowDistanceVec_compatible.size(); tmp++)
			tmpSquareSum += pow((windowDistanceVec_compatible[tmp] - tmpMean), 2);
		tmpSd = (double)sqrt((double)tmpSquareSum)/(double)windowDistanceVec_compatible.size();		
	}	

	void get_KmerCount_compatible_1_mean_sd(int& tmpMean, double& tmpSd)
	{
		int tmpSum = 0;
		for(int tmp = 0; tmp < KmerCountVec_compatible_1.size(); tmp++)
			tmpSum += KmerCountVec_compatible_1[tmp];
		tmpMean = (double)tmpSum/(double)KmerCountVec_compatible_1.size();
		int tmpSquareSum = 0;
		for(int tmp = 0; tmp < KmerCountVec_compatible_1.size(); tmp++)
			tmpSquareSum += pow((KmerCountVec_compatible_1[tmp] - tmpMean), 2);
		tmpSd = (double)sqrt((double)tmpSquareSum)/(double)KmerCountVec_compatible_1.size();
	}

	void get_KmerCount_compatible_2_mean_sd(int& tmpMean, double& tmpSd)
	{
		int tmpSum = 0;
		for(int tmp = 0; tmp < KmerCountVec_compatible_2.size(); tmp++)
			tmpSum += KmerCountVec_compatible_2[tmp];
		tmpMean = tmpSum/KmerCountVec_compatible_2.size();
		
		int tmpSquareSum = 0;
		for(int tmp = 0; tmp < KmerCountVec_compatible_2.size(); tmp++)
			tmpSquareSum += pow((KmerCountVec_compatible_2[tmp] - tmpMean), 2);
		tmpSd = (double)sqrt((double)tmpSquareSum)/(double)KmerCountVec_compatible_2.size();	
	}

	void summmarize(int& readCount, 
		int& tmpMean_KmerCount_discriminative_1, double& tmpSd_KmerCount_discriminative_1,
		int& tmpMean_KmerCount_discriminative_2, double& tmpSd_KmerCount_discriminative_2,
		int& tmpMean_fusionSite_discriminative, double& tmpSd_fusionSite_discriminative,
		int& tmpMean_windowDistanceVec_discriminative, double& tmpSd_windowDistanceVec_discriminative,
		int& tmpMean_KmerCount_compatible_1, double& tmpSd_KmerCount_compatible_1,
		int& tmpMean_KmerCount_compatible_2, double& tmpSd_KmerCount_compatible_2,
		int& tmpMean_fusionSite_compatible, double& tmpSd_fusionSite_compatible,
		int& tmpMean_windowDistanceVec_compatible, double& tmpSd_windowDistanceVec_compatible)
	{
		readCount = this->get_readCount();

		this->get_KmerCount_discriminative_1_mean_sd(tmpMean_KmerCount_discriminative_1, tmpSd_KmerCount_discriminative_1);
		this->get_KmerCount_discriminative_2_mean_sd(tmpMean_KmerCount_discriminative_2, tmpSd_KmerCount_discriminative_2);
		this->get_fusionSite_discriminative_mean_sd(tmpMean_fusionSite_discriminative, tmpSd_fusionSite_discriminative);
		this->get_windowDistanceVec_discriminative_mean_sd(tmpMean_windowDistanceVec_discriminative, tmpSd_windowDistanceVec_discriminative);

		this->get_KmerCount_compatible_1_mean_sd(tmpMean_KmerCount_compatible_1, tmpSd_KmerCount_compatible_1);
		this->get_KmerCount_compatible_2_mean_sd(tmpMean_KmerCount_compatible_2, tmpSd_KmerCount_compatible_2);
		this->get_fusionSite_compatible_mean_sd(tmpMean_fusionSite_compatible, tmpSd_fusionSite_compatible);
		this->get_windowDistanceVec_compatible_mean_sd(tmpMean_windowDistanceVec_compatible, tmpSd_windowDistanceVec_compatible);
	}
};

class FusionDetectionEvaluation_Vec_Info
{
private:
	map<string, int> geneId2indexMap;
	vector<string> geneIdVec;
	map< int, map<int,int> > geneIndexPair2fusionIndexMap;

	vector<FusionDetectionEvaluation_Info> fusionInfoVec;
public:
	FusionDetectionEvaluation_Vec_Info()
	{}

	void initiate_geneIdList(string& geneIdListFile)
	{
		ifstream geneIdList_ifs(geneIdListFile.c_str());
		int tmpIndex = 0;
		while(!geneIdList_ifs.eof())
		{
			string tmpStr;
			getline(geneIdList_ifs, tmpStr);
			if(tmpStr == "")
				break;
			tmpIndex ++;
			geneId2indexMap.insert(pair<string,int>(tmpStr, tmpIndex));
			geneIdVec.push_back(tmpStr);
		}
		geneIdList_ifs.close();
	}

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

	void load_fusionReadInfoFile(string& tmpFusionReadInfoFile, int tmpPairReadLengthSum)
	{
		ifstream results_ifs(tmpFusionReadInfoFile.c_str());
		while(!results_ifs.eof())
		{
			string tmpStr; 
			getline(results_ifs, tmpStr);
			if(tmpStr == "")
				break;
			this->update_fusionReadInfo(tmpStr, tmpPairReadLengthSum);
		}
		results_ifs.close();
	}

	void update_fusionReadInfo(string& tmpFusionReadInfoStr, int tmpPairReadLengthSum)
	{
		vector<string> tmpFieldVec;
		this->parseStr2fieldVec_tab(tmpFieldVec, tmpFusionReadInfoStr);
		string tmpReadId = tmpFieldVec[0];
		//cout << "tmpReadId: " << tmpReadId << endl;
		map<string, int>::iterator tmpIter_1 = geneId2indexMap.find(tmpFieldVec[1]);
		map<string, int>::iterator tmpIter_2 = geneId2indexMap.find(tmpFieldVec[2]);
		if((tmpIter_1 == geneId2indexMap.end())||(tmpIter_2 == geneId2indexMap.end()))
		{
			cout << "tmpGene not found: " << tmpFieldVec[1] << "\t" << tmpFieldVec[2] << endl;
			exit(1);
		}

		string tmp_gene_id_1, tmp_gene_id_2;
		int tmp_gene_index_1, tmp_gene_index_2,
			tmp_KmerCount_discriminative_1, tmp_KmerCount_discriminative_2, 
			tmp_fusionSite_discriminative, tmp_windowDistance_discriminative,
			tmp_KmerCount_compatible_1, tmp_KmerCount_compatible_2,
			tmp_fusionSite_compatible, tmp_windowDistance_compatible;
		
		if((tmpIter_1->second) < (tmpIter_2->second))
		{
			tmp_gene_index_1 = tmpIter_1->second;
			tmp_gene_index_2 = tmpIter_2->second;

			tmp_gene_id_1 = tmpFieldVec[1];
			tmp_gene_id_2 = tmpFieldVec[2];

			tmp_KmerCount_discriminative_1 = atoi(tmpFieldVec[3].c_str());
			tmp_KmerCount_discriminative_2 = atoi(tmpFieldVec[4].c_str());
			tmp_fusionSite_discriminative = atoi(tmpFieldVec[5].c_str());
			tmp_windowDistance_discriminative = atoi(tmpFieldVec[6].c_str());

			tmp_KmerCount_compatible_1 = atoi(tmpFieldVec[7].c_str());
			tmp_KmerCount_compatible_2 = atoi(tmpFieldVec[8].c_str());
			tmp_fusionSite_compatible = atoi(tmpFieldVec[9].c_str());
			tmp_windowDistance_compatible = atoi(tmpFieldVec[10].c_str());
		}
		else
		{
			tmp_gene_index_1 = tmpIter_2->second;
			tmp_gene_index_2 = tmpIter_1->second;

			tmp_gene_id_1 = tmpFieldVec[2];
			tmp_gene_id_2 = tmpFieldVec[1];

			tmp_KmerCount_discriminative_1 = atoi(tmpFieldVec[4].c_str());
			tmp_KmerCount_discriminative_2 = atoi(tmpFieldVec[3].c_str());
			tmp_fusionSite_discriminative = tmpPairReadLengthSum - atoi(tmpFieldVec[5].c_str());
			tmp_windowDistance_discriminative = atoi(tmpFieldVec[6].c_str());

			tmp_KmerCount_compatible_1 = atoi(tmpFieldVec[8].c_str());
			tmp_KmerCount_compatible_2 = atoi(tmpFieldVec[7].c_str());
			tmp_fusionSite_compatible = tmpPairReadLengthSum - atoi(tmpFieldVec[9].c_str());
			tmp_windowDistance_compatible = atoi(tmpFieldVec[10].c_str());
		}
		//cout << "tmp_gene_index_1: " << tmp_gene_index_1 << endl << "tmp_gene_index_2: " << tmp_gene_index_2 << endl;
		//cout << "tmp_gene_id_1: " << tmp_gene_id_1 << endl << "tmp_gene_id_2: " << tmp_gene_id_2 << endl;
		map< int, map<int,int> >::iterator tmpFusionId2indexMapIter = geneIndexPair2fusionIndexMap.find(tmp_gene_index_1);
		if(tmpFusionId2indexMapIter != geneIndexPair2fusionIndexMap.end()) // gene_1 found
		{
			map<int,int>::iterator tmpMapIter = (tmpFusionId2indexMapIter->second).find(tmp_gene_index_2);
			if(tmpMapIter != (tmpFusionId2indexMapIter->second).end()) // gene_1 found, gene_2 found
			{
				int tmpFusionIndex = tmpMapIter->second;
				fusionInfoVec[tmpFusionIndex].addNewRead(tmpReadId, tmp_KmerCount_discriminative_1, tmp_KmerCount_discriminative_2, 
					tmp_fusionSite_discriminative, tmp_windowDistance_discriminative, tmp_KmerCount_compatible_1, 
					tmp_KmerCount_compatible_2, tmp_fusionSite_compatible, tmp_windowDistance_compatible);
			}
			else // gene_1 found, gene_2 not found
			{
				int tmp_new_fusionVec_index = fusionInfoVec.size();
				(tmpFusionId2indexMapIter->second).insert(pair<int,int>(tmp_gene_index_2, tmp_new_fusionVec_index));
				FusionDetectionEvaluation_Info tmpFusionInfo;
				tmpFusionInfo.initiate_1stRead(tmp_gene_index_1, tmp_gene_index_2, tmp_gene_id_1, tmp_gene_id_2, tmpReadId, 
					tmp_KmerCount_discriminative_1, tmp_KmerCount_discriminative_2, tmp_fusionSite_discriminative, tmp_windowDistance_discriminative,
					tmp_KmerCount_compatible_1, tmp_KmerCount_compatible_2, tmp_fusionSite_compatible, tmp_windowDistance_compatible);
				fusionInfoVec.push_back(tmpFusionInfo);				
			}
		}
		else
		{
			map<int,int> tmpMap;
			int tmp_new_fusionVec_index = fusionInfoVec.size();
			tmpMap.insert(pair<int,int>(tmp_gene_index_2, tmp_new_fusionVec_index));
			geneIndexPair2fusionIndexMap.insert(pair<int, map<int,int> >(tmp_gene_index_1, tmpMap));
			FusionDetectionEvaluation_Info tmpFusionInfo;
			tmpFusionInfo.initiate_1stRead(tmp_gene_index_1, tmp_gene_index_2, tmp_gene_id_1, tmp_gene_id_2, tmpReadId, 
				tmp_KmerCount_discriminative_1, tmp_KmerCount_discriminative_2, tmp_fusionSite_discriminative, tmp_windowDistance_discriminative,
				tmp_KmerCount_compatible_1, tmp_KmerCount_compatible_2, tmp_fusionSite_compatible, tmp_windowDistance_compatible);
			fusionInfoVec.push_back(tmpFusionInfo);
		}
	}

	void print_fusion_summary(string& tmpFusionSummaryFile_discriminative, string& tmpFusionSummaryFile_compatible)
	{
		ofstream sum_ofs_discriminative(tmpFusionSummaryFile_discriminative.c_str());
		ofstream sum_ofs_compatible(tmpFusionSummaryFile_compatible.c_str());
		for(int tmp = 0; tmp < fusionInfoVec.size(); tmp++)
		{
			int tmpReadCount;
			int tmpMean_KmerCount_discriminative_1; double tmpSd_KmerCount_discriminative_1;
			int tmpMean_KmerCount_discriminative_2; double tmpSd_KmerCount_discriminative_2;
			int tmpMean_fusionSite_discriminative; double tmpSd_fusionSite_discriminative;
			int tmpMean_windowDistanceVec_discriminative; double tmpSd_windowDistanceVec_discriminative;
			int tmpMean_KmerCount_compatible_1; double tmpSd_KmerCount_compatible_1;
			int tmpMean_KmerCount_compatible_2; double tmpSd_KmerCount_compatible_2;
			int tmpMean_fusionSite_compatible; double tmpSd_fusionSite_compatible;
			int tmpMean_windowDistanceVec_compatible; double tmpSd_windowDistanceVec_compatible;
			fusionInfoVec[tmp].summmarize(tmpReadCount,
				tmpMean_KmerCount_discriminative_1, tmpSd_KmerCount_discriminative_1,
				tmpMean_KmerCount_discriminative_2, tmpSd_KmerCount_discriminative_2,
				tmpMean_fusionSite_discriminative, tmpSd_fusionSite_discriminative,
				tmpMean_windowDistanceVec_discriminative, tmpSd_windowDistanceVec_discriminative,
				tmpMean_KmerCount_compatible_1, tmpSd_KmerCount_compatible_1,
				tmpMean_KmerCount_compatible_2, tmpSd_KmerCount_compatible_2,
				tmpMean_fusionSite_compatible, tmpSd_fusionSite_compatible,
				tmpMean_windowDistanceVec_compatible, tmpSd_windowDistanceVec_compatible);
			int tmp_gene_index_1 = fusionInfoVec[tmp].return_gene_index_1();
			int tmp_gene_index_2 = fusionInfoVec[tmp].return_gene_index_2();
			string tmp_gene_id_1 = fusionInfoVec[tmp].return_gene_id_1();
			string tmp_gene_id_2 = fusionInfoVec[tmp].return_gene_id_2();
			sum_ofs_discriminative << tmp_gene_index_1 << "\t" << tmp_gene_index_2 << "\t"
				<< tmp_gene_id_1 << "\t" << tmp_gene_id_2 << "\t" << tmpReadCount << "\t" 
				<< tmpMean_KmerCount_discriminative_1 << "-" << tmpSd_KmerCount_discriminative_1 << "\t"
				<< tmpMean_KmerCount_discriminative_2 << "-" << tmpSd_KmerCount_discriminative_2 << "\t"
				<< tmpMean_fusionSite_discriminative << "-" << tmpSd_fusionSite_discriminative << "\t"
				<< tmpMean_windowDistanceVec_discriminative << "-" << tmpSd_windowDistanceVec_discriminative << endl;
			sum_ofs_compatible << tmp_gene_index_1 << "\t" << tmp_gene_index_2 << "\t"
				<< tmp_gene_id_1 << "\t" << tmp_gene_id_2 << "\t" << tmpReadCount << "\t"
				<< tmpMean_KmerCount_compatible_1 << "-" << tmpSd_KmerCount_compatible_1 << "\t"
				<< tmpMean_KmerCount_compatible_2 << "-" << tmpSd_KmerCount_compatible_2 << "\t"
				<< tmpMean_fusionSite_compatible << "-" << tmpSd_fusionSite_compatible << "\t"
				<< tmpMean_windowDistanceVec_compatible << "-" << tmpSd_windowDistanceVec_compatible << endl;				
		}
		sum_ofs_discriminative.close();
		sum_ofs_compatible.close();
	}
};
#endif