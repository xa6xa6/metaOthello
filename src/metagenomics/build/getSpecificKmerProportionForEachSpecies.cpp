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
#include "../../mps3Lib/read_block_test.h"
#include "../../mps3Lib/otherFunc.h"
#include "../../mps3Lib/index_info.h"
using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "#0 Executable" << endl;
		// /scratch/lcph222/Xinan/lothelloClassifier/metagenomics/varyLength/20mer/merged_specific_repetitive/merged.Kmer.setId.distribution
		cout << "#1 inputMergedKmerSetIdDistributionFIle" << endl;
		cout << "#2 inputJfHistoDir" << endl;
		// /scratch/lcph222/Xinan/lothelloClassifier/metagenomics/varyLength/20mer/merged_specific_repetitive/taxo_info/species.txt
		cout << "#3 SpeciesInfoFile" << endl;
		cout << "#4 outputSpeciesDiscriminativeTotalKmerCountFile" << endl;
		exit(1);
	}
	string inputMergedKmerSetIdDistributionFIle = argv[1];
	string inputJfHistoDir = argv[2];
	string SpeciesInfoFile = argv[3];
	string outputSpeciesDiscriminativeTotalKmerCountFile = argv[4];

	vector<unsigned long long> KmerFreqVec_specific;
	ifstream mergedKmerDistribution_ifs(inputMergedKmerSetIdDistributionFIle.c_str());
	string tmp1stLineInMergedKmerDistributionFile;
	getline(mergedKmerDistribution_ifs, tmp1stLineInMergedKmerDistributionFile);
	while(!mergedKmerDistribution_ifs.eof())
	{
		string tmpStr;
		getline(mergedKmerDistribution_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpFreqStr = tmpStr.substr(tabLoc + 1);
		//string tmpIdStr = tmpStr.substr(0, tabLoc);
		//int tmpId = atoi(tmpIdStr.c_str());
		unsigned long long tmpFreq = atoll(tmpFreqStr.c_str());
		KmerFreqVec_specific.push_back(tmpFreq);
	}
	mergedKmerDistribution_ifs.close();
	cout << "KmerFreqVec_specific.size(): " << KmerFreqVec_specific.size() << endl;

	vector<int> speciesIdVec;
	ifstream speciesInfo_ifs(SpeciesInfoFile.c_str());
	string tmp1stLineInSpeciesInfoFile;
	getline(speciesInfo_ifs, tmp1stLineInSpeciesInfoFile);
	while(!speciesInfo_ifs.eof())
	{
		string tmpStr;
		getline(speciesInfo_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		string tmpSpeciesIdStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		int tmpSpeciesIdInt = atoi(tmpSpeciesIdStr.c_str());
		speciesIdVec.push_back(tmpSpeciesIdInt);
	}
	speciesInfo_ifs.close();
	cout << "speciesIdVec.size(): " << speciesIdVec.size() << endl;

	vector<unsigned long long> KmerFreqVec_total;
	inputJfHistoDir += "/";
	for(int tmp = 0; tmp < speciesIdVec.size(); tmp++)
	{
		vector<unsigned long long> tmpKmerCountVec;
		string tmpFile = inputJfHistoDir + int_to_str(speciesIdVec[tmp]) + ".jfHisto";
		ifstream tmp_ifs(tmpFile.c_str());
		while(!tmp_ifs.eof())
		{
			string tmpStr;
			getline(tmp_ifs, tmpStr);
			if(tmpStr == "")
				break;
			int blankLoc = tmpStr.find(" ");
			string tmpCountStr = tmpStr.substr(blankLoc + 1);
			unsigned long long tmpCountInt = atoll(tmpCountStr.c_str());
			tmpKmerCountVec.push_back(tmpCountInt);
		}
		tmp_ifs.close();
		unsigned long long tmpKmerCount_total = 0;
		for(int tmp2 = 0; tmp2 < tmpKmerCountVec.size(); tmp2++)
			tmpKmerCount_total += tmpKmerCountVec[tmp2];
		KmerFreqVec_total.push_back(tmpKmerCount_total);
	}
	cout << "KmerFreqVec_total.size(): " << KmerFreqVec_total.size() << endl;
	ofstream output_ofs(outputSpeciesDiscriminativeTotalKmerCountFile.c_str());
	for(int tmp = 0; tmp < speciesIdVec.size(); tmp++)
	{
		unsigned long long tmpKmerCount_specific = KmerFreqVec_specific[tmp];
		unsigned long long tmpKmerCount_total = KmerFreqVec_total[tmp];
		double tmpPerc = (double)tmpKmerCount_specific/(double)tmpKmerCount_total;
		int tmpSpeciesId = speciesIdVec[tmp];
		output_ofs << tmp << "\t" << tmpSpeciesId << "\t" << tmpKmerCount_specific 
			<< "\t" << tmpKmerCount_total << "\t" << tmpPerc << endl;
	}
	output_ofs.close();
	return 0;
}