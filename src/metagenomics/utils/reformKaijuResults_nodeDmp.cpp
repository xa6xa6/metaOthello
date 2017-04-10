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

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputNodeDmpFile" << endl;
		cout << "#2 inputKaijuResults_raw" << endl;
		cout << "#3 outputKaijuResults_reformed" << endl;
		exit(1);
	}
	string inputKaijuResults_raw = argv[2];
	string outputKaijuResults_reformed = argv[3];

	string input_node_dmp_file = argv[1];
	map<int, pair<int,string> > nodeDmpMap;
	cout << "start to load node dmp file" << endl;
	ifstream node_dmp_ifs(input_node_dmp_file.c_str());
	while(!node_dmp_ifs.eof())
	{
		string tmpStr;
		getline(node_dmp_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		if(tabLoc_1 == string::npos)
			continue;
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		if(tabLoc_2 == string::npos)
			continue;
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		if(tabLoc_3 == string::npos)
			continue;
		int tabLoc_4 = tmpStr.find("\t", tabLoc_3 + 1);
		if(tabLoc_4 == string::npos)
			continue;
		int tabLoc_5 = tmpStr.find("\t", tabLoc_4 + 1);
		if(tabLoc_5 == string::npos)
			continue;
		string tmp_id_1_str = tmpStr.substr(0, tabLoc_1);
		int tmp_id_1 = atoi(tmp_id_1_str.c_str());
		string tmp_id_2_str = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		int tmp_id_2 = atoi(tmp_id_2_str.c_str());
		string tmp_rank_str = tmpStr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
		nodeDmpMap.insert(pair<int, pair<int,string> >(tmp_id_1, pair<int,string>(tmp_id_2, tmp_rank_str)));
	}
	node_dmp_ifs.close();

	cout << "start to generate taxoIdRankPairVecVec and taxoId2vecIndexMap" << endl;
	map<int, int> taxoId2vecIndexMap;
	vector< vector< pair<int,string> > > taxoIdRankPairVecVec;
	int tmpVecIndex = 0;
	for(map<int, pair<int,string> >::iterator tmpIter = nodeDmpMap.begin(); tmpIter != nodeDmpMap.end(); tmpIter++)
	{
		vector< pair<int,string> > tmpTaxoIdRankPairVec;
		int tmpIter_id_1 = tmpIter->first;
		int tmpIter_id_2 = (tmpIter->second).first;
		string tmpIter_id_1_rank = (tmpIter->second).second;
		taxoId2vecIndexMap.insert(pair<int,int>(tmpIter_id_1, tmpVecIndex));
		tmpVecIndex ++;
		tmpTaxoIdRankPairVec.push_back(pair<int,string>(tmpIter_id_1, tmpIter_id_1_rank));
		while(nodeDmpMap.find(tmpIter_id_2) != nodeDmpMap.end())
		{
			map<int, pair<int,string> >::iterator tmpIter_2 = nodeDmpMap.find(tmpIter_id_2);
			if(tmpIter_2 == nodeDmpMap.end()) // not found
			{
				cout << "error! (nodeDmpMap.find(tmpIter_id_2) != nodeDmpMap.end()) but (tmpIter_2 == nodeDmpMap.end())" << endl;
				exit(1);
			}
			// found
			tmpIter_id_1 = tmpIter_2->first;
			tmpIter_id_2 = (tmpIter_2->second).first;
			tmpIter_id_1_rank = (tmpIter_2->second).second;
			tmpTaxoIdRankPairVec.push_back(pair<int,string>(tmpIter_id_1, tmpIter_id_1_rank));		
			if(tmpIter_id_1 == tmpIter_id_2)
				break;		
		}
		taxoIdRankPairVecVec.push_back(tmpTaxoIdRankPairVec);
	}
	cout << "end of generating taxoIdRankPairVecVec and taxoId2vecIndexMap" << endl;
	ifstream raw_ifs(inputKaijuResults_raw.c_str());
	ofstream reformed_ofs(outputKaijuResults_reformed.c_str());
	int tmpReadNO = 0;
	while(!raw_ifs.eof())
	{
		string tmpStr;
		getline(raw_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		string tmpRead_mappedOrNot_str = tmpStr.substr(0, tabLoc_1);
		string tmpRead_name = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string tmpTaxoIdStr = tmpStr.substr(tabLoc_2 + 1);
		int tmpTaxoId = atoi(tmpTaxoIdStr.c_str());
		//cout << "tmpTaxoId: " << tmpTaxoId << endl;
		tmpReadNO ++;
		// if(tmpReadNO > 10)
		// 	break;
		if(tmpTaxoId < 0)
		{
			reformed_ofs << ">" << tmpRead_name << "\t-1\t-1\t-1\t-1\t-1\t-1" << endl;
			continue;
		}

		if(tmpRead_mappedOrNot_str == "C")
		{}
		else if(tmpRead_mappedOrNot_str == "U")
		{
			reformed_ofs << ">" << tmpRead_name << "\t-1\t-1\t-1\t-1\t-1\t-1" << endl;
			continue;
		}
		else
		{
			cout << "invalid tmpRead_mappedOrNot_str: " << tmpRead_mappedOrNot_str << endl;
			exit(1);
		}
		reformed_ofs << ">" << tmpRead_name << "\t";

		map<int,int>::iterator tmpTaxoId2vecIndexMapIter = taxoId2vecIndexMap.find(tmpTaxoId);
		if(tmpTaxoId2vecIndexMapIter == taxoId2vecIndexMap.end())
		{
			reformed_ofs << "-1\t-1\t-1\t-1\t-1\t-1" << endl;
			continue;
		}			
		else
		{
			int tmpFoundVecIndex = tmpTaxoId2vecIndexMapIter->second;
			int tmpSpeciesId = -1, tmpGenusId = -1, tmpFamilyId = -1, tmpOrderId = -1, tmpClassId = -1, tmpPhylumId = -1;
			int tmpTaxoVecSize = taxoIdRankPairVecVec[tmpFoundVecIndex].size();
			for(int tmp2 = 0; tmp2 < tmpTaxoVecSize; tmp2 ++)
			{
				if(((taxoIdRankPairVecVec[tmpFoundVecIndex])[tmp2]).second == "species")
					tmpSpeciesId = ((taxoIdRankPairVecVec[tmpFoundVecIndex])[tmp2]).first;
				else if(((taxoIdRankPairVecVec[tmpFoundVecIndex])[tmp2]).second == "genus")
					tmpGenusId = ((taxoIdRankPairVecVec[tmpFoundVecIndex])[tmp2]).first;
				else if(((taxoIdRankPairVecVec[tmpFoundVecIndex])[tmp2]).second == "family")
					tmpFamilyId = ((taxoIdRankPairVecVec[tmpFoundVecIndex])[tmp2]).first;
				else if(((taxoIdRankPairVecVec[tmpFoundVecIndex])[tmp2]).second == "order")
					tmpOrderId = ((taxoIdRankPairVecVec[tmpFoundVecIndex])[tmp2]).first;
				else if(((taxoIdRankPairVecVec[tmpFoundVecIndex])[tmp2]).second == "class")
					tmpClassId = ((taxoIdRankPairVecVec[tmpFoundVecIndex])[tmp2]).first;
				else if(((taxoIdRankPairVecVec[tmpFoundVecIndex])[tmp2]).second == "phylum")
					tmpPhylumId = ((taxoIdRankPairVecVec[tmpFoundVecIndex])[tmp2]).first;
				else
				{}																								
			}
			reformed_ofs << tmpSpeciesId << "\t" << tmpGenusId << "\t" << tmpFamilyId << "\t" 
				<< tmpOrderId << "\t" << tmpClassId << "\t" << tmpPhylumId << endl;
			continue;		
		}
	}
	raw_ifs.close();
	reformed_ofs.close();
	return 0;
}