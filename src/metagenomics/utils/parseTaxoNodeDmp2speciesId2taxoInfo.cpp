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

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 input_node_dmp_file" << endl;
		cout << "#2 output_speciesId2taxoInfo_complete_file" << endl;
		exit(1);
	}
	string input_node_dmp_file = argv[1];
	string output_speciesId2taxoInfo_complete_file = argv[2];
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
	cout << "start to generate taxoIdRankPairVecVec" << endl;
	vector< vector< pair<int,string> > > taxoIdRankPairVecVec;
	for(map<int, pair<int,string> >::iterator tmpIter = nodeDmpMap.begin(); tmpIter != nodeDmpMap.end(); tmpIter++)
	{
		vector< pair<int,string> > tmpTaxoIdRankPairVec;
		int tmpIter_id_1 = tmpIter->first;
		int tmpIter_id_2 = (tmpIter->second).first;
		string tmpIter_id_1_rank = (tmpIter->second).second;
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
	cout << "start to print taxo id rank pair info" << endl;
	ofstream speciesId2taxoInfo_complete_ofs(output_speciesId2taxoInfo_complete_file.c_str());
	for(int tmp = 0; tmp < taxoIdRankPairVecVec.size(); tmp++)
	{
		int tmpTaxoVecSize = taxoIdRankPairVecVec[tmp].size();
		if(tmpTaxoVecSize < 1)
		{
			cout << "error! tmpTaxoVecSize < 1" << endl;
			exit(1);
		}
		speciesId2taxoInfo_complete_ofs << ((taxoIdRankPairVecVec[tmp])[0]).first << "\t" << ((taxoIdRankPairVecVec[tmp])[0]).second;
		if(tmpTaxoVecSize > 1)
		{	
			for(int tmp2 = 1; tmp2 < tmpTaxoVecSize; tmp2++)
				speciesId2taxoInfo_complete_ofs << "\t" << (taxoIdRankPairVecVec[tmp])[tmp2].first 
					<< "\t" << (taxoIdRankPairVecVec[tmp])[tmp2].second;
		}
		speciesId2taxoInfo_complete_ofs << endl;
	}
	speciesId2taxoInfo_complete_ofs.close();
	return 0;
}