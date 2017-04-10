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
		cout << "#2 inputFastaFile" << endl;
		cout << "#3 outputNonNodeFile" << endl;
		exit(1);
	}
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

	string inputFastaFile = argv[2];
	string outputNonNodeFile = argv[3];
	ifstream fa_ifs(inputFastaFile.c_str());
	ofstream nonNode_ofs(outputNonNodeFile.c_str());
	while(!fa_ifs.eof())
	{
		string tmpIdStr, tmpSeq;
		getline(fa_ifs, tmpIdStr);
		if(tmpIdStr == "")
			break;
		getline(fa_ifs, tmpSeq);
		int tabLoc_1 = tmpIdStr.find("_");
		int tabLoc_2 = tmpIdStr.find("_", tabLoc_1 + 1);
		int tabLoc_3 = tmpIdStr.find("_", tabLoc_2 + 1);
		int tabLoc_4 = tmpIdStr.find("_", tabLoc_3 + 1);
		int tabLoc_5 = tmpIdStr.find("_", tabLoc_4 + 1);
		int tabLoc_6 = tmpIdStr.find("_", tabLoc_5 + 1);
		string tmpIdStr_1 = tmpIdStr.substr(1, tabLoc_1 - 1);
		string tmpIdStr_2 = tmpIdStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string tmpIdStr_3 = tmpIdStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		string tmpIdStr_4 = tmpIdStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
		string tmpIdStr_5 = tmpIdStr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
		string tmpIdStr_6 = tmpIdStr.substr(tabLoc_5 + 1, tabLoc_6 - tabLoc_5 - 1);
		int tmpId_1 = atoi(tmpIdStr_1.c_str());
		int tmpId_2 = atoi(tmpIdStr_2.c_str());
		int tmpId_3 = atoi(tmpIdStr_3.c_str());
		int tmpId_4 = atoi(tmpIdStr_4.c_str());
		int tmpId_5 = atoi(tmpIdStr_5.c_str());
		int tmpId_6 = atoi(tmpIdStr_6.c_str());
		// cout << "tmpId_1: " << tmpId_1 << endl;
		// cout << "tmpId_2: " << tmpId_2 << endl;
		// cout << "tmpId_3: " << tmpId_3 << endl;
		// cout << "tmpId_4: " << tmpId_4 << endl;
		// cout << "tmpId_5: " << tmpId_5 << endl;
		// cout << "tmpId_6: " << tmpId_6 << endl;
		if(((taxoId2vecIndexMap.find(tmpId_1) == taxoId2vecIndexMap.end())&&(tmpId_1 > 0))
			||((taxoId2vecIndexMap.find(tmpId_2) == taxoId2vecIndexMap.end())&&(tmpId_2 > 0))
			||((taxoId2vecIndexMap.find(tmpId_3) == taxoId2vecIndexMap.end())&&(tmpId_3 > 0))
			||((taxoId2vecIndexMap.find(tmpId_4) == taxoId2vecIndexMap.end())&&(tmpId_4 > 0))
			||((taxoId2vecIndexMap.find(tmpId_5) == taxoId2vecIndexMap.end())&&(tmpId_5 > 0))
			||((taxoId2vecIndexMap.find(tmpId_6) == taxoId2vecIndexMap.end())&&(tmpId_6 > 0)))
		{
			nonNode_ofs << tmpIdStr;
			if((taxoId2vecIndexMap.find(tmpId_1) == taxoId2vecIndexMap.end())&&(tmpId_1 > 0))
				nonNode_ofs << "\t" << tmpId_1;
			if((taxoId2vecIndexMap.find(tmpId_2) == taxoId2vecIndexMap.end())&&(tmpId_2 > 0))
				nonNode_ofs << "\t" << tmpId_2;
			if((taxoId2vecIndexMap.find(tmpId_3) == taxoId2vecIndexMap.end())&&(tmpId_3 > 0))
				nonNode_ofs << "\t" << tmpId_3;
			if((taxoId2vecIndexMap.find(tmpId_4) == taxoId2vecIndexMap.end())&&(tmpId_4 > 0))
				nonNode_ofs << "\t" << tmpId_4;
			if((taxoId2vecIndexMap.find(tmpId_5) == taxoId2vecIndexMap.end())&&(tmpId_5 > 0))
				nonNode_ofs << "\t" << tmpId_5;
			if((taxoId2vecIndexMap.find(tmpId_6) == taxoId2vecIndexMap.end())&&(tmpId_6 > 0))
				nonNode_ofs << "\t" << tmpId_6;																			
			nonNode_ofs << endl;
		}
	}
	fa_ifs.close();
	nonNode_ofs.close();
	return 0;
}