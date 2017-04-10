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
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputSpeciesId2taxoIdFile" << endl;
		cout << "#2 inputKaijuResults_raw" << endl;
		cout << "#3 outputKaijuResults_reformed" << endl;
		exit(1);
	}
	string inputSpeciesId2taxoIdFile = argv[1];
	string inputKaijuResults_raw = argv[2];
	string outputKaijuResults_reformed = argv[3];

	vector<int> speciesIdVec;
	vector<int> genusIdVec;
	vector<int> familyIdVec;
	vector<int> orderIdVec;
	vector<int> classIdVec;
	vector<int> phylumIdVec;
	map<int,int> speciesId2vecIndexMap;
	map<int,int> genusId2vecIndexMap;
	map<int,int> familyId2vecIndexMap;
	map<int,int> orderId2vecIndexMap;
	map<int,int> classId2vecIndexMap;
	map<int,int> phylumId2vecIndexMap;

	int tmpVecIndex = 0;
	ifstream speciesId2taxoId_ifs(inputSpeciesId2taxoIdFile.c_str());
	while(!speciesId2taxoId_ifs.eof())
	{
		string tmpStr;
		getline(speciesId2taxoId_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		int tabLoc_4 = tmpStr.find("\t", tabLoc_3 + 1);
		int tabLoc_5 = tmpStr.find("\t", tabLoc_4 + 1);
		string tmpSpeciesIdStr = tmpStr.substr(0, tabLoc_1);
		string tmpGenusIdStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string tmpFamilyIdStr = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		string tmpOrderIdStr = tmpStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
		string tmpClassIdStr = tmpStr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
		string tmpPhylumIdStr = tmpStr.substr(tabLoc_5 + 1);
		int tmpSpeciesId = atoi(tmpSpeciesIdStr.c_str());
		int tmpGenusId = atoi(tmpGenusIdStr.c_str());
		int tmpFamilyId = atoi(tmpFamilyIdStr.c_str());
		int tmpOrderId = atoi(tmpOrderIdStr.c_str());
		int tmpClassId = atoi(tmpClassIdStr.c_str());
		int tmpPhylumId = atoi(tmpPhylumIdStr.c_str());
		speciesIdVec.push_back(tmpSpeciesId);
		genusIdVec.push_back(tmpGenusId);
		familyIdVec.push_back(tmpFamilyId);
		orderIdVec.push_back(tmpOrderId);
		classIdVec.push_back(tmpClassId);
		phylumIdVec.push_back(tmpPhylumId);
		// cout << "tmpSpeciesId: " << tmpSpeciesId << endl;
		// cout << "tmpGenusId: " << tmpGenusId << endl;
		// cout << "tmpFamilyId: " << tmpFamilyId << endl;
		// cout << "tmpOrderId: " << tmpOrderId << endl;
		// cout << "tmpClassId: " << tmpClassId << endl;
		// cout << "tmpPhylumId: " << tmpPhylumId << endl;
		if((speciesId2vecIndexMap.find(tmpSpeciesId) == speciesId2vecIndexMap.end())&&(tmpSpeciesId > 0))
			speciesId2vecIndexMap.insert(pair<int,int>(tmpSpeciesId, tmpVecIndex));
		if((genusId2vecIndexMap.find(tmpGenusId) == genusId2vecIndexMap.end())&&(tmpGenusId > 0))
			genusId2vecIndexMap.insert(pair<int,int>(tmpGenusId, tmpVecIndex));
		if((familyId2vecIndexMap.find(tmpFamilyId) == familyId2vecIndexMap.end())&&(tmpFamilyId > 0))
			familyId2vecIndexMap.insert(pair<int,int>(tmpFamilyId, tmpVecIndex));		
		if((orderId2vecIndexMap.find(tmpOrderId) == orderId2vecIndexMap.end())&&(tmpOrderId > 0))
			orderId2vecIndexMap.insert(pair<int,int>(tmpOrderId, tmpVecIndex));
		if((classId2vecIndexMap.find(tmpClassId) == classId2vecIndexMap.end())&&(tmpClassId > 0))
			classId2vecIndexMap.insert(pair<int,int>(tmpClassId, tmpVecIndex));
		if((phylumId2vecIndexMap.find(tmpPhylumId) == phylumId2vecIndexMap.end())&&(tmpPhylumId > 0))
			phylumId2vecIndexMap.insert(pair<int,int>(tmpPhylumId, tmpVecIndex));
		tmpVecIndex ++;
	}
	speciesId2taxoId_ifs.close();

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

		map<int,int>::iterator speciesMapIter = speciesId2vecIndexMap.find(tmpTaxoId);
		if(speciesMapIter != speciesId2vecIndexMap.end()) // found in speciesMap
		{
			//cout << "found in species map" << endl;
			int tmpVecIndex = speciesMapIter->second;
			reformed_ofs << speciesIdVec[tmpVecIndex] << "\t" << genusIdVec[tmpVecIndex] << "\t"
				<< familyIdVec[tmpVecIndex] << "\t" << orderIdVec[tmpVecIndex] << "\t" 
				<< classIdVec[tmpVecIndex] << "\t" << phylumIdVec[tmpVecIndex] << endl;
			continue;
		}
		else
		{
			//cout << "not found in species map" << endl;
			map<int,int>::iterator genusMapIter = genusId2vecIndexMap.find(tmpTaxoId);
			if(genusMapIter != genusId2vecIndexMap.end()) // found in genusMap
			{
				//cout << "found in genus map" << endl;
				int tmpVecIndex = genusMapIter->second;
				reformed_ofs << "-1\t" << genusIdVec[tmpVecIndex] << "\t"
					<< familyIdVec[tmpVecIndex] << "\t" << orderIdVec[tmpVecIndex] << "\t" 
					<< classIdVec[tmpVecIndex] << "\t" << phylumIdVec[tmpVecIndex] << endl;
				continue;				
			}
			else
			{
				//cout << "not found in genus map" << endl;
				map<int,int>::iterator familyMapIter = familyId2vecIndexMap.find(tmpTaxoId);
				if(familyMapIter != familyId2vecIndexMap.end()) // found in genusMap
				{
					//cout << "found in family map" << endl;
					int tmpVecIndex = familyMapIter->second;
					reformed_ofs << "-1\t-1\t"
						<< familyIdVec[tmpVecIndex] << "\t" << orderIdVec[tmpVecIndex] << "\t" 
						<< classIdVec[tmpVecIndex] << "\t" << phylumIdVec[tmpVecIndex] << endl;
					continue;				
				}
				else
				{
					//cout << "not found in family map" << endl;
					map<int,int>::iterator orderMapIter = orderId2vecIndexMap.find(tmpTaxoId);
					if(orderMapIter != orderId2vecIndexMap.end()) // found in genusMap
					{
						//cout << "found in order map" << endl;
						int tmpVecIndex = orderMapIter->second;
						reformed_ofs << "-1\t-1\t-1\t" << orderIdVec[tmpVecIndex] << "\t" 
							<< classIdVec[tmpVecIndex] << "\t" << phylumIdVec[tmpVecIndex] << endl;
						continue;				
					}
					else
					{
						//cout << "not found in order map" << endl;
						map<int,int>::iterator classMapIter = classId2vecIndexMap.find(tmpTaxoId);
						if(classMapIter != classId2vecIndexMap.end()) // found in genusMap
						{
							//cout << "found in class map" << endl;
							int tmpVecIndex = classMapIter->second;
							reformed_ofs << "-1\t-1\t-1\t-1\t" 
								<< classIdVec[tmpVecIndex] << "\t" << phylumIdVec[tmpVecIndex] << endl;
							continue;				
						}
						else
						{
							//cout << "not found in class map" << endl;
							map<int,int>::iterator phylumMapIter = phylumId2vecIndexMap.find(tmpTaxoId);
							if(phylumMapIter != phylumId2vecIndexMap.end()) // found in genusMap
							{
								//cout << "found in phylum map" << endl;
								int tmpVecIndex = phylumMapIter->second;
								reformed_ofs << "-1\t-1\t-1\t-1\t-1\t" << phylumIdVec[tmpVecIndex] << endl;
								continue;				
							}
							else
							{
								//cout << "not found in phylum map" << endl;
								reformed_ofs << "-1\t-1\t-1\t-1\t-1\t-1" << endl;
								continue;	
							}
						}
					}
				}
			}
		}

	}
	raw_ifs.close();
	reformed_ofs.close();
	return 0;
}