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

void parseStr2fieldVec(vector<string>& tmpFieldVec, string& tmpStr)
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
		cout << "#1 inputSpecieCountFile" << endl;
		cout << "#2 inputSpeciesId2taxoIdFile" << endl;
		cout << "#3 outputPrefix" << endl;
		exit(1);
	}
	string inputSpecieCountFile = argv[1];
	string inputSpeciesId2taxoIdFile = argv[2];
	string outputPrefix = argv[3];

	vector<string> speciesIdVec;
	vector<string> genusIdVec;
	vector<string> familyIdVec;
	vector<string> orderIdVec;
	vector<string> classIdVec;
	vector<string> phylumIdVec;

	vector< pair<string, int> > idCountPairVec_genus;
	vector< pair<string, int> > idCountPairVec_family;
	vector< pair<string, int> > idCountPairVec_order;
	vector< pair<string, int> > idCountPairVec_class;
	vector< pair<string, int> > idCountPairVec_phylum;
	ifstream speciesId2taxoId_ifs(inputSpeciesId2taxoIdFile.c_str());
	while(!speciesId2taxoId_ifs.eof())
	{
		string tmpStr;
		getline(speciesId2taxoId_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpFieldVec;
		parseStr2fieldVec(tmpFieldVec, tmpStr);
		string tmpSpeciesStr = tmpFieldVec[0];
		string tmpGenusStr = tmpFieldVec[1];
		string tmpFamilyStr = tmpFieldVec[2];
		string tmpOrderStr = tmpFieldVec[3];
		string tmpClassStr = tmpFieldVec[4];
		string tmpPhylumStr = tmpFieldVec[5];

		speciesIdVec.push_back(tmpSpeciesStr);
		genusIdVec.push_back(tmpGenusStr);
		familyIdVec.push_back(tmpFamilyStr);
		orderIdVec.push_back(tmpOrderStr);
		classIdVec.push_back(tmpClassStr);
		phylumIdVec.push_back(tmpPhylumStr);	

		bool tmpAlreadyExists_bool_genus = false;
		for(int tmp = 0; tmp < idCountPairVec_genus.size(); tmp++)
		{
			if(tmpGenusStr == idCountPairVec_genus[tmp].first)
			{
				tmpAlreadyExists_bool_genus = true;
				break;
			}
		}
		if(!tmpAlreadyExists_bool_genus)
			idCountPairVec_genus.push_back(pair<string,int>(tmpGenusStr, 0));

		bool tmpAlreadyExists_bool_family = false;
		for(int tmp = 0; tmp < idCountPairVec_family.size(); tmp++)
		{
			if(tmpFamilyStr == idCountPairVec_family[tmp].first)
			{
				tmpAlreadyExists_bool_family = true;
				break;
			}
		}
		if(!tmpAlreadyExists_bool_family)
			idCountPairVec_family.push_back(pair<string,int>(tmpFamilyStr, 0));

		bool tmpAlreadyExists_bool_order = false;
		for(int tmp = 0; tmp < idCountPairVec_order.size(); tmp++)
		{
			if(tmpOrderStr == idCountPairVec_order[tmp].first)
			{
				tmpAlreadyExists_bool_order = true;
				break;
			}
		}
		if(!tmpAlreadyExists_bool_order)
			idCountPairVec_order.push_back(pair<string,int>(tmpOrderStr, 0));

		bool tmpAlreadyExists_bool_class = false;
		for(int tmp = 0; tmp < idCountPairVec_class.size(); tmp++)
		{
			if(tmpClassStr == idCountPairVec_class[tmp].first)
			{
				tmpAlreadyExists_bool_class = true;
				break;
			}
		}
		if(!tmpAlreadyExists_bool_class)
			idCountPairVec_class.push_back(pair<string,int>(tmpClassStr, 0));

		bool tmpAlreadyExists_bool_phylum = false;
		for(int tmp = 0; tmp < idCountPairVec_phylum.size(); tmp++)
		{
			if(tmpPhylumStr == idCountPairVec_phylum[tmp].first)
			{
				tmpAlreadyExists_bool_phylum = true;
				break;
			}
		}
		if(!tmpAlreadyExists_bool_phylum)
			idCountPairVec_phylum.push_back(pair<string,int>(tmpPhylumStr, 0));					
	}
	speciesId2taxoId_ifs.close();
	cout << "inputSpecieCountFile: " << inputSpecieCountFile << endl;
	ifstream speciesCount_ifs(inputSpecieCountFile.c_str());
	while(!speciesCount_ifs.eof())
	{
		string tmpStr;
		getline(speciesCount_ifs, tmpStr);
		//cout << "tmpStr: " << tmpStr << endl;
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpSpeciesId = tmpStr.substr(0, tabLoc);
		//cout << "tmpSpeciesId: " << tmpSpeciesId << endl;
		string tmpSpeciesCountStr = tmpStr.substr(tabLoc + 1);
		int tmpSpeciesCount = atoi(tmpSpeciesCountStr.c_str());
		//cout << "tmpSpeciesCount: " << tmpSpeciesCount << endl;
		for(int tmp = 0; tmp < speciesIdVec.size(); tmp++)
		{
			if(tmpSpeciesId == speciesIdVec[tmp])
			{
				string tmpGenusId = genusIdVec[tmp];
				//cout << "tmpGenusId: " << tmpGenusId << endl;
				for(int tmpTaxoIndex = 0; tmpTaxoIndex < idCountPairVec_genus.size(); tmpTaxoIndex ++)
				{
					if(tmpGenusId == idCountPairVec_genus[tmpTaxoIndex].first)
						(idCountPairVec_genus[tmpTaxoIndex].second) += tmpSpeciesCount;
				}
				string tmpFamilyId = familyIdVec[tmp];
				for(int tmpTaxoIndex = 0; tmpTaxoIndex < idCountPairVec_family.size(); tmpTaxoIndex ++)
				{
					if(tmpFamilyId == idCountPairVec_family[tmpTaxoIndex].first)
						(idCountPairVec_family[tmpTaxoIndex].second) += tmpSpeciesCount;
				}
				string tmpOrderId = orderIdVec[tmp];
				for(int tmpTaxoIndex = 0; tmpTaxoIndex < idCountPairVec_order.size(); tmpTaxoIndex ++)
				{
					if(tmpOrderId == idCountPairVec_order[tmpTaxoIndex].first)
						(idCountPairVec_order[tmpTaxoIndex].second) += tmpSpeciesCount;
				}
				string tmpClassId = classIdVec[tmp];
				for(int tmpTaxoIndex = 0; tmpTaxoIndex < idCountPairVec_class.size(); tmpTaxoIndex ++)
				{
					if(tmpClassId == idCountPairVec_class[tmpTaxoIndex].first)
						(idCountPairVec_class[tmpTaxoIndex].second) += tmpSpeciesCount;
				}
				string tmpPhylumId = phylumIdVec[tmp];
				for(int tmpTaxoIndex = 0; tmpTaxoIndex < idCountPairVec_phylum.size(); tmpTaxoIndex ++)
				{
					if(tmpPhylumId == idCountPairVec_phylum[tmpTaxoIndex].first)
						(idCountPairVec_phylum[tmpTaxoIndex].second) += tmpSpeciesCount;
				}
			}
		}
	}
	speciesCount_ifs.close();

	string outputFile_genus = outputPrefix + ".genus";
	string outputFile_family = outputPrefix + ".family";
	string outputFile_order = outputPrefix + ".order";
	string outputFile_class = outputPrefix + ".class";
	string outputFile_phylum = outputPrefix + ".phylum";
	ofstream genus_ofs(outputFile_genus.c_str());
	ofstream family_ofs(outputFile_family.c_str());
	ofstream order_ofs(outputFile_order.c_str());
	ofstream class_ofs(outputFile_class.c_str());
	ofstream phylum_ofs(outputFile_phylum.c_str());
	for(int tmp = 0; tmp < idCountPairVec_genus.size(); tmp++)
		genus_ofs << idCountPairVec_genus[tmp].first << "\t" << idCountPairVec_genus[tmp].second << endl;
	for(int tmp = 0; tmp < idCountPairVec_family.size(); tmp++)
		family_ofs << idCountPairVec_family[tmp].first << "\t" << idCountPairVec_family[tmp].second << endl;
	for(int tmp = 0; tmp < idCountPairVec_order.size(); tmp++)
		order_ofs << idCountPairVec_order[tmp].first << "\t" << idCountPairVec_order[tmp].second << endl;
	for(int tmp = 0; tmp < idCountPairVec_class.size(); tmp++)
		class_ofs << idCountPairVec_class[tmp].first << "\t" << idCountPairVec_class[tmp].second << endl;
	for(int tmp = 0; tmp < idCountPairVec_phylum.size(); tmp++)
		phylum_ofs << idCountPairVec_phylum[tmp].first << "\t" << idCountPairVec_phylum[tmp].second << endl;							
	genus_ofs.close();
	family_ofs.close();
	order_ofs.close();
	class_ofs.close();
	phylum_ofs.close();
	return 0;
}