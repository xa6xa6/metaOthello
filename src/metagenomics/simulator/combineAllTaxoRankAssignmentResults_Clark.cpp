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

void parseStr2fieldVec_comma(vector<string>& tmpFieldVec, string& tmpStr)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find(",", startLoc);
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
	if(argc != 8)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputReadAssignmentFile_Species or NULL" << endl;
		cout << "#2 inputReadAssignmentFile_Genus or NULL" << endl;
		cout << "#3 inputReadAssignmentFile_Family or NULL" << endl;
		cout << "#4 inputReadAssignmentFile_Order or NULL" << endl;
		cout << "#5 inputReadAssignmentFile_Class or NULL" << endl;
		cout << "#6 inputReadAssignmentFile_Phylum or NULL" << endl;
		cout << "#7 outputReadAssignmentFile_combined" << endl;
		exit(1);
	}
	string inputReadAssignmentFile_Species = argv[1];
	string inputReadAssignmentFile_Genus = argv[2];
	string inputReadAssignmentFile_Family = argv[3];
	string inputReadAssignmentFile_Order = argv[4];
	string inputReadAssignmentFile_Class = argv[5];
	string inputReadAssignmentFile_Phylum = argv[6];
	string outputReadAssignmentFile_combined = argv[7];

	bool inputReadAssignmentFile_Species_valid_bool = true;
	bool inputReadAssignmentFile_Genus_valid_bool = true;
	bool inputReadAssignmentFile_Family_valid_bool = true;
	bool inputReadAssignmentFile_Order_valid_bool = true;
	bool inputReadAssignmentFile_Class_valid_bool = true;
	bool inputReadAssignmentFile_Phylum_valid_bool = true;
	if(inputReadAssignmentFile_Species == "NULL")
		inputReadAssignmentFile_Species_valid_bool = false;
	if(inputReadAssignmentFile_Genus == "NULL")
		inputReadAssignmentFile_Genus_valid_bool = false;
	if(inputReadAssignmentFile_Family == "NULL")
		inputReadAssignmentFile_Family_valid_bool = false;
	if(inputReadAssignmentFile_Order == "NULL")
		inputReadAssignmentFile_Order_valid_bool = false;
	if(inputReadAssignmentFile_Class == "NULL")
		inputReadAssignmentFile_Class_valid_bool = false;
	if(inputReadAssignmentFile_Phylum == "NULL")
		inputReadAssignmentFile_Phylum_valid_bool = false;									

	ifstream species_ifs; 
	if(inputReadAssignmentFile_Species_valid_bool)
		species_ifs.open(inputReadAssignmentFile_Species.c_str());
	ifstream genus_ifs;
	if(inputReadAssignmentFile_Genus_valid_bool)
		genus_ifs.open(inputReadAssignmentFile_Genus.c_str());
	ifstream family_ifs;
	if(inputReadAssignmentFile_Family_valid_bool)
		family_ifs.open(inputReadAssignmentFile_Family.c_str());
	ifstream order_ifs;
	if(inputReadAssignmentFile_Order_valid_bool)
		order_ifs.open(inputReadAssignmentFile_Order.c_str());
	ifstream class_ifs;
	if(inputReadAssignmentFile_Class_valid_bool)
		class_ifs.open(inputReadAssignmentFile_Class.c_str());
	ifstream phylum_ifs;
	if(inputReadAssignmentFile_Phylum_valid_bool)
		phylum_ifs.open(inputReadAssignmentFile_Phylum.c_str());
	ofstream combined_ofs;
	combined_ofs.open(outputReadAssignmentFile_combined.c_str());

	string tmpStr;
	if(inputReadAssignmentFile_Species_valid_bool)
		getline(species_ifs, tmpStr);
	if(inputReadAssignmentFile_Genus_valid_bool)
		getline(genus_ifs, tmpStr);
	if(inputReadAssignmentFile_Family_valid_bool)
		getline(family_ifs, tmpStr);			
	if(inputReadAssignmentFile_Order_valid_bool)
		getline(order_ifs, tmpStr);
	if(inputReadAssignmentFile_Class_valid_bool)
		getline(class_ifs, tmpStr);
	if(inputReadAssignmentFile_Phylum_valid_bool)
		getline(phylum_ifs, tmpStr);	

	while(1)
	{
		if(((species_ifs.eof())&&inputReadAssignmentFile_Species_valid_bool)
			||((genus_ifs.eof())&&inputReadAssignmentFile_Genus_valid_bool)
			||((family_ifs.eof())&&inputReadAssignmentFile_Family_valid_bool)
			||((order_ifs.eof())&&inputReadAssignmentFile_Order_valid_bool)
			||((class_ifs.eof())&&inputReadAssignmentFile_Class_valid_bool)
			||((phylum_ifs.eof())&&inputReadAssignmentFile_Phylum_valid_bool))
			break;
		
		string tmpStr_species, tmpStr_genus, tmpStr_family, tmpStr_order, tmpStr_class, tmpStr_phylum;
		if(inputReadAssignmentFile_Species_valid_bool)
			getline(species_ifs, tmpStr_species);
		if(inputReadAssignmentFile_Genus_valid_bool)
			getline(genus_ifs, tmpStr_genus);
		if(inputReadAssignmentFile_Family_valid_bool)
			getline(family_ifs, tmpStr_family);
		if(inputReadAssignmentFile_Order_valid_bool)
			getline(order_ifs, tmpStr_order);
		if(inputReadAssignmentFile_Class_valid_bool)
			getline(class_ifs, tmpStr_class);
		if(inputReadAssignmentFile_Phylum_valid_bool)
			getline(phylum_ifs, tmpStr_phylum);
		
		if(((tmpStr_species == "")&&inputReadAssignmentFile_Species_valid_bool)
			||((tmpStr_genus == "")&&inputReadAssignmentFile_Genus_valid_bool)
			||((tmpStr_family == "")&&inputReadAssignmentFile_Family_valid_bool)
			||((tmpStr_order == "")&&inputReadAssignmentFile_Order_valid_bool)
			||((tmpStr_class == "")&&inputReadAssignmentFile_Class_valid_bool)
			||((tmpStr_phylum == "")&&inputReadAssignmentFile_Phylum_valid_bool))
			break;
		
		string tmpReadId_species, tmpReadId_genus, tmpReadId_family, 
			tmpReadId_order, tmpReadId_class, tmpReadId_phylum;
		string tmpReadAssignmentStr_species, tmpReadAssignmentStr_genus, tmpReadAssignmentStr_family,
			tmpReadAssignmentStr_order, tmpReadAssignmentStr_class, tmpReadAssignmentStr_phylum;
		int tmpReadAssignmentId_species, tmpReadAssignmentId_genus, tmpReadAssignmentId_family,
			tmpReadAssignmentId_order, tmpReadAssignmentId_class, tmpReadAssignmentId_phylum;
		
		if(inputReadAssignmentFile_Species_valid_bool)
		{
			vector<string> tmpCommaFieldVec;
			parseStr2fieldVec_comma(tmpCommaFieldVec, tmpStr_species);
			tmpReadId_species = tmpCommaFieldVec[0] + "," + tmpCommaFieldVec[1];
			tmpReadAssignmentStr_species = tmpCommaFieldVec[2];
			if(tmpReadAssignmentStr_species == "NA")
				tmpReadAssignmentId_species = -1;
			else
				tmpReadAssignmentId_species = atoi(tmpReadAssignmentStr_species.c_str());
		}
		else
			tmpReadAssignmentId_species = -1;

		if(inputReadAssignmentFile_Genus_valid_bool)
		{
			vector<string> tmpCommaFieldVec;
			parseStr2fieldVec_comma(tmpCommaFieldVec, tmpStr_genus);
			tmpReadId_genus = tmpCommaFieldVec[0] + "," + tmpCommaFieldVec[1];
			tmpReadAssignmentStr_genus = tmpCommaFieldVec[2];
			if(tmpReadAssignmentStr_genus == "NA")
				tmpReadAssignmentId_genus = -1;
			else
				tmpReadAssignmentId_genus = atoi(tmpReadAssignmentStr_genus.c_str());
		}
		else
			tmpReadAssignmentId_genus = -1;

		if(inputReadAssignmentFile_Family_valid_bool)
		{
			vector<string> tmpCommaFieldVec;
			parseStr2fieldVec_comma(tmpCommaFieldVec, tmpStr_family);
			tmpReadId_family = tmpCommaFieldVec[0] + "," + tmpCommaFieldVec[1];
			tmpReadAssignmentStr_family = tmpCommaFieldVec[2];
			if(tmpReadAssignmentStr_family == "NA")
				tmpReadAssignmentId_family = -1;
			else
				tmpReadAssignmentId_family = atoi(tmpReadAssignmentStr_family.c_str());
		}
		else
			tmpReadAssignmentId_family = -1;

		if(inputReadAssignmentFile_Order_valid_bool)
		{
			vector<string> tmpCommaFieldVec;
			parseStr2fieldVec_comma(tmpCommaFieldVec, tmpStr_order);
			tmpReadId_order = tmpCommaFieldVec[0] + "," + tmpCommaFieldVec[1];
			tmpReadAssignmentStr_order = tmpCommaFieldVec[2];
			if(tmpReadAssignmentStr_order == "NA")
				tmpReadAssignmentId_order = -1;
			else
				tmpReadAssignmentId_order = atoi(tmpReadAssignmentStr_order.c_str());
		}
		else
			tmpReadAssignmentId_order = -1;

		if(inputReadAssignmentFile_Class_valid_bool)
		{
			vector<string> tmpCommaFieldVec;
			parseStr2fieldVec_comma(tmpCommaFieldVec, tmpStr_class);
			tmpReadId_class = tmpCommaFieldVec[0] + "," + tmpCommaFieldVec[1];
			tmpReadAssignmentStr_class = tmpCommaFieldVec[2];
			if(tmpReadAssignmentStr_class == "NA")
				tmpReadAssignmentId_class = -1;
			else
				tmpReadAssignmentId_class = atoi(tmpReadAssignmentStr_class.c_str());
		}
		else
			tmpReadAssignmentId_class = -1;		

		if(inputReadAssignmentFile_Phylum_valid_bool)
		{
			vector<string> tmpCommaFieldVec;
			parseStr2fieldVec_comma(tmpCommaFieldVec, tmpStr_phylum);
			tmpReadId_phylum = tmpCommaFieldVec[0] + "," + tmpCommaFieldVec[1];
			tmpReadAssignmentStr_phylum = tmpCommaFieldVec[2];
			if(tmpReadAssignmentStr_phylum == "NA")
				tmpReadAssignmentId_phylum = -1;
			else
				tmpReadAssignmentId_phylum = atoi(tmpReadAssignmentStr_phylum.c_str());
		}
		else
			tmpReadAssignmentId_phylum = -1;

		// if((tmpReadId_species != tmpReadId_genus)&&
		// 	||(tmpReadId_species != tmpReadId_family)
		// 	||(tmpReadId_species != tmpReadId_order)
		// 	||(tmpReadId_species != tmpReadId_class)
		// 	||(tmpReadId_species != tmpReadId_phylum))
		// {
		// 	cout << "((tmpReadId_species != tmpReadId_genus)||(tmpReadId_species != tmpReadId_family)" << endl;
		// 	cout << "||(tmpReadId_species != tmpReadId_order)||(tmpReadId_species != tmpReadId_class)" << endl;
		// 	cout << "||(tmpReadId_species != tmpReadId_phylum))" << endl;
		// 	cout << "tmpReadId_species: " << tmpReadId_species << endl;
		// 	cout << "tmpReadId_genus: " << tmpReadId_genus << endl;
		// 	cout << "tmpReadId_family: " << tmpReadId_family << endl;
		// 	cout << "tmpReadId_order: " << tmpReadId_order << endl;
		// 	cout << "tmpReadId_class: " << tmpReadId_class << endl;
		// 	cout << "tmpReadId_phylum: " << tmpReadId_phylum << endl;
		// 	exit(1);
		// }
		combined_ofs << tmpReadId_species << "\t" << tmpReadAssignmentId_species << "\t"
			<< tmpReadAssignmentId_genus << "\t" << tmpReadAssignmentId_family << "\t"
			<< tmpReadAssignmentId_order << "\t" << tmpReadAssignmentId_class << "\t"
			<< tmpReadAssignmentId_phylum << endl;

	}

	species_ifs.close();
	genus_ifs.close();
	family_ifs.close();
	order_ifs.close();
	class_ifs.close();
	phylum_ifs.close();
	return 0;
}