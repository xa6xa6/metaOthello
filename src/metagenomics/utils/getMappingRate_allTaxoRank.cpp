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

void parseStr2fieldVec_line(vector<string>& tmpFieldVec, string& tmpStr)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find("_", startLoc);
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
	if(argc != 3)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputReadAssignmentFile" << endl;
		cout << "#2 outputStatsFile" << endl;
		//cout << "#3 taxo_name" << endl;
		exit(1);
	}
	string inputReadAssignmentFile = argv[1];
	string outputStatsFile = argv[2];

	int total_num = 0;
	int mapped_num_total = 0;
	int unmapped_num_total = 0;
	int mapped_num_species = 0;
	int mapped_num_genus = 0;
	int mapped_num_family = 0;
	int mapped_num_order = 0;
	int mapped_num_class = 0;
	int mapped_num_phylum = 0;
	int unmapped_num_species = 0;
	int unmapped_num_genus = 0;
	int unmapped_num_family = 0;
	int unmapped_num_order = 0;
	int unmapped_num_class = 0;
	int unmapped_num_phylum = 0;	

	ifstream readAssignment_ifs(inputReadAssignmentFile.c_str());
	ofstream stats_ofs(outputStatsFile.c_str());

	unsigned long long tmpLineNO = 0;
	while(!readAssignment_ifs.eof())
	{
		string tmpStr;
		getline(readAssignment_ifs, tmpStr);
		if(tmpStr == "")
			break;
        tmpLineNO ++;
        unsigned long long tmpThousandIndex = tmpLineNO / 5000000;
        if(tmpLineNO == tmpThousandIndex * 5000000)
            cout << "Processed Line #: " << tmpLineNO << endl;

        total_num ++;
        vector<string> tmpFieldVec_tab;
        parseStr2fieldVec_tab(tmpFieldVec_tab, tmpStr);
        
        string tmpAssignment_species_str = tmpFieldVec_tab[1];
        string tmpAssignment_genus_str = tmpFieldVec_tab[2];
        string tmpAssignment_family_str = tmpFieldVec_tab[3];
        string tmpAssignment_order_str = tmpFieldVec_tab[4];
        string tmpAssignment_class_str = tmpFieldVec_tab[5];
        string tmpAssignment_phylum_str = tmpFieldVec_tab[6];

        int tmpAssignment_species_int = atoi(tmpAssignment_species_str.c_str());
        int tmpAssignment_genus_int = atoi(tmpAssignment_genus_str.c_str());
        int tmpAssignment_family_int = atoi(tmpAssignment_family_str.c_str());
        int tmpAssignment_order_int = atoi(tmpAssignment_order_str.c_str());
        int tmpAssignment_class_int = atoi(tmpAssignment_class_str.c_str());
        int tmpAssignment_phylum_int = atoi(tmpAssignment_phylum_str.c_str());
		
		if(tmpAssignment_species_int >= 0)
			mapped_num_species ++;
		else
			unmapped_num_species ++;

		if(tmpAssignment_genus_int >= 0)
			mapped_num_genus ++;
		else
			unmapped_num_genus ++;

		if(tmpAssignment_family_int >= 0)
			mapped_num_family ++;
		else
			unmapped_num_family ++;

		if(tmpAssignment_order_int >= 0)
			mapped_num_order ++;
		else
			unmapped_num_order ++;	

		if(tmpAssignment_class_int >= 0)
			mapped_num_class ++;
		else
			unmapped_num_class ++;	

		if(tmpAssignment_phylum_int >= 0)
			mapped_num_phylum ++;
		else
			unmapped_num_phylum ++;

		if((tmpAssignment_species_int >= 0)||(tmpAssignment_genus_int >= 0)
			||(tmpAssignment_family_int >= 0)||(tmpAssignment_order_int >= 0)
			||(tmpAssignment_class_int >= 0)||(tmpAssignment_phylum_int >= 0))			
			mapped_num_total ++;
		else
			unmapped_num_total ++;
	}
	double mappingRate_total = ((double)mapped_num_total/(double)total_num) * 100;
	double mappingRate_phylum = ((double)mapped_num_phylum/(double)total_num) * 100;
	double mappingRate_class = ((double)mapped_num_class/(double)total_num) * 100;
	double mappingRate_order = ((double)mapped_num_order/(double)total_num) * 100;
	double mappingRate_family = ((double)mapped_num_family/(double)total_num) * 100;
	double mappingRate_genus = ((double)mapped_num_genus/(double)total_num) * 100;
	double mappingRate_species = ((double)mapped_num_species/(double)total_num) * 100;

	stats_ofs << "Total     #:\t" << total_num << endl << endl;
	stats_ofs << "Taxo_rank\tMapping_Rate\tMapped_read_#\tUnmapped_read_#" << endl;
	stats_ofs << "Total\t" << mappingRate_total << "%\t" << mapped_num_total << "\t" << unmapped_num_total << endl 
			<< "Phylum\t" << mappingRate_phylum << "%\t" << mapped_num_phylum << "\t" << unmapped_num_phylum << endl
			<< "Class\t" << mappingRate_class << "%\t" << mapped_num_class << "\t" << unmapped_num_class << endl
			<< "Order\t" << mappingRate_order << "%\t" << mapped_num_order << "\t" << unmapped_num_order << endl
			<< "Family\t" << mappingRate_family << "%\t" << mapped_num_family << "\t" << unmapped_num_family << endl
			<< "Genus\t" << mappingRate_genus << "%\t" << mapped_num_genus << "\t" << unmapped_num_genus << endl
			<< "Species\t" << mappingRate_species << "%\t" << mapped_num_species << "\t" << unmapped_num_species << endl;

	stats_ofs.close();
	readAssignment_ifs.close();
	return 0;
}