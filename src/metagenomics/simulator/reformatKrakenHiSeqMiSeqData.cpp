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

string get_formattedSpeciesName_from_readId(string& tmpReadIdStr)
{
	int lineLoc_1 = tmpReadIdStr.find("_");
	int lineLoc_2 = tmpReadIdStr.find("_", lineLoc_1 + 1);
	return tmpReadIdStr.substr(1, lineLoc_2 - 1);
}

string get_taxoInfoStr_from_formattedSpeciesName(string& tmp_formattedSpeciesName, 
	vector<string>& species_formattedNameVec, vector<string>& species_taxoInfoVec)
{
	string toReturn_taxoInfoStr = "NULL";
	for(int tmp = 0; tmp < species_formattedNameVec.size(); tmp++)
	{
		if(species_formattedNameVec[tmp] == tmp_formattedSpeciesName)
			return species_taxoInfoVec[tmp];
	}
	return toReturn_taxoInfoStr;
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

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 taxoInfo_overall_file" << endl;
		cout << "#2 inputFa" << endl;
		cout << "#3 outputFile_prefix" << endl;
 		exit(1);
 	}
 	string taxoInfo_overall_file = argv[1];
 	string inputFa = argv[2];
 	string outputFile_prefix = argv[3];
 	string outputFa = outputFile_prefix + ".fa";
 	vector<string> species_formattedNameVec;
 	vector<string> species_taxoInfoVec;
 	ifstream taxoInfo_overall_ifs(taxoInfo_overall_file.c_str());
 	string tmp1stLine;
 	getline(taxoInfo_overall_ifs, tmp1stLine);
 	while(!taxoInfo_overall_ifs.eof())
 	{
 		string tmpStr;
 		getline(taxoInfo_overall_ifs, tmpStr);
 		if(tmpStr == "")
 			break;
 		vector<string> tmpTaxoInfoFieldVec;
 		parseStr2fieldVec_tab(tmpTaxoInfoFieldVec, tmpStr);
 		string species_id = tmpTaxoInfoFieldVec[1];
 		string species_name = tmpTaxoInfoFieldVec[2];
 		string genus_id = tmpTaxoInfoFieldVec[4];
 		string genus_name = tmpTaxoInfoFieldVec[5];
  		string family_id = tmpTaxoInfoFieldVec[7];
 		string family_name = tmpTaxoInfoFieldVec[8];
  		string order_id = tmpTaxoInfoFieldVec[10];
 		string order_name = tmpTaxoInfoFieldVec[11];
  		string class_id = tmpTaxoInfoFieldVec[13];
 		string class_name = tmpTaxoInfoFieldVec[14];
  		string phylum_id = tmpTaxoInfoFieldVec[16];
 		string phylum_name = tmpTaxoInfoFieldVec[17];

 		int blankLoc = species_name.find(" ");
 		string tmp_species_name_formatted = species_name.substr(0, 1) + "_" + species_name.substr(blankLoc + 1);
 		string tmpTaxoInfoStr = species_id + "_" + genus_id + "_" + family_id 
 			+ "_" + order_id + "_" + class_id + "_" + phylum_id;
 		species_formattedNameVec.push_back(tmp_species_name_formatted);
 		species_taxoInfoVec.push_back(tmpTaxoInfoStr);
 	}
 	taxoInfo_overall_ifs.close();

 	ifstream rawfa_ifs(inputFa.c_str());
 	string outputFa_twoLineFa = outputFa + ".inter2Line";
 	ofstream inter2Line_ofs(outputFa_twoLineFa.c_str());

 	string tmpReadId = "";
 	string tmpReadSeq = "";
 	getline(rawfa_ifs, tmpReadId);
 	while(!rawfa_ifs.eof())
 	{
 		string tmpStr;
 		getline(rawfa_ifs, tmpStr);
 		if(tmpStr == "")
 			break;
 		if(tmpStr.at(0) == '>')
 		{

 			inter2Line_ofs << tmpReadId << endl << tmpReadSeq << endl;
 			tmpReadId = tmpStr;
 			tmpReadSeq = "";
 		}
 		else//
 			tmpReadSeq += tmpStr;
 	}
 	inter2Line_ofs << tmpReadId << endl << tmpReadSeq << endl;
 	rawfa_ifs.close();
 	inter2Line_ofs.close();

 	string lastReadId_raw_formattedSpeciesName = "";
 	string lastReadId_raw_taxoInfoStr = ""; 	
 	// lastReadId_raw_formattedSpeciesName = get_formattedSpeciesName_from_readId(tmpReadId);
 	// lastReadId_raw_taxoInfoStr = get_taxoInfoStr_from_formattedSpeciesName(
 	// 	lastReadId_raw_formattedSpeciesName, species_formattedNameVec, species_taxoInfoVec);
 	ifstream inter2Line_ifs(outputFa_twoLineFa.c_str());
 	string outputFa_invalid = outputFile_prefix + ".invalid.fa";
 	ofstream fa_invalid_ofs(outputFa_invalid.c_str());
 	ofstream fa_ofs(outputFa.c_str());
 	unsigned long long readNum = 0;
 	while(!inter2Line_ifs.eof())
 	{
 		string tmpReadIdStr, tmpReadSeq;
 		getline(inter2Line_ifs, tmpReadIdStr);
 		if(tmpReadIdStr == "")
 			break;
 		getline(inter2Line_ifs, tmpReadSeq);
 		readNum ++;
 		string tmp_formattedSpeciesName = get_formattedSpeciesName_from_readId(tmpReadIdStr);
 		if(tmp_formattedSpeciesName != lastReadId_raw_formattedSpeciesName)
 		{
 			lastReadId_raw_formattedSpeciesName = tmp_formattedSpeciesName;
 			lastReadId_raw_taxoInfoStr = get_taxoInfoStr_from_formattedSpeciesName(
 				lastReadId_raw_formattedSpeciesName, species_formattedNameVec, species_taxoInfoVec);
 		}
 		if(lastReadId_raw_taxoInfoStr == "NULL")
 			fa_invalid_ofs << tmpReadIdStr << "_" << readNum << endl << tmpReadSeq << endl;
 		else
 			fa_ofs << ">" << lastReadId_raw_taxoInfoStr << "_" << readNum << endl << tmpReadSeq << endl;
 	}
 	fa_ofs.close();
 	fa_invalid_ofs.close();
 	inter2Line_ifs.close();
	return 0;
}