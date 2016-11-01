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

int return_read_taxo_id(string& tmpReadId, int taxo_rank)
{
	string tmpTrimmedReadId = tmpReadId.substr(1);
	vector<string> tmpFieldVec_line;	
	parseStr2fieldVec_line(tmpFieldVec_line, tmpTrimmedReadId);
	string tmpTaxoIdStr = tmpFieldVec_line[8-taxo_rank];
	int tmpId = atoi(tmpTaxoIdStr.c_str());
	if(tmpId < 0)
		return -1;
	else
		return tmpId;
}

int return_lothello_assigned_taxo_id(string& tmpLothelloStr, int taxo_rank)
{
	vector<string> tmpFieldVec_tab;
	parseStr2fieldVec_tab(tmpFieldVec_tab, tmpLothelloStr);
	string tmpTaxoIdStr = tmpFieldVec_tab[1];
	int tmpId = atoi(tmpTaxoIdStr.c_str());
	if(tmpId < 0)
		return -1;
	else
		return tmpId;
}

int return_otherTool_assigned_taxo_id(string& tmpOtherToolStr, int taxo_rank)
{
	vector<string> tmpFieldVec_comma;
	parseStr2fieldVec_comma(tmpFieldVec_comma, tmpOtherToolStr);
	string tmpTaxoIdStr = tmpFieldVec_comma[2];
	if(tmpTaxoIdStr == "NA")
		return -1;
	else
	{
		int tmpId = atoi(tmpTaxoIdStr.c_str());
		if(tmpId < 0)
			return -1;
		else
			return tmpId;
	}
}

int main(int argc, char** argv)
{
	if((argc != 9)&&(argc != 11))
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputLothelloResults" << endl;
		cout << "#2 inputOtherToolResults" << endl;
		cout << "#3 incorrectButOtherToolCorrectAssignmentFile" << endl;
		cout << "#4 Fa_or_Fq" << endl;
		cout << "#5 SE_or_PE" << endl;
		cout << "#6 Taxo_name" << endl;
		cout << "#7 inputReadFile_SE or inputReadFile_PE_1" << endl;
		cout << "#8 outputReadFile_SE or inputReadFile_PE_2" << endl;
		cout << "(#9 outputReadFile_PE_1)" << endl;
		cout << "(#10 outputReadFile_PE_2)" << endl;
		exit(1);
	}
	string inputLothelloResults = argv[1];
	string inputOtherToolResults = argv[2];
	string incorrectButOtherToolCorrectAssignmentFile = argv[3];
	string Fa_or_Fq_str = argv[4];
	string SE_or_PE_str = argv[5];
	string Taxo_name = argv[6];

	int taxo_rank;
	if(Taxo_name == "Species")
		taxo_rank = 8;
	else if(Taxo_name == "Genus")
		taxo_rank = 7;
	else if(Taxo_name == "Phylum")
		taxo_rank = 3;
	else
	{
		cout << "tmp Taxo_name is not supported!" << endl << "Taxo_name: " << Taxo_name << endl;
		exit(1);
	}

	string inputReadFile_SE, inputReadFile_PE_1, inputReadFile_PE_2,
		outputReadFile_SE, outputReadFile_PE_1, outputReadFile_PE_2;

	bool Fa_or_Fq_bool;
	if(Fa_or_Fq_str == "Fa")
		Fa_or_Fq_bool = true;
	else if(Fa_or_Fq_str == "Fq")
		Fa_or_Fq_bool = false;
	else
	{
		cout << "invalid Fa_or_Fq_str: " << Fa_or_Fq_str << endl;
		exit(1);
	}

	bool SE_or_PE_bool;
	if((SE_or_PE_str == "SE")&&(argc == 9))
		SE_or_PE_bool = true;
	else if((SE_or_PE_str == "PE")&&(argc == 11))
		SE_or_PE_bool = false;
	else
	{
		cout << "invalid SE_or_PE_str or argc: " << endl << "SE_or_PE_str: " 
			<< SE_or_PE_str << endl << "argc: " << argc << endl;; 
		exit(1);
	}

	if(SE_or_PE_bool)
	{
		inputReadFile_SE = argv[7];
		outputReadFile_SE = argv[8];
	}
	else
	{
		inputReadFile_PE_1 = argv[7];
		inputReadFile_PE_2 = argv[8];
		outputReadFile_PE_1 = argv[9];
		outputReadFile_PE_2 = argv[10];
	}

	ofstream incorrect_ofs(incorrectButOtherToolCorrectAssignmentFile.c_str());
	ifstream lothello_ifs(inputLothelloResults.c_str());
	ifstream otherTool_ifs(inputOtherToolResults.c_str());
	string tmp1stLineInOtherToolResults;
	getline(otherTool_ifs, tmp1stLineInOtherToolResults);
	if(SE_or_PE_bool) // SE
	{
		ifstream SE_ifs(inputReadFile_SE.c_str());
		ofstream SE_ofs(outputReadFile_SE.c_str());
		while((!lothello_ifs.eof())&&(!otherTool_ifs.eof())&&(!SE_ifs.eof()))
		{
			string tmpReadId_SE, tmpLothelloStr, tmpOtherToolStr;
			getline(SE_ifs, tmpReadId_SE);
			getline(lothello_ifs, tmpLothelloStr);
			getline(otherTool_ifs, tmpOtherToolStr);
			if((tmpReadId_SE == "")||(tmpLothelloStr == "")||(tmpOtherToolStr == ""))
				break;
			string tmpReadSeq_SE;
			getline(SE_ifs, tmpReadSeq_SE);
			if(!Fa_or_Fq_bool)
			{
				string tmpReadComment_SE;
				getline(SE_ifs, tmpReadComment_SE);
				string tmpReadQual_SE;
				getline(SE_ifs, tmpReadQual_SE);
			}
			int tmpTruthId = return_read_taxo_id(tmpReadId_SE, taxo_rank);
			int tmpLothelloId = return_lothello_assigned_taxo_id(tmpLothelloStr, taxo_rank);
			int tmpOtherToolId = return_otherTool_assigned_taxo_id(tmpOtherToolStr, taxo_rank);
			if((tmpTruthId >= 0)&&(tmpTruthId == tmpOtherToolId)&&(tmpLothelloId != tmpTruthId))
			{
				if(Fa_or_Fq_bool)
					SE_ofs << tmpReadId_SE << endl;
				else
					SE_ofs << ">" << tmpReadId_SE.substr(1) << endl;
				SE_ofs << tmpReadSeq_SE << endl;
				incorrect_ofs << tmpReadId_SE << "\t" << tmpTruthId << "\t" << tmpLothelloId << "\t" << tmpOtherToolId << endl;
			}
		}
		SE_ifs.close();
		SE_ofs.close();
		incorrect_ofs.close();
	}
	else // PE
	{
		ifstream PE_ifs_1(inputReadFile_PE_1.c_str());
		ifstream PE_ifs_2(inputReadFile_PE_2.c_str());
		ofstream PE_ofs_1(outputReadFile_PE_1.c_str());
		ofstream PE_ofs_2(outputReadFile_PE_2.c_str());
		while((!lothello_ifs.eof())&&(!otherTool_ifs.eof())&&(!PE_ifs_1.eof())&&(!PE_ifs_2.eof()))
		{
			string tmpReadId_PE_1, tmpReadId_PE_2, tmpLothelloStr, tmpOtherToolStr;
			getline(PE_ifs_1, tmpReadId_PE_1);
			getline(PE_ifs_2, tmpReadId_PE_2);
			getline(lothello_ifs, tmpLothelloStr);
			getline(otherTool_ifs, tmpOtherToolStr);
			if((tmpReadId_PE_1 == "")||(tmpReadId_PE_2 == "")||(tmpLothelloStr == "")||(tmpOtherToolStr == ""))
				break;
			string tmpReadSeq_PE_1, tmpReadSeq_PE_2;
			getline(PE_ifs_1, tmpReadSeq_PE_1);
			getline(PE_ifs_2, tmpReadSeq_PE_2);
			if(!Fa_or_Fq_bool)
			{
				string tmpReadComment_PE_1, tmpReadComment_PE_2, tmpReadQual_PE_1, tmpReadQual_PE_2;
				getline(PE_ifs_1, tmpReadComment_PE_1);
				getline(PE_ifs_1, tmpReadQual_PE_1);
				getline(PE_ifs_2, tmpReadComment_PE_2);
				getline(PE_ifs_2, tmpReadQual_PE_2);
			}
			int tmpTruthId = return_read_taxo_id(tmpReadId_PE_1, taxo_rank);
			int tmpLothelloId = return_lothello_assigned_taxo_id(tmpLothelloStr, taxo_rank);
			int tmpOtherToolId = return_otherTool_assigned_taxo_id(tmpOtherToolStr, taxo_rank);
			if((tmpTruthId >= 0)&&(tmpTruthId == tmpOtherToolId)&&(tmpLothelloId != tmpTruthId))
			{
				if(Fa_or_Fq_bool)
				{
					PE_ofs_1 << tmpReadId_PE_1 << endl;
					PE_ofs_2 << tmpReadId_PE_2 << endl;
				}
				else
				{
					PE_ofs_1 << ">" << tmpReadId_PE_1.substr(1) << endl;
					PE_ofs_2 << ">" << tmpReadId_PE_2.substr(1) << endl;
				}
				PE_ofs_1 << tmpReadSeq_PE_1 << endl;
				PE_ofs_2 << tmpReadSeq_PE_2 << endl;
				incorrect_ofs << tmpReadId_PE_1 << "\t" << tmpTruthId << "\t" << tmpLothelloId << "\t" << tmpOtherToolId << endl;
			}								
		}
		PE_ifs_1.close();
		PE_ifs_2.close();
		PE_ofs_1.close();
		PE_ofs_2.close();
		incorrect_ofs.close();
	}
	return 0;
}