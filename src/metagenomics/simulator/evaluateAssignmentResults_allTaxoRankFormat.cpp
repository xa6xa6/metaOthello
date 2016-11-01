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
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputReadAssignmentFile" << endl;
		cout << "#2 outputDir" << endl;
		cout << "#3 taxo_name" << endl;
		exit(1);
	}
	string inputReadAssignmentFile = argv[1];
	string outputDir = argv[2];
	int taxo_rank;
	string taxo_level_name = argv[3];
	if(taxo_level_name == "Phylum")
		taxo_rank = 3;
	else if(taxo_level_name == "Class")
		taxo_rank = 4;
	else if(taxo_level_name == "Order")
		taxo_rank = 5;
	else if(taxo_level_name == "Family")
		taxo_rank = 6;
	else if(taxo_level_name == "Genus")
		taxo_rank = 7;
	else if(taxo_level_name == "Species")
		taxo_rank = 8;
	else
	{
		cout << "invalid taxo_level_name: " << taxo_level_name << endl;
		exit(1);
	}


	outputDir += "/";
	string cmd_mkdir = "mkdir " + outputDir;
	system(cmd_mkdir.c_str());
	string log_file = outputDir + "log.txt";
	ofstream log_ofs(log_file.c_str());
	log_ofs << "start to evaluate read assignment results" << endl;

	ifstream readAssignment_ifs(inputReadAssignmentFile.c_str());
	string unlabeledReadFile = outputDir + "unlabeled.readAssign.txt";
	string correctlyAssignedReadFile = outputDir + "correct.readAssign.txt";
	string incorrectlyAssignedReadFile = outputDir + "incorrect.readAssign.txt";
	string unmappedReadFile = outputDir + "unmapped.readAssign.txt";
	string statsFile = outputDir + "stats.txt";
	
	ofstream unlabeled_ofs(unlabeledReadFile.c_str());
	ofstream correct_ofs(correctlyAssignedReadFile.c_str());
	ofstream incorrect_ofs(incorrectlyAssignedReadFile.c_str());
	ofstream unmapped_ofs(unmappedReadFile.c_str());
	ofstream stats_ofs(statsFile.c_str());

	int total_num = 0;
	int labeled_num = 0;
	int unlabeled_num = 0;
	int labeled_correct_num = 0;
	int labeled_incorrect_num = 0;
	int labeled_unmapped_num = 0;

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
        
        string tmpAssignmentStr = tmpFieldVec_tab[8 - taxo_rank + 1];
        int tmpAssignmentId = atoi(tmpAssignmentStr.c_str());
        
        string tmpReadIdStr = tmpFieldVec_tab[0];
        string tmpReadIdStr_trimmed;
        if((tmpReadIdStr.at(0) == '>')||(tmpReadIdStr.at(0) == '@'))
        	tmpReadIdStr_trimmed = tmpReadIdStr.substr(1);
        else
        	tmpReadIdStr_trimmed = tmpReadIdStr;
        vector<string> tmpFieldVec_line;
        parseStr2fieldVec_line(tmpFieldVec_line, tmpReadIdStr_trimmed);
        string tmpReadTrueTaxoIdStr = tmpFieldVec_line[8 - taxo_rank];
        int tmpReadTrueTaxoId = atoi(tmpReadTrueTaxoIdStr.c_str());

        if(tmpReadTrueTaxoId < 0)
        {	
        	unlabeled_num ++;
        	unlabeled_ofs << tmpStr << endl;
        }
        else
        {
        	labeled_num ++;
        	if(tmpAssignmentId < 0)
        	{	
        		labeled_unmapped_num ++;
        		unmapped_ofs << tmpStr << endl;
        	}
        	else if(tmpAssignmentId == tmpReadTrueTaxoId)
        	{
        		labeled_correct_num ++;
        		correct_ofs << tmpStr << endl;
        	}
        	else
        	{	
        		labeled_incorrect_num ++;
        		incorrect_ofs << tmpStr << endl;
        	}
        }
	}

	stats_ofs << "Total     #:\t" << total_num << endl << endl;
	
	stats_ofs << "unlabeled #:\t" << unlabeled_num << endl << endl;

	stats_ofs << "Labeled   #:\t" << labeled_num << endl;
	stats_ofs << "Correct   #:\t" << labeled_correct_num << endl;
	stats_ofs << "Incorrect #:\t" << labeled_incorrect_num << endl;
	stats_ofs << "Unmapped  #:\t" << labeled_unmapped_num << endl << endl;

	double sensi_perc = ((double)labeled_correct_num/(double)labeled_num) * 100;
	double preci_perc = ((double)labeled_correct_num/(double)(labeled_correct_num + labeled_incorrect_num)) * 100;
	double fscore = 2/(100/sensi_perc + 100/preci_perc);

	stats_ofs << "Sensi_perc:\t" << sensi_perc << "%" << endl;
	stats_ofs << "Preci_perc:\t" << preci_perc << "%" << endl;
	stats_ofs << "F-score   :\t" << fscore << endl;

	stats_ofs.close();
	unlabeled_ofs.close();
	correct_ofs.close();
	incorrect_ofs.close();
	unmapped_ofs.close();
	readAssignment_ifs.close();
	log_ofs.close();
	return 0;
}