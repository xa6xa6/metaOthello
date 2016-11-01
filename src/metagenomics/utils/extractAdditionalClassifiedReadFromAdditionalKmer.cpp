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

#define CORRECT_ASSIGNMENT 64320

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputResults_ori" << endl;
		cout << "#2 inputResults_new" << endl;
		cout << "#3 fq_1" << endl;
		cout << "#4 fq_2" << endl;
		cout << "#5 results_dir" << endl;
		exit(1);
	}
	string inputResults_ori = argv[1];
	string inputResults_new = argv[2];
	string fq_1_file = argv[3];
	string fq_2_file = argv[4];
	string results_dir = argv[5]; results_dir += "/";
	string cmd_mkdir_results = "mkdir " + results_dir;
	system(cmd_mkdir_results.c_str());
	string additional_correct_results_file = results_dir + "additional_correct_read_classification.txt";
	string additional_correct_fq_file_1 = results_dir + "additional_correct_read.1.fq";
	string additional_correct_fq_file_2 = results_dir + "additional_correct_read.2.fq";
	string additional_correct_fq_file_merged = results_dir + "additional_correct_read.merged.fq";

	ifstream fq_1_ifs(fq_1_file.c_str());
	ifstream fq_2_ifs(fq_2_file.c_str());
	ifstream results_ori_ifs(inputResults_ori.c_str());
	ifstream results_new_ifs(inputResults_new.c_str());
	ofstream additional_correct_ofs(additional_correct_results_file.c_str());
	ofstream additional_correct_fq_1_ofs(additional_correct_fq_file_1.c_str());
	ofstream additional_correct_fq_2_ofs(additional_correct_fq_file_2.c_str());
	ofstream additional_correct_fq_merged_ofs(additional_correct_fq_file_merged.c_str());

	while((!fq_1_ifs.eof())&&(!fq_2_ifs.eof())&&(!results_ori_ifs.eof())&&(!results_new_ifs.eof()))
	{
		string tmpFqStr_1_id, tmpFqStr_2_id, tmpFqStr_1_seq, tmpFqStr_2_seq, 
			tmpFqStr_1_comment, tmpFqStr_2_comment, tmpFqStr_1_qual, tmpFqStr_2_qual,
			tmpAssignStr_ori, tmpAssignStr_new;
		getline(fq_1_ifs, tmpFqStr_1_id);
		getline(fq_2_ifs, tmpFqStr_2_id);
		getline(results_ori_ifs, tmpAssignStr_ori);
		getline(results_new_ifs, tmpAssignStr_new);
		if((tmpFqStr_1_id == "")||(tmpFqStr_2_id == "")||(tmpAssignStr_ori == "")||(tmpAssignStr_new == ""))
			break;
		getline(fq_1_ifs, tmpFqStr_1_seq);
		getline(fq_1_ifs, tmpFqStr_1_comment);
		getline(fq_1_ifs, tmpFqStr_1_qual);
		getline(fq_2_ifs, tmpFqStr_2_seq);
		getline(fq_2_ifs, tmpFqStr_2_comment);
		getline(fq_2_ifs, tmpFqStr_2_qual);

		int assignStr_ori_tabLoc_1 = tmpAssignStr_ori.find("\t");
		int assignStr_ori_tabLoc_2 = tmpAssignStr_ori.find("\t", assignStr_ori_tabLoc_1 + 1);
		int assignStr_new_tabLoc_1 = tmpAssignStr_new.find("\t");
		int assignStr_new_tabLoc_2 = tmpAssignStr_new.find("\t", assignStr_new_tabLoc_1 + 1);
		string tmpAssignIdStr_ori = tmpAssignStr_ori.substr(assignStr_ori_tabLoc_1 + 1, assignStr_ori_tabLoc_2 - assignStr_ori_tabLoc_1 - 1);
		string tmpAssignIdStr_new = tmpAssignStr_new.substr(assignStr_new_tabLoc_1 + 1, assignStr_new_tabLoc_2 - assignStr_new_tabLoc_1 - 1);
		int tmpAssignId_ori = atoi(tmpAssignIdStr_ori.c_str());
		int tmpAssignId_new = atoi(tmpAssignIdStr_new.c_str());
		if((tmpAssignId_new = CORRECT_ASSIGNMENT)&&(tmpAssignId_ori != CORRECT_ASSIGNMENT))
		{
			additional_correct_ofs << tmpAssignStr_ori << "\t" << tmpAssignStr_new << endl;
			additional_correct_fq_1_ofs << tmpFqStr_1_id << endl << tmpFqStr_1_seq << endl 
				<< tmpFqStr_1_comment << endl << tmpFqStr_1_qual << endl;
			additional_correct_fq_2_ofs << tmpFqStr_2_id << endl << tmpFqStr_2_seq << endl 
				<< tmpFqStr_2_comment << endl << tmpFqStr_2_qual << endl;				
			additional_correct_fq_merged_ofs << tmpFqStr_1_id << endl << tmpFqStr_1_seq << tmpFqStr_2_seq << endl;
		}
	}


	fq_1_ifs.close();
	fq_2_ifs.close();
	results_ori_ifs.close();
	results_new_ifs.close();
	additional_correct_ofs.close();
	additional_correct_fq_1_ofs.close();
	additional_correct_fq_2_ofs.close();
	additional_correct_fq_merged_ofs.close();
	return 0;
}