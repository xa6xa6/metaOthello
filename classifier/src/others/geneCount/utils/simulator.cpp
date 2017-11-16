// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <bitset>
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

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputArtBin" << endl;
		// /scratch/xli262/gcOthello/gene2faFile.txt
		cout << "#2 gene2faFileReadCountListFile" << endl;
		cout << "#3 outputDir" << endl;
		exit(1);
	}
	string inputArtBin = argv[1];
	string gene2faFileReadCountListFile = argv[2];
	string outputDir = argv[3];
	inputArtBin += "/";
	outputDir += "/";
	string mkdir_output = "mkdir " + outputDir;
	system(mkdir_output.c_str());
	string fqFile_end1 = outputDir + "simu.1.fq";
	string fqFile_end2 = outputDir + "simu.2.fq";
	ofstream simu_1_ofs(fqFile_end1.c_str());
	ofstream simu_2_ofs(fqFile_end2.c_str());
	string log_file = outputDir + "log.txt";
	ofstream log_ofs(log_file.c_str());
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    log_ofs << endl << "[" << asctime(local) << "start to simulate reads ..." << endl;	
	log_ofs << "Command Line:" << endl;
	for(int tmp = 0; tmp < argc; tmp++)
		log_ofs << "#" << tmp << "\t" << argv[tmp] << endl; 
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    log_ofs << endl << "[" << asctime(local) << "start to initiate options for art" << endl;
	string opt_sequencing_platform = "HS25";
	string opt_length = "100 -m 200 -s 10 -p";

	vector<string> geneFaPathVec;
	vector<string> geneIdVec;
	vector<int> geneCountVec;
	ifstream gene2faFile_ifs(gene2faFileReadCountListFile.c_str());
	while(!gene2faFile_ifs.eof())
	{
		string tmpStr;
		getline(gene2faFile_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc + 1);
		string tmpGeneId = tmpStr.substr(0, tabLoc);
		string tmpGeneFaFile = tmpStr.substr(tabLoc + 1, tabLoc_2 - tabLoc - 1);
		string tmpGeneCountStr = tmpStr.substr(tabLoc_2 + 1);
		int tmpGeneCount = atoi(tmpGeneCountStr.c_str());
		geneFaPathVec.push_back(tmpGeneFaFile);
		geneIdVec.push_back(tmpGeneId);
		geneCountVec.push_back(tmpGeneCount);
	}
	gene2faFile_ifs.close();

	string outputDir_rawSimulatedFa = outputDir + "raw_simulated_fa/";
	string mkdir_rawSimulatedFa = "mkdir " + outputDir_rawSimulatedFa;
	system(mkdir_rawSimulatedFa.c_str());
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    log_ofs << endl << "[" << asctime(local) << "start to simulate reads for each fa" << endl;
    cout << endl << "[" << asctime(local) << "start to simulate reads for each fa" << endl;
    for(int tmp = 0; tmp < geneFaPathVec.size(); tmp++)
    {
    	//if(tmp >= 10)
    	//	break;
		string tmpGeneFa = geneFaPathVec[tmp];
		string tmpGeneId = geneIdVec[tmp];
		int tmpGeneCount = geneCountVec[tmp];
		log_ofs << "tmpGeneIndex:\t" << tmp << endl; 
		if(tmpGeneCount == 0)
		{
			log_ofs << "tmpGeneCount:\t0" << endl;
			log_ofs << "tmpGeneFa:\t" << tmpGeneFa << endl;
			log_ofs << "tmpGeneId:\t" << tmpGeneId << endl;
			continue;
		}
		else
		{
			log_ofs << "tmpGeneCount:\t" << tmpGeneCount << endl;
			log_ofs << "tmpGeneFa:\t" << tmpGeneFa << endl;
			log_ofs << "tmpGeneId:\t" << tmpGeneId << endl;
		}
		string tmpTargetFa = outputDir_rawSimulatedFa + int_to_str(tmp+1) + "_" + tmpGeneId;
		string tmp_cmd_simulate_read = inputArtBin + "art_illumina -ss " + opt_sequencing_platform
			+ " -na -l " + opt_length + " -i " + tmpGeneFa + " -c " + int_to_str(tmpGeneCount)
			+ " -o " + tmpTargetFa + ".";
		//cout << "tmp_cmd_simulate_read: " << tmp_cmd_simulate_read << endl;
		system(tmp_cmd_simulate_read.c_str());
		string tmpFqFile_end1 = tmpTargetFa + ".1.fq";
		string tmpFqFile_end2 = tmpTargetFa + ".2.fq";
		cout << "tmpFqFile_end1: " << tmpFqFile_end1 << endl;
		cout << "tmpFqFile_end2: " << tmpFqFile_end2 << endl;
		
		ifstream tmpFq_1_ifs(tmpFqFile_end1.c_str());
		int tmpReadNum = 0;
		while(!tmpFq_1_ifs.eof())
		{
			string tmpId;
			getline(tmpFq_1_ifs, tmpId);
			if(tmpId == "")
				break;
			string tmpSeq, tmpComm, tmpQual;
			getline(tmpFq_1_ifs, tmpSeq);
			getline(tmpFq_1_ifs, tmpComm);
			getline(tmpFq_1_ifs, tmpQual);
			tmpReadNum ++;
			simu_1_ofs << "@" << tmp + 1 << "_" << tmpGeneId << "_" << tmpReadNum << "/1" 
				<< endl << tmpSeq << endl << tmpComm << endl << tmpQual << endl;
		}
		tmpFq_1_ifs.eof();
		
		ifstream tmpFq_2_ifs(tmpFqFile_end2.c_str());
		tmpReadNum = 0;
		while(!tmpFq_2_ifs.eof())
		{
			string tmpId;
			getline(tmpFq_2_ifs, tmpId);
			if(tmpId == "")
				break;
			string tmpSeq, tmpComm, tmpQual;
			getline(tmpFq_2_ifs, tmpSeq);
			getline(tmpFq_2_ifs, tmpComm);
			getline(tmpFq_2_ifs, tmpQual);
			tmpReadNum ++;
			simu_2_ofs << "@" << tmp + 1 << "_" << tmpGeneId << "_" << tmpReadNum << "/2" 
				<< endl << tmpSeq << endl << tmpComm << endl << tmpQual << endl;
		}
		tmpFq_2_ifs.eof();

		string rm_tmpSimuSeq_1 = "rm " + tmpFqFile_end1;
		string rm_tmpSimuSeq_2 = "rm " + tmpFqFile_end2;
		//log_ofs << "rm_tmpSimuSeq_1:\t" << rm_tmpSimuSeq_1 << endl;
		//log_ofs << "rm_tmpSimuSeq_2:\t" << rm_tmpSimuSeq_2 << endl; 
		system(rm_tmpSimuSeq_1.c_str());
		system(rm_tmpSimuSeq_2.c_str());
	}
	log_ofs.close();
	simu_1_ofs.close();
	simu_2_ofs.close();
	return 0;
}