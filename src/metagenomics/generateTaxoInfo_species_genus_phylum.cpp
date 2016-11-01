#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>
#include "../mps3Lib/read_block_test.h"
#include "../mps3Lib/otherFunc.h"
#include "../mps3Lib/index_info.h"
#include "general/chromosomeSeq_info_vec.h"
#include "general/species_info.h"
#include "general/genus_info.h"
#include "general/phylum_info.h"
time_t nowtime;
struct tm *local;

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputChromosomeSeq2taxoInfoFile output_species_file output_genus_file output_phylum_file" << endl;
		exit(1);
	}
	string inputChromosomeSeq2taxoInfoFile = argv[1];

	string output_species_file = argv[2];
	string output_genus_file = argv[3];
	string output_phylum_file = argv[4];
	
	cout << "start to initiate ChromosomeSeq_Info_Vec" << endl;
	ChromosomeSeq_Info_Vec tmpChrSeqInfoVec;
	tmpChrSeqInfoVec.initiate(inputChromosomeSeq2taxoInfoFile);

	cout << "start to initiate Spcies_Info" << endl;
	Species_Info tmpSpeciesInfo;
	tmpSpeciesInfo.initiate(tmpChrSeqInfoVec);

	cout << "start to initiate Genus_Info" << endl;
	Genus_Info tmpGenusInfo;
	tmpGenusInfo.initiate(tmpChrSeqInfoVec);

	cout << "start to initiate Species_Info" << endl;
	Phylum_Info tmpPhylumInfo;
	tmpPhylumInfo.initiate(tmpChrSeqInfoVec);	

	tmpSpeciesInfo.printSpeciesInfo(output_species_file);
	tmpGenusInfo.printGenusInfo(output_genus_file);
	tmpPhylumInfo.printPhylumInfo(output_phylum_file);
	return 0;
}