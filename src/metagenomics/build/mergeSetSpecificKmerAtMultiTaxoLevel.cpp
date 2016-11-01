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
#include "../general/NCBIfullTaxoID2Name_info.h"
#include "../general/bacterialTaxo_info.h"

using namespace std;

time_t nowtime;
struct tm *local;

int getMergedKmerId(int species_num, int tmpId_species, int genus_num, int tmpId_genus,
	int family_num, int tmpId_family, int order_num, int tmpId_order, 
	int class_num, int tmpId_class, int phylum_num, int tmpId_phylum)
{
	int species_repetitive_id = species_num + 1;
	int genus_repetitive_id = genus_num + 1;
	int family_repetitive_id = family_num + 1;
	int order_repetitive_id = order_num + 1;
	int class_repetitive_id = class_num + 1;
	int phylum_repetitive_id = phylum_num + 1;
	if((tmpId_species < 1)||(tmpId_species > species_repetitive_id)
		||(tmpId_genus < 1)||(tmpId_genus > genus_repetitive_id)
		||(tmpId_family < 1)||(tmpId_family > family_repetitive_id)
		||(tmpId_order < 1)||(tmpId_order > order_repetitive_id)
		||(tmpId_class < 1)||(tmpId_class > class_repetitive_id)
		||(tmpId_phylum < 1)||(tmpId_phylum > phylum_repetitive_id))
	{
		cout << "((tmpId_species < 1)||(tmpId_species > species_repetitive_id)";
		cout <<	"||(tmpId_genus < 1)||(tmpId_genus > genus_repetitive_id)";
		cout <<	"||(tmpId_family < 1)||(tmpId_family > family_repetitive_id)";
		cout <<	"||(tmpId_order < 1)||(tmpId_order > order_repetitive_id)";
		cout <<	"||(tmpId_class < 1)||(tmpId_class > class_repetitive_id)";
		cout <<	"||(tmpId_phylum < 1)||(tmpId_phylum > phylum_repetitive_id))" << endl << endl;

		cout << "tmpId_species: " << tmpId_species << endl << "tmpId_genus: " << tmpId_genus << endl
			<< "tmpId_family: " << tmpId_family << endl << "tmpId_order: " << tmpId_order << endl
			<< "tmpId_class: " << tmpId_class << endl << "tmpId_phylum: " << tmpId_phylum << endl;
		exit(1);
	}

	if(tmpId_species < species_repetitive_id) // species specific
		return tmpId_species;
	else if(tmpId_genus < genus_repetitive_id) // genus specific 
		return tmpId_genus + species_num;
	else if(tmpId_family < family_repetitive_id)
		return tmpId_family + genus_num + species_num;
	else if(tmpId_order < order_repetitive_id)
		return tmpId_order + family_num + genus_num + species_num;
	else if(tmpId_class < class_repetitive_id)
		return tmpId_class + order_num + family_num + genus_num + species_num;
	else if(tmpId_phylum < phylum_repetitive_id)
		return tmpId_phylum + class_num + order_num + family_num + genus_num + species_num;
	else
		return 1 + phylum_num + class_num + order_num + family_num + genus_num + species_num;
}

int main(int argc, char** argv)
{
	if(argc != 11)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputSpeciesTaxoIdFile" << endl;
		cout << "#2 NCBIfullTaxoId2NameFile" << endl;		
		cout << "#3 inputKmerClassFile_species" << endl;
		cout << "#4 inputKmerClassFile_genus" << endl;
		cout << "#5 inputKmerClassFile_family" << endl;
		cout << "#6 inputKmerClassFile_order" << endl;
		cout << "#7 inputKmerClassFile_class" << endl;
		cout << "#8 inputKmerClassFile_phylum" << endl;
		cout << "#9 outputDir" << endl;
		cout << "#10 Kmer_length" << endl;
		exit(1);
	}
	string outputDir = argv[9];
	outputDir += "/";
	string mkdir_dir = "mkdir " + outputDir;
	system(mkdir_dir.c_str());
	string log_file = outputDir + "log.txt";
	ofstream log_ofs(log_file.c_str());

    nowtime = time(NULL);
    local = localtime(&nowtime);        
    log_ofs << endl << "[" << asctime(local) << "start to do mergeSetSpecificKemrAtMultiTaxoLevel" << endl;
    cout << endl << "[" << asctime(local) << "start to do mergeSetSpecificKemrAtMultiTaxoLevel" << endl;
	log_ofs << "Command Line:" << endl;
	for(int tmp = 0; tmp < 10; tmp++)
		log_ofs << "#" << tmp << "\t" << argv[tmp] << endl;

	string taxo_dir = outputDir + "taxo_info/";
	string mkdir_taxo_dir = "mkdir " + taxo_dir;
	system(mkdir_taxo_dir.c_str());

    nowtime = time(NULL);
    local = localtime(&nowtime);        
    log_ofs << endl << "[" << asctime(local) << "start to initiate bacterialTaxoInfo" << endl;
    cout << endl << "[" << asctime(local) << "start to initiate bacterialTaxoInfo" << endl;
	string inputSpeciesTaxoIdFile = argv[1];
	string NCBIfullTaxoId2NameFile = argv[2];
	BacterialTaxo_Info bacterialTaxoInfo;
	bacterialTaxoInfo.initiate_bacterialTaxoFile_NCBIfullTaxoId2NameFile(
		inputSpeciesTaxoIdFile, NCBIfullTaxoId2NameFile);
	bacterialTaxoInfo.reissueTaxoIdName_all();
	bacterialTaxoInfo.print(taxo_dir);
	int species_num = bacterialTaxoInfo.return_taxo_num(8);
	int genus_num = bacterialTaxoInfo.return_taxo_num(7);
	int family_num = bacterialTaxoInfo.return_taxo_num(6);
	int order_num = bacterialTaxoInfo.return_taxo_num(5);
	int class_num = bacterialTaxoInfo.return_taxo_num(4);
	int phylum_num = bacterialTaxoInfo.return_taxo_num(3);
	log_ofs << endl << "species_num: " << species_num << endl << "genus_num: " << genus_num << endl
		<< "family_num: " << family_num << endl << "order_num: " << order_num << endl
		<< "class_num: " << class_num << endl << "phylum_num: " << phylum_num << endl;

    nowtime = time(NULL);
    local = localtime(&nowtime);        
    log_ofs << endl << "[" << asctime(local) << "start to scan each Kmer Id file and get merged Kmer id" << endl;
    cout << endl << "[" << asctime(local) << "start to scan each Kmer Id file and get merged Kmer id" << endl;
	string inputKmerClassFile_species = argv[3];
	ifstream species_ifs(inputKmerClassFile_species.c_str());
	string inputKmerClassFile_genus = argv[4];
	ifstream genus_ifs(inputKmerClassFile_genus.c_str());
	string inputKmerClassFile_family = argv[5];
	ifstream family_ifs(inputKmerClassFile_family.c_str());
	string inputKmerClassFile_order = argv[6];
	ifstream order_ifs(inputKmerClassFile_order.c_str());
	string inputKmerClassFile_class = argv[7];
	ifstream class_ifs(inputKmerClassFile_class.c_str());
	string inputKmerClassFile_phylum = argv[8];
	ifstream phylum_ifs(inputKmerClassFile_phylum.c_str());					

	string outputMergedKmerSetFile = outputDir + "mergedKmerSetId.txt";
	string Kmer_length_str = argv[10];
	int Kmer_length = atoi(Kmer_length_str.c_str());
	ofstream mergedKmerSet_ofs(outputMergedKmerSetFile.c_str());
	unsigned long long tmpLineNO = 0;
	while((!species_ifs.eof())&&(!genus_ifs.eof())&&(!family_ifs.eof())&&
		(!order_ifs.eof())&&(!class_ifs.eof())&&(!phylum_ifs.eof()))
	{
		string tmpSpeciesStr, tmpGenusStr, tmpFamilyStr, 
			tmpOrderStr, tmpClassStr, tmpPhylumStr;
		getline(species_ifs, tmpSpeciesStr);
		getline(genus_ifs, tmpGenusStr);
		getline(family_ifs, tmpFamilyStr);
		getline(order_ifs, tmpOrderStr);
		getline(class_ifs, tmpClassStr);
		getline(phylum_ifs, tmpPhylumStr);
		if((tmpSpeciesStr == "")||(tmpGenusStr == "")||(tmpFamilyStr == "")
			||(tmpOrderStr == "")||(tmpClassStr == "")||(tmpPhylumStr == ""))							
		{
			cout << "((tmpSpeciesStr == "")||(tmpGenusStr == "")||(tmpFamilyStr == "")";
			cout <<	"||(tmpOrderStr == "")||(tmpClassStr == "")||(tmpPhylumStr == ""))" << endl;
			break;
		}

        tmpLineNO ++;
        unsigned long long tmpThousandIndex = tmpLineNO / 1000000;
        if(tmpLineNO == tmpThousandIndex * 1000000)          
            cout << "Processed Line #: " << tmpLineNO << endl;

		string tmpKmer_species = tmpSpeciesStr.substr(0, Kmer_length);
		int tmpId_species = atoi((tmpSpeciesStr.substr(Kmer_length + 1)).c_str());
		string tmpKmer_genus = tmpGenusStr.substr(0, Kmer_length);
		int tmpId_genus = atoi((tmpGenusStr.substr(Kmer_length + 1)).c_str());
		string tmpKmer_family = tmpFamilyStr.substr(0, Kmer_length);
		int tmpId_family = atoi((tmpFamilyStr.substr(Kmer_length + 1)).c_str());
		string tmpKmer_order = tmpOrderStr.substr(0, Kmer_length);
		int tmpId_order = atoi((tmpOrderStr.substr(Kmer_length + 1)).c_str());
		string tmpKmer_class = tmpClassStr.substr(0, Kmer_length);
		int tmpId_class = atoi((tmpClassStr.substr(Kmer_length + 1)).c_str());
		string tmpKmer_phylum = tmpPhylumStr.substr(0, Kmer_length);
		int tmpId_phylum = atoi((tmpPhylumStr.substr(Kmer_length + 1)).c_str());										

		string tmpKmerStr;
		if((tmpKmer_species != tmpKmer_genus)||(tmpKmer_species != tmpKmer_family)
			||(tmpKmer_species != tmpKmer_order)||(tmpKmer_species != tmpKmer_class)
			||(tmpKmer_species != tmpKmer_phylum))
		{
			cout << " ((tmpKmer_species != tmpKmer_genus)||(tmpKmer_species != tmpKmer_family)";
			cout << "||(tmpKmer_species != tmpKmer_order)||(tmpKmer_species != tmpKmer_class)";
			cout << "||(tmpKmer_species != tmpKmer_phylum)) " << endl;
			exit(1);
		}
		else
			tmpKmerStr = tmpKmer_species;
		int tmpMergedKmerId = getMergedKmerId(species_num, tmpId_species, genus_num, tmpId_genus,
			family_num, tmpId_family, order_num, tmpId_order, class_num, tmpId_class, phylum_num, tmpId_phylum);
		mergedKmerSet_ofs << tmpKmerStr << "\t" << tmpMergedKmerId << endl;
	}

	log_ofs << "All jobs done!" << endl;
	cout << "All jobs done!" << endl;
	mergedKmerSet_ofs.close();
	species_ifs.close();
	genus_ifs.close();
	family_ifs.close();
	order_ifs.close();
	class_ifs.close();
	phylum_ifs.close();
	log_ofs.close();
	return 0;
}