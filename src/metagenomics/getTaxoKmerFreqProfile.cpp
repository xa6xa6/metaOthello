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

#include "../metagenomics/general/reissuedGenomeID2TaxoID_info.h"
#include "../mps3Lib/read_block_test.h"
#include "../mps3Lib/otherFunc.h"
#include "../mps3Lib/index_info.h"

using namespace std;

void getJfPathVecFromJfFileListFile(
	vector<string>& JfPathVec, string& Jf2taxoInfoFile)
{
	ifstream Jf2taxoInfo_ifs(Jf2taxoInfoFile.c_str());
	while(!Jf2taxoInfo_ifs.eof())
	{
		string tmpStr;
		getline(Jf2taxoInfo_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpJfFile = tmpStr.substr(0, tabLoc);
		JfPathVec.push_back(tmpJfFile);
	}
	Jf2taxoInfo_ifs.close();
}

int main(int argc, char** argv)
{
	if(argc != 8)
	{
		cout << "Executable#0" << endl;
		cout << "inputJellyFishBin#1 " << endl;
		cout << "NCBIfullTaxoId2NameFile#2" << endl;
		cout << "InputJfFileListFile_withTaxoId#3" << endl;
		cout << "outputDir#4" << endl;
		cout << "Kmer_length#5" << endl; 
		cout << "threads_num#6" << endl;
		cout << "bf_size#7" << endl;
		exit(1);
	}
	cout << "new" << endl;
	string inputJellyFishBin = argv[1];
	inputJellyFishBin += "/";
	string NCBIfullTaxoId2NameFile = argv[2];
	string InputJfFileListFile_withTaxoId = argv[3];
	string outputDir = argv[4];
	outputDir += "/";
	string mkdir_output = "mkdir " + outputDir;
	system(mkdir_output.c_str());
	string output_log = outputDir + "log.txt";
	ofstream log_ofs(output_log.c_str());
	string Kmer_length_str = argv[5];
	int Kmer_length = atoi(Kmer_length_str.c_str());
	string threads_num_str = argv[6];
	int threads_num = atoi(threads_num_str.c_str());
	string bf_size_str = argv[7];
	//int bf_size = atoi(bf_size_str.c_str());		
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    //////// STEP 1: start to initiate_Jf2taxoInfoFile_NCBIfullTaxoId2NameFile ////////////////
    ///////////////////////////////////////////////////////////////////////////////////////  
    ///////////////////////////////////////////////////////////////////////////////////////
    cout << "start to initiate_Jf2taxoInfoFile_NCBIfullTaxoId2NameFile" << endl;
	ReissuedGenomeID2TaxoID_Info reissuedGenomeId2TaxoIdInfo;
	cout << "start to do reissuedGenomeId2TaxoIdInfo.initiate_Jf2taxoInfoFile_NCBIfullTaxoId2NameFile(" << endl;
	reissuedGenomeId2TaxoIdInfo.initiate_Jf2taxoInfoFile_NCBIfullTaxoId2NameFile(
		InputJfFileListFile_withTaxoId, NCBIfullTaxoId2NameFile);
	reissuedGenomeId2TaxoIdInfo.generate_reissuedTaxoIdVec_taxoIdNamePairVec();						
	string outputDir_taxoInfo = outputDir + "taxoInfo/";
	reissuedGenomeId2TaxoIdInfo.print_taxoInfoDir_Jflist_reissuedId_oriId_genome_taxo(
		outputDir_taxoInfo);

    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    //////////////////// STEP 2: start to get dump jf file to fa file /////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////  
    ///////////////////////////////////////////////////////////////////////////////////////
	cout << "start to getJfPathVecFromJfFileListFile" << endl;
	vector<string> JfPathVec;
	getJfPathVecFromJfFileListFile(JfPathVec, InputJfFileListFile_withTaxoId);
	
	cout << "start to get dumpFaPathVec" << endl;
	vector<string> dumpFaPathVec;
	string outputDir_dump_fa = outputDir + "dump_fa/";
	string mkdir_dump_fa = "mkdir " + outputDir_dump_fa;
	string outputFile_dump_fa = outputDir + "jf_to_dump_fa.fileList.txt";
	ofstream dump_fa_log_ofs(outputFile_dump_fa.c_str());
	system(mkdir_dump_fa.c_str());
	cout << "JfPathVec.size(): " << JfPathVec.size() << endl;
	for(int tmp = 0; tmp < JfPathVec.size(); tmp++)
	{
		string tmp_dump_fa_file = outputDir_dump_fa + int_to_str(tmp) + ".dump_fa";
		dumpFaPathVec.push_back(tmp_dump_fa_file);
		dump_fa_log_ofs << JfPathVec[tmp] << "\t" << tmp_dump_fa_file << endl;
		string tmp_cmd_dump_fa = inputJellyFishBin + "jellyfish dump "
			+ JfPathVec[tmp] + " > " + tmp_dump_fa_file;
		cout << "tmp_cmd_dump_fa: " << tmp_cmd_dump_fa << endl;
		system(tmp_cmd_dump_fa.c_str());
	}
	dump_fa_log_ofs.close();

    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////// STEP 3: start to get genome-level Kmer frequency profile ////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////  
    ///////////////////////////////////////////////////////////////////////////////////////
	string outputDir_genome_level = outputDir + "genome_level";
	string mkdir_genome_level = "mkdir " + outputDir_genome_level;
	system(mkdir_genome_level.c_str());
	
	// cat
	string cat_genome_dump_fa = "cat";
	for(int tmp = 0; tmp < dumpFaPathVec.size(); tmp++)
	{
		cat_genome_dump_fa += " ";
		cat_genome_dump_fa += dumpFaPathVec[tmp];
	}
	cat_genome_dump_fa += " > ";
	string merged_genome_dump_fa = outputDir + "genome_level.merged.dumpFa";
	cat_genome_dump_fa += merged_genome_dump_fa;
	cout << endl << "start to cat fa file for each genome and get merged genome.dumpFa" << endl;
	cout << "cat_genome_dump_fa: " << cat_genome_dump_fa << endl;
	system(cat_genome_dump_fa.c_str());

	// count
	string cmd_jellyFish_countCanoKmer_genome = inputJellyFishBin + "jellyfish count -o "
		+ outputDir + "genome_level.merged.dumpFa.jf -m " + Kmer_length_str + " -t " + threads_num_str
		+ " -s " + bf_size_str + " -C " + merged_genome_dump_fa;
	cout << endl << "start to do Kmer counting for merged genome.dumpFa" << endl;
	cout << "cmd_jellyFish_countCanoKmer_genome: " << cmd_jellyFish_countCanoKmer_genome << endl;
	system(cmd_jellyFish_countCanoKmer_genome.c_str());

	// dump
	string cmd_jellyFish_dump_to_KmerFreq_genome = inputJellyFishBin + "jellyfish dump -t -c -o " +
		outputDir + "genome_level.merged.dumpFa.KmerFreq " + outputDir + "genome_level.merged.dumpFa.jf";
	cout << endl << "start to dump jf in human readable Kmer Fa file" << endl;
	cout << "cmd_jellyFish_dump_to_KmerFreq_genome: " << cmd_jellyFish_dump_to_KmerFreq_genome << endl;
	system(cmd_jellyFish_dump_to_KmerFreq_genome.c_str());
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////// STEP 4: start to get other taxo-level Kmer frequency profile ////////////////
    ///////////////////////////////////////////////////////////////////////////////////////  
    ///////////////////////////////////////////////////////////////////////////////////////
    cout << "start to initiate jfFileVecVec" << endl;
	vector< vector<string> > jfFileVecVec_species;
	vector< vector<string> > jfFileVecVec_genus;
	vector< vector<string> > jfFileVecVec_phylum;

	int total_species_num = reissuedGenomeId2TaxoIdInfo.return_species_num();
	cout << "total_species_num: " << total_species_num << endl;
	for(int tmpSpeciesIndex = 0; tmpSpeciesIndex < total_species_num; tmpSpeciesIndex++)
	{
		vector<string> tmpVec;
		jfFileVecVec_species.push_back(tmpVec);
	}

	int total_genus_num = reissuedGenomeId2TaxoIdInfo.return_genus_num();
	cout << "total_genus_num: " << total_genus_num << endl;
	for(int tmpGenusIndex = 0; tmpGenusIndex < total_genus_num; tmpGenusIndex++)
	{
		vector<string> tmpVec;
		jfFileVecVec_genus.push_back(tmpVec);
	}

	int total_phylum_num = reissuedGenomeId2TaxoIdInfo.return_phylum_num();
	cout << "total_phylum_num: " << total_phylum_num << endl;
	for(int tmpPhylumIndex = 0; tmpPhylumIndex < total_phylum_num; tmpPhylumIndex++)
	{
		vector<string> tmpVec;
		jfFileVecVec_phylum.push_back(tmpVec);
	}	

	int total_genome_num = reissuedGenomeId2TaxoIdInfo.return_genome_num();
	cout << "total_genome_num: " << total_genome_num << endl;
	for(int tmpGenomeIndex = 0; tmpGenomeIndex < total_genome_num; tmpGenomeIndex++)
	{
		string tmpJfFile = JfPathVec[tmpGenomeIndex];
		cout << "tmp: " << tmpGenomeIndex << " : " << tmpJfFile << endl;
		int tmpReissuedSpeciesId = reissuedGenomeId2TaxoIdInfo.return_reissuedTaxoId_from_genomeId(tmpGenomeIndex + 1, 8);
		cout << "tmpReissuedSpeciesId: " << tmpReissuedSpeciesId << endl;
		jfFileVecVec_species[tmpReissuedSpeciesId].push_back(tmpJfFile);
		int tmpReissuedGenusId = reissuedGenomeId2TaxoIdInfo.return_reissuedTaxoId_from_genomeId(tmpGenomeIndex + 1, 7);
		cout << "tmpReissuedGenusId: " << tmpReissuedGenusId << endl;
		jfFileVecVec_genus[tmpReissuedGenusId].push_back(tmpJfFile);
		int tmpReissuedPhylumId = reissuedGenomeId2TaxoIdInfo.return_reissuedTaxoId_from_genomeId(tmpGenomeIndex + 1, 3);
		cout << "tmpReissuedPhylumId: " << tmpReissuedPhylumId << endl;
		jfFileVecVec_phylum[tmpReissuedPhylumId].push_back(tmpJfFile);
	}

	//////////////////////////// SPECIES ///////////////////////////////
	cout << "start to merge genome jf files to generate species level jf files" << endl;
	log_ofs << "start to merge genome jf files to generate species level jf files" << endl;
	string outputDir_species_level = outputDir + "species_level/";
	string mkdir_species_level = "mkdir " + outputDir_species_level;
	system(mkdir_species_level.c_str());
	
	// merge jf for each species
	vector<string> JfPathVec_species;
	vector<string> dumpFaPathVec_species;
	for(int tmpSpeciesIndex = 0; tmpSpeciesIndex < total_species_num; tmpSpeciesIndex++)
	{	
		// merge jfs of genomes from the same species
		string merge2jf_path = outputDir_species_level + int_to_str(tmpSpeciesIndex) + ".merged.jf";
		JfPathVec_species.push_back(merge2jf_path);
		string mergeJf_cmd;
		if(jfFileVecVec_species[tmpSpeciesIndex].size() == 0)
		{
			cout << "error! jfFileVecVec_species[tmpSpeciesIndex].size() == 0" << endl;
			exit(1);
		}
		else if(jfFileVecVec_species[tmpSpeciesIndex].size() == 1)
			mergeJf_cmd = "cp " + (jfFileVecVec_species[tmpSpeciesIndex])[0] + " " + merge2jf_path;
		else
		{	
			mergeJf_cmd = inputJellyFishBin + "jellyfish merge -o " + merge2jf_path;
			for(int tmpGenomeIndex = 0; tmpGenomeIndex < jfFileVecVec_species[tmpSpeciesIndex].size(); tmpGenomeIndex ++)
			{
				mergeJf_cmd += " ";
				mergeJf_cmd += (jfFileVecVec_species[tmpSpeciesIndex])[tmpGenomeIndex];
			}
		}
		cout << "tmp mergeJf_cmd: " << mergeJf_cmd << endl;
		log_ofs << "tmp mergeJf_cmd: " << mergeJf_cmd << endl;
		system(mergeJf_cmd.c_str());
		// print Kmer frequency in fasta format
		string jfDump2fa_path = outputDir_species_level + int_to_str(tmpSpeciesIndex) + ".merged.dumpFa";
		dumpFaPathVec_species.push_back(jfDump2fa_path);
		string cmd_dump_fa = inputJellyFishBin + "jellyfish dump " + merge2jf_path + " > " + jfDump2fa_path;
		cout << "tmp cmd_dump_fa: " << cmd_dump_fa << endl;
		log_ofs << "tmp cmd_dump_fa: " << cmd_dump_fa << endl;
		system(cmd_dump_fa.c_str());		
	}
	// cat
	string cat_species_dump_fa = "cat";
	for(int tmp = 0; tmp < dumpFaPathVec_species.size(); tmp++)
	{
		cat_species_dump_fa += " ";
		cat_species_dump_fa += dumpFaPathVec_species[tmp];
	}
	cat_species_dump_fa += " > ";
	string merged_species_dump_fa = outputDir + "species_level.merged.dumpFa";
	cat_species_dump_fa += merged_species_dump_fa;
	cout << endl << "start to cat fa file for each species and get merged species.dumpFa" << endl;
	cout << "cat_species_dump_fa: " << cat_species_dump_fa << endl;
	system(cat_species_dump_fa.c_str());

	// count
	string cmd_jellyFish_countCanoKmer_species = inputJellyFishBin + "jellyfish count -o "
		+ outputDir + "species_level.merged.dumpFa.jf -m " + Kmer_length_str + " -t " + threads_num_str
		+ " -s " + bf_size_str + " -C " + merged_species_dump_fa;
	cout << endl << "start to do Kmer counting for merged species.dumpFa" << endl;
	cout << "cmd_jellyFish_countCanoKmer_species: " << cmd_jellyFish_countCanoKmer_species << endl;
	system(cmd_jellyFish_countCanoKmer_species.c_str());

	// dump
	string cmd_jellyFish_dump_to_KmerFreq_species = inputJellyFishBin + "jellyfish dump -t -c -o " 
		+ outputDir + "species_level.merged.dumpFa.KmerFreq " + outputDir + "species_level.merged.dumpFa.jf";
	cout << endl << "start to dump jf in human readable Kmer Fa file" << endl;
	cout << "cmd_jellyFish_dump_to_KmerFreq_species: " << cmd_jellyFish_dump_to_KmerFreq_species << endl;
	system(cmd_jellyFish_dump_to_KmerFreq_species.c_str());

	////////////////////////////  GENUS  ///////////////////////////////
	cout << "start to generate genus level fa files" << endl;
	log_ofs << "start to generate genus level fa files" << endl;
	string outputDir_genus_level = outputDir + "genus_level/";
	string mkdir_genus_level = "mkdir " + outputDir_genus_level;
	system(mkdir_genus_level.c_str());

	// merge jf for each genus
	vector<string> JfPathVec_genus;
	vector<string> dumpFaPathVec_genus;
	for(int tmpGenusIndex = 0; tmpGenusIndex < total_genus_num; tmpGenusIndex++)
	{
		string merge2jf_path = outputDir_genus_level + int_to_str(tmpGenusIndex) + ".merged.jf";
		JfPathVec_genus.push_back(merge2jf_path);
		string mergeJf_cmd;

		if(jfFileVecVec_genus[tmpGenusIndex].size() == 0)
		{
			cout << "error! jfFileVecVec_genus[tmpGenusIndex].size() == 0" << endl;
			exit(1);
		}
		else if(jfFileVecVec_genus[tmpGenusIndex].size() == 1)
			mergeJf_cmd = "cp " + (jfFileVecVec_genus[tmpGenusIndex])[0] + " " + merge2jf_path;
		else
		{	
			mergeJf_cmd = inputJellyFishBin + "jellyfish merge -o " + merge2jf_path;
			for(int tmpGenomeIndex = 0; tmpGenomeIndex < jfFileVecVec_genus[tmpGenusIndex].size(); tmpGenomeIndex ++)
			{
				mergeJf_cmd += " ";
				mergeJf_cmd += (jfFileVecVec_genus[tmpGenusIndex])[tmpGenomeIndex];
			}
		}

		cout << "tmp mergeJf_cmd: " << mergeJf_cmd << endl;
		log_ofs << "tmp mergeJf_cmd: " << mergeJf_cmd << endl;
		system(mergeJf_cmd.c_str());
		// print Kmer frequency in fasta format
		string jfDump2fa_path = outputDir_genus_level + int_to_str(tmpGenusIndex) + ".merged.dumpFa";
		dumpFaPathVec_genus.push_back(jfDump2fa_path);
		string cmd_dump_fa = inputJellyFishBin + "jellyfish dump " + merge2jf_path + " > " + jfDump2fa_path;
		cout << "tmp cmd_dump_fa: " << cmd_dump_fa << endl;
		log_ofs << "tmp cmd_dump_fa: " << cmd_dump_fa << endl;
		system(cmd_dump_fa.c_str());
	}

	// cat
	string cat_genus_dump_fa = "cat";
	for(int tmp = 0; tmp < dumpFaPathVec_genus.size(); tmp++)
	{
		cat_genus_dump_fa += " ";
		cat_genus_dump_fa += dumpFaPathVec_genus[tmp];
	}
	cat_genus_dump_fa += " > ";
	string merged_genus_dump_fa = outputDir + "genus_level.merged.dumpFa";
	cat_genus_dump_fa += merged_genus_dump_fa;
	cout << endl << "start to cat fa file for each genus and get merged genus.dumpFa" << endl;
	cout << "cat_genus_dump_fa: " << cat_genus_dump_fa << endl;
	system(cat_genus_dump_fa.c_str());
	
	// count
	string cmd_jellyFish_countCanoKmer_genus = inputJellyFishBin + "jellyfish count -o "
		+ outputDir + "genus_level.merged.dumpFa.jf -m " + Kmer_length_str + " -t " + threads_num_str
		+ " -s " + bf_size_str + " -C " + merged_genus_dump_fa;
	cout << endl << "start to do Kmer counting for merged genus.dumpFa" << endl;
	cout << "cmd_jellyFish_countCanoKmer_genus: " << cmd_jellyFish_countCanoKmer_genus << endl;
	system(cmd_jellyFish_countCanoKmer_genus.c_str());

	// dump
	string cmd_jellyFish_dump_to_KmerFreq_genus = inputJellyFishBin + "jellyfish dump -t -c -o " 
		+ outputDir + "genus_level.merged.dumpFa.KmerFreq " + outputDir + "genus_level.merged.dumpFa.jf";
	cout << endl << "start to dump jf in human readable Kmer Fa file" << endl;
	cout << "cmd_jellyFish_dump_to_KmerFreq_genus: " << cmd_jellyFish_dump_to_KmerFreq_genus << endl;
	system(cmd_jellyFish_dump_to_KmerFreq_genus.c_str());

	////////////////////////////  PHYLUM ///////////////////////////////
	cout << "start to generate phylum level fa files" << endl;
	log_ofs << "start to generate phylum level fa files" << endl;
	string outputDir_phylum_level = outputDir + "phylum_level/";
	string mkdir_phylum_level = "mkdir " + outputDir_phylum_level;
	system(mkdir_phylum_level.c_str());
	
	// merge jf for each phylum
	vector<string> JfPathVec_phylum;
	vector<string> dumpFaPathVec_phylum;
	for(int tmpPhylumIndex = 0; tmpPhylumIndex < total_phylum_num; tmpPhylumIndex++)
	{
		string merge2jf_path = outputDir_phylum_level + int_to_str(tmpPhylumIndex) + ".merged.jf";
		JfPathVec_phylum.push_back(merge2jf_path);
		string mergeJf_cmd;
		if(jfFileVecVec_phylum[tmpPhylumIndex].size() == 0)
		{
			cout << "error! jfFileVecVec_phylum[tmpPhylumIndex].size() == 0" << endl;
			exit(1);
		}
		else if(jfFileVecVec_phylum[tmpPhylumIndex].size() == 1)
			mergeJf_cmd = "cp " + (jfFileVecVec_phylum[tmpPhylumIndex])[0] + " " + merge2jf_path;
		else
		{
			mergeJf_cmd = inputJellyFishBin + "jellyfish merge -o " + merge2jf_path;
			for(int tmpGenomeIndex = 0; tmpGenomeIndex < jfFileVecVec_phylum[tmpPhylumIndex].size(); tmpGenomeIndex ++)
			{
				mergeJf_cmd += " ";
				mergeJf_cmd += (jfFileVecVec_phylum[tmpPhylumIndex])[tmpGenomeIndex];
			}
		}

		cout << "tmp mergeJf_cmd: " << mergeJf_cmd << endl;
		log_ofs << "tmp mergeJf_cmd: " << mergeJf_cmd << endl;
		system(mergeJf_cmd.c_str());
		// print Kmer frequency in fasta format
		string jfDump2fa_path = outputDir_phylum_level + int_to_str(tmpPhylumIndex) + ".merged.dumpFa";
		dumpFaPathVec_phylum.push_back(jfDump2fa_path);
		string cmd_dump_fa = inputJellyFishBin + "jellyfish dump " + merge2jf_path + " > " + jfDump2fa_path;
		cout << "tmp cmd_dump_fa: " << cmd_dump_fa << endl;
		log_ofs << "tmp cmd_dump_fa: " << cmd_dump_fa << endl;
		system(cmd_dump_fa.c_str());		
	}

	// cat
	string cat_phylum_dump_fa = "cat";
	for(int tmp = 0; tmp < dumpFaPathVec_phylum.size(); tmp++)
	{
		cat_phylum_dump_fa += " ";
		cat_phylum_dump_fa += dumpFaPathVec_phylum[tmp];
	}
	cat_phylum_dump_fa += " > ";
	string merged_phylum_dump_fa = outputDir + "phylum_level.merged.dumpFa";
	cat_phylum_dump_fa += merged_phylum_dump_fa;
	cout << endl << "start to cat fa file for each phylum and get merged phylum.dumpFa" << endl;
	cout << "cat_phylum_dump_fa: " << cat_phylum_dump_fa << endl;
	system(cat_phylum_dump_fa.c_str());
	
	// count
	string cmd_jellyFish_countCanoKmer_phylum = inputJellyFishBin + "jellyfish count -o "
		+ outputDir + "phylum_level.merged.dumpFa.jf -m " + Kmer_length_str + " -t " + threads_num_str
		+ " -s " + bf_size_str + " -C " + merged_phylum_dump_fa;
	cout << endl << "start to do Kmer counting for merged phylum.dumpFa" << endl;
	cout << "cmd_jellyFish_countCanoKmer_phylum: " << cmd_jellyFish_countCanoKmer_phylum << endl;
	system(cmd_jellyFish_countCanoKmer_phylum.c_str());

	// dump
	string cmd_jellyFish_dump_to_KmerFreq_phylum = inputJellyFishBin + "jellyfish dump -t -c -o "
		+ outputDir + "phylum_level.merged.dumpFa.KmerFreq " + outputDir + "phylum_level.merged.dumpFa.jf";
	cout << endl << "start to dump jf in human readable Kmer Fa file" << endl;
	cout << "cmd_jellyFish_dump_to_KmerFreq_phylum: " << cmd_jellyFish_dump_to_KmerFreq_phylum << endl;
	system(cmd_jellyFish_dump_to_KmerFreq_phylum.c_str());

	cout << "All jobs done!" << endl;
	log_ofs << "All jobs done!" << endl;
	log_ofs.close();
	return 0;
}