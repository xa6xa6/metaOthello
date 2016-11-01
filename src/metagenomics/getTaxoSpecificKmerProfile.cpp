// InputKmerSetFileListFile_withTaxoId
// 1. Kmer_file 
// 2. genome_id (9)
// 3. species_id (8)
// 4. genus_id (7)
// 5. family_id (6)
// 6. order_id (5)
// 7. class_id (4)
// 8. phylum_id (3)
// 9. kindom_id (2)
// 10. domain_id (1)
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
#include "general/reissuedGenomeID2TaxoID_info.h"
#include "../general/sortKmerBin_info.h"
#include "../Ye_implementation/othello.h"
#include "../Ye_implementation/mulothindex.h"
#include <chrono>
#include <inttypes.h>
#include "../Ye_implementation/io_helper.h"
#include "../mps3Lib/read_block_test.h"
#include "../mps3Lib/otherFunc.h"
#include "../mps3Lib/index_info.h"
using namespace std;

typedef unsigned long long keyT;
typedef uint64_t valueT;
typedef uint8_t freqT;
typedef uint8_t assignedFlagT;

keyT convertKmer2Key(char*s, int kmerLength)
{
    int tmpLength = 0;
    switch (*s) 
    {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            keyT ret = 0;
            while (*s == 'A' || *s == 'C' || *s =='T' || *s =='G')
            {
                ret <<=2;
                switch (*s) 
                {
                    case 'T':
                        ret ++;
                    case 'G':
                        ret ++;
                    case 'C':
                        ret ++;
                }
                s++;
                tmpLength ++;
                if(tmpLength == kmerLength)
                    return ret;
            }
    }
}

void KmerFreqUnionFile2freqArray(string& inputKmerFreqUnionFile, freqT* freqArray,
	IOHelper<keyT,valueT> *helper, MulOthIndex<keyT>* moth, int Kmer_length)
{
	ifstream KmerFreqUnion_ifs(inputKmerFreqUnionFile.c_str());
	unsigned long long tmpLineNO = 0;
    while(!KmerFreqUnion_ifs.eof())
    {
    	string tmpStr;
    	getline(KmerFreqUnion_ifs, tmpStr);
    	if(tmpStr == "")
    		break;
        tmpLineNO ++;
        unsigned long long tmpThousandIndex = tmpLineNO / 1000000;
        if(tmpLineNO == tmpThousandIndex * 1000000)          
            cout << "Processed Line #: " << tmpLineNO << endl; 
    	string tmpKmerStr = tmpStr.substr(0, Kmer_length);
	    char tmpKmerCharArray[Kmer_length];
		std::copy(tmpKmerStr.begin(), tmpKmerStr.end(), tmpKmerCharArray);
		keyT tmpKeyT = convertKmer2Key(tmpKmerCharArray, Kmer_length);
		valueT tmpValueT = moth->query(tmpKeyT);
    	string tmpFreqStr = tmpStr.substr(Kmer_length + 1);
    	unsigned long long tmpFreq = atoll(tmpFreqStr.c_str());
        if(tmpFreq <= 255)
 		   	freqArray[tmpValueT] = tmpFreq;
    	else
    		freqArray[tmpValueT] = 255;
    }
    KmerFreqUnion_ifs.close();
}

void getKmerFilePathVecFromKmerFileListFile(
	vector<string>& KmerFilePathVec, string& Kmer2taxoInfoFile)
{
	ifstream Kmer2taxoInfo_ifs(Kmer2taxoInfoFile.c_str());
	while(!Kmer2taxoInfo_ifs.eof())
	{
		string tmpStr;
		getline(Kmer2taxoInfo_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpKmerFile = tmpStr.substr(0, tabLoc);
		KmerFilePathVec.push_back(tmpKmerFile);
	}
	Kmer2taxoInfo_ifs.close();
}

void getKmerFreqUnionFilePathVecFromKmerFreqUnionFileListFile(
	vector<string>& KmerFreqUnionFilePathVec, string& KmerFreqUnionFileListFile)
{
	ifstream KmerFreqUnionFileList_ifs(KmerFreqUnionFileListFile.c_str());
	while(!KmerFreqUnionFileList_ifs.eof())
	{
		string tmpStr;
		getline(KmerFreqUnionFileList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		KmerFreqUnionFilePathVec.push_back(tmpStr);
	}
	KmerFreqUnionFileList_ifs.close();
}

int main(int argc, char** argv)
{
	if(argc != 9)
	{
		cout << "Executable#0" << endl;
		cout << "inputKmerIndexBitArrayFile#1" << endl;
		cout << "NCBIfullTaxoId2NameFile#2" << endl;
		cout << "InputKmerSetFileListFile_withTaxoId#3" << endl;
		cout << "inputKmerFreqUnionFileListFile#4" << endl;
		cout << "Kmer_length#5" << endl;
		cout << "split_bit_num#6" << endl;
		cout << "Kmer_num#7" << endl;
		cout << "outputDir#8" << endl;
		exit(1);
	}
	string NCBIfullTaxoId2NameFile = argv[2];
	string InputKmerSetFileListFile_withTaxoId = argv[3];
	string inputKmerFreqUnionFileListFile = argv[4];
	string KmerLengthStr = argv[5];
	string splitbitStr = argv[6];
    string KmerNumMaxStr = argv[7];
	string outputDir = argv[8];


	// KmerFreqUnionFileVec[0]--phylum; KmerFreqUnionFileVec[1]--genus; 
	// KmerFreqUnionFileVec[2]--species; KmerFreqUnionFileVec[3]--genome	
	vector<string> KmerFreqUnionFileVec; 
	getKmerFreqUnionFilePathVecFromKmerFreqUnionFileListFile(
		KmerFreqUnionFileVec, inputKmerFreqUnionFileListFile);

	int taxo_level_num = KmerFreqUnionFileVec.size(); // 
	if(taxo_level_num != 4)
	{
		cout << "currently, only 4 taxo ranks are supported: phylum, genus, species and  genome" << endl;
		exit(1);
	}

	cout << "start to create output dir" << endl;
	outputDir += "/";
	string cmd_mkdir = "mkdir " + outputDir;
	system(cmd_mkdir.c_str());
	string output_stats = outputDir + "stats.txt";
	string output_log = outputDir + "log.txt";
	ofstream stats_ofs(output_stats.c_str());
	ofstream log_ofs(output_log.c_str());

	string Kmer2taxoInfoFile = InputKmerSetFileListFile_withTaxoId;
	ReissuedGenomeID2TaxoID_Info reissuedGenomeId2TaxoIdInfo;
	reissuedGenomeId2TaxoIdInfo.initiate_Kmer2taxoInfoFile_NCBIfullTaxoId2NameFile(Kmer2taxoInfoFile, NCBIfullTaxoId2NameFile);
	reissuedGenomeId2TaxoIdInfo.generate_reissuedTaxoIdVec_taxoIdNamePairVec();						
	string outputDir_taxoInfo = outputDir + "taxoInfo";
	reissuedGenomeId2TaxoIdInfo.print_taxoInfoDir_KmerSet_reissuedId_oriId_genome_taxo(outputDir_taxoInfo);

	int total_genome_num = reissuedGenomeId2TaxoIdInfo.return_genome_num();
	int total_species_num = reissuedGenomeId2TaxoIdInfo.return_species_num();
	int total_genus_num = reissuedGenomeId2TaxoIdInfo.return_genus_num();
	int total_phylum_num = reissuedGenomeId2TaxoIdInfo.return_phylum_num();
	cout << "total_genome_num:\t" << total_genome_num << endl << "total_species_num:\t" << total_species_num << endl 
		<< "total_genus_num:\t" << total_genus_num << endl << "total_phylum_num:\t" << total_phylum_num << endl;
	stats_ofs << "total_genome_num:\t" << total_genome_num << endl << "total_species_num:\t" << total_species_num 
		<< endl << "total_genus_num:\t" << total_genus_num << endl << "total_phylum_num:\t" << total_phylum_num << endl;		

	
	int Kmer_length = atoi(KmerLengthStr.c_str());
    int splitbit = atoi(splitbitStr.c_str());
    unsigned long long int KmerNumMax = atoll(KmerNumMaxStr.c_str());

    cout << "start to load union kmer index bit array file" << endl;
    log_ofs << "start to load union kmer index bit array file" << endl;
    IOHelper<keyT,valueT> *helper;
    helper = new ConstantLengthKmerHelper<keyT,valueT>(Kmer_length,splitbit);
    MulOthIndex<keyT> * moth;
    moth = new MulOthIndex<keyT>(argv[1], helper);

    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ////////////// STEP 1: get Kmer Frequence array at each taxo level ////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////  
    ///////////////////////////////////////////////////////////////////////////////////////
    cout << "start to load union kmer frequency file and output repetitive Kmer" << endl;
    log_ofs << "start to load union kmer frequency file and output repetitive Kmer" << endl;
	// freqArrayVec[0]--phylum; freqArrayVec[1]--genus; 
	// freqArrayVec[2]--species; freqArrayVec[3]--genome	    
    vector<freqT*> freqArrayVec; 
    for(int tmp = 0; tmp < taxo_level_num; tmp++)
	{
		cout << "tmp taxo_level: " << tmp << endl;
		log_ofs << "tmp taxo_level: " << tmp << endl;
	    freqT* tmpFreqArray = new freqT[KmerNumMax];
	    KmerFreqUnionFile2freqArray(KmerFreqUnionFileVec[tmp], tmpFreqArray, helper, moth, Kmer_length);
	    freqArrayVec.push_back(tmpFreqArray);
	}

	cout << "start to set Kmer assignedOrNot Flag" << endl;
	log_ofs << "start to set Kmer assignedOrNot Flag" << endl;
	//bitset<KmerNumMax> bitvec;
	assignedFlagT* KmerAssignedFlagArray = new assignedFlagT[KmerNumMax];

	cout << "start to initiate allKmerFilePathVec" << endl;
	log_ofs << "start to initiate allKmerFilePathVec" << endl;
	string outputDir_KmerClass = outputDir + "KmerClass/";
	string mkdir_KmerClass = "mkdir " + outputDir_KmerClass;
	system(mkdir_KmerClass.c_str());
	int repetitiveClassId = total_genome_num + total_species_num + total_genus_num + total_phylum_num + 1;
	string outputFile_repetitiveKmerClass = outputDir_KmerClass + int_to_str(repetitiveClassId) + ".Kmer";
	ofstream repetitiveKmerClass_ofs(outputFile_repetitiveKmerClass.c_str());

    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    //////// STEP 2: scan each genome Kmer file to get taxo-specific Kmers ////////////////
    ///////////////////////////////////////////////////////////////////////////////////////  
    ///////////////////////////////////////////////////////////////////////////////////////
	vector<string> KmerFilePathVec;
	getKmerFilePathVecFromKmerFileListFile(KmerFilePathVec, Kmer2taxoInfoFile);
	int KmerSetNum = KmerFilePathVec.size();
	for(int tmp = 0; tmp < KmerSetNum; tmp++)
	{
		string tmpKmerSetPath = KmerFilePathVec[tmp];
		int tmpKmerSetSpecificClassId_genome, tmpKmerSetSpecificClassId_species,
			tmpKmerSetSpecificClassId_genus, tmpKmerSetSpecificClassId_phylum;
		reissuedGenomeId2TaxoIdInfo.get_KmerSetSpecificClass(tmp,
			tmpKmerSetSpecificClassId_genome, tmpKmerSetSpecificClassId_species,
			tmpKmerSetSpecificClassId_genus, tmpKmerSetSpecificClassId_phylum);
		vector<int> tmpKmerSetSpecificClassIdVec;
		tmpKmerSetSpecificClassIdVec.push_back(tmpKmerSetSpecificClassId_phylum);
		tmpKmerSetSpecificClassIdVec.push_back(tmpKmerSetSpecificClassId_genus);
		tmpKmerSetSpecificClassIdVec.push_back(tmpKmerSetSpecificClassId_species);
		tmpKmerSetSpecificClassIdVec.push_back(tmpKmerSetSpecificClassId_genome);
		vector<int> tmpKmerFreqVec;
		ifstream tmpKmerSet_ifs(tmpKmerSetPath.c_str());
		string outputFile_tmpSetSpecificKmer = outputDir_KmerClass + int_to_str(tmp) + ".taxoSpecific.Kmer";
		ofstream tmpSetSpecificKmer_ofs(outputFile_tmpSetSpecificKmer.c_str());
		while(!tmpKmerSet_ifs.eof())
		{
			string tmpStr;
			getline(tmpKmerSet_ifs, tmpStr);
			if(tmpStr == "")
				break;
			string tmpKmerStr = tmpStr.substr(0, Kmer_length);
			char tmpSeqChar[Kmer_length];
			std::copy(tmpKmerStr.begin(), tmpKmerStr.end(), tmpSeqChar);
			keyT tmpKeyT = convertKmer2Key(tmpSeqChar, Kmer_length);
			valueT tmpValueT = moth->query(tmpKeyT);
			if(!KmerAssignedFlagArray[tmpValueT]) //not assigned
			{
				KmerAssignedFlagArray[tmpValueT] = 1;
				bool setSpecific_or_repetitive_bool = false;
				for(int tmpTaxoLevel = taxo_level_num - 1; tmpTaxoLevel >= 0; tmpTaxoLevel--)
				{
					if((freqArrayVec[tmpTaxoLevel])[tmpValueT] == 1)
					{
						int tmpKmerSetSpecificClassId = tmpKmerSetSpecificClassIdVec[tmpTaxoLevel];
						tmpSetSpecificKmer_ofs << tmpKmerStr << "\t" << tmpKmerSetSpecificClassId << endl;
						setSpecific_or_repetitive_bool = true;
						break;
					}
				}
				if(!setSpecific_or_repetitive_bool)
					repetitiveKmerClass_ofs << tmpKmerStr << "\t" << int_to_str(repetitiveClassId) << endl;
			}
		}
		tmpKmerSet_ifs.close();
		tmpSetSpecificKmer_ofs.close();
	}
    for(int tmp = 0; tmp < taxo_level_num; tmp++)
		delete freqArrayVec[tmp];
	delete KmerAssignedFlagArray;
	delete helper;
	delete moth;
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    //////// STEP 3: start to cat taxo-specific Kmer file of each genome// ////////////////
    ///////////////////////////////////////////////////////////////////////////////////////  
    ///////////////////////////////////////////////////////////////////////////////////////
	cout << "start to create 'cat' command" << endl;
	log_ofs << "start to create 'cat' command" << endl;
	string outputFile_taxoSpecificKmerFile = outputDir + "taxoSpecific.Kmer";
	string cat_taxoSpecificKmerFile = "cat";
	for(int tmp = 0; tmp < KmerSetNum; tmp++)
	{
		string tmpSetSpecificKmerFile = outputDir_KmerClass + int_to_str(tmp) + ".taxoSpecific.Kmer";
		cat_taxoSpecificKmerFile += " ";
		cat_taxoSpecificKmerFile += tmpSetSpecificKmerFile;
	}
	cat_taxoSpecificKmerFile += " ";
	cat_taxoSpecificKmerFile += outputFile_repetitiveKmerClass;
	cat_taxoSpecificKmerFile += " > ";
	cat_taxoSpecificKmerFile += outputFile_taxoSpecificKmerFile;
	cout << "cat_taxoSpecificKmerFile: " << endl << cat_taxoSpecificKmerFile << endl;
	log_ofs << "cat_taxoSpecificKmerFile: " << endl << cat_taxoSpecificKmerFile << endl;
	system(cat_taxoSpecificKmerFile.c_str());
	repetitiveKmerClass_ofs.close();
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    // STEP 3: start to sort each catted taxo-specific kmer file according to prefix 3 char
    ///////////////////////////////////////////////////////////////////////////////////////  
    ///////////////////////////////////////////////////////////////////////////////////////
	cout << "start to sort catted taxo-specific kmer file according to prefix 3 char" << endl;
	log_ofs << "start to sort catted taxo-specific kmer file according to prefix 3 char" << endl;
	string outputDir_tmpSort = outputDir + "tmpSortDir";
	string outputFile_taxoSpecificKmerFile_sortedByPrefix3 = outputFile_taxoSpecificKmerFile + ".sortedByPrefix3";
	SortKmerBin_Info tmpSortKmerBinInfo;
	tmpSortKmerBinInfo.initiate_mkTmpDir(outputFile_taxoSpecificKmerFile, 
		outputFile_taxoSpecificKmerFile_sortedByPrefix3, outputDir_tmpSort);
	tmpSortKmerBinInfo.sortByPrefix3charOnly();	


	cout << "All jobs done!" << endl;
	log_ofs << "All jobs done!" << endl;
	stats_ofs.close();
	log_ofs.close();

	return 0;
}