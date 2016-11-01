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
#include "../mps3Lib/read_block_test.h"
#include "../mps3Lib/otherFunc.h"
#include "../mps3Lib/index_info.h"
#include "general/bacterialTaxo_info.h"
#include "../general/sortKmerBin_info.h"
#include "../Ye_implementation/othello.h"
#include "../Ye_implementation/mulothindex.h"
#include <chrono>
#include <inttypes.h>
#include "../Ye_implementation/io_helper.h"
using namespace std;

typedef unsigned long long keyT;
typedef uint64_t valueT;
typedef uint16_t setIdT;
typedef uint8_t freqT;

//vector<keyT> keys;
//vector<valueT> values;

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

int main(int argc, char** argv)
{
	if(argc != 11)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputKmerIndexBitArrayFile" << endl;
		cout << "#2 inputSortedKmerUnionFile" << endl;
		cout << "#3 inputKmerDir" << endl;
		cout << "#4 outputDir" << endl;
		cout << "#5 KmerLength" << endl;
		cout << "#6 splitBit" << endl;
		cout << "#7 KmerNumMax" << endl;
		cout << "#8 inputSpeciesTaxoIdFile" << endl;
		cout << "#9 NCBIfullTaxoId2NameFile" << endl;
		cout << "#10 taxo_level_name" << endl;
		exit(1);
	}
    string inputKmerIndexBitArrayFile = argv[1];
	string inputSortedKmerUnionFile = argv[2];
	string inputKmerDir = argv[3];
	string outputDir = argv[4];
	string KmerLengthStr = argv[5];
	int Kmer_length = atoi(KmerLengthStr.c_str());
    string splitbitStr = argv[6];
    int splitbit = atoi(splitbitStr.c_str());
    string KmerNumMaxStr = argv[7];
    unsigned long long int KmerNumMax = atoll(KmerNumMaxStr.c_str());	
	string inputSpeciesTaxoIdFile = argv[8];
	string NCBIfullTaxoId2NameFile = argv[9];
	string taxo_level_name = argv[10];
	int taxo_rank;
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
	string output_log = outputDir + "log.txt";
	ofstream log_ofs(output_log.c_str());
	log_ofs << "inputKmerIndexBitArrayFile: " << endl << inputKmerIndexBitArrayFile << endl;
	log_ofs << "inputSortedKmerUnionFile: " << endl << inputSortedKmerUnionFile << endl;
	log_ofs << "inputKmerDir: " << endl << inputKmerDir << endl;
	log_ofs << "outputDir: " << endl << outputDir << endl;
	log_ofs << "Kmer_length: " << endl << Kmer_length << endl;
	log_ofs << "splitbit: " << endl << splitbit << endl;
	log_ofs << "KmerNumMax: " << endl << KmerNumMax << endl;
	log_ofs << "inputSpeciesTaxoIdFile: " << endl << inputSpeciesTaxoIdFile << endl;
	log_ofs << "NCBIfullTaxoId2NameFile: " << endl << NCBIfullTaxoId2NameFile << endl;
	log_ofs << "taxo_level_name: " << endl << taxo_level_name << endl;

	cout << "taxo_level_name: " << taxo_level_name << endl;
	cout << "taxo_rank: " << taxo_rank << endl;
	log_ofs << endl << "taxo_level_name: " << taxo_level_name << endl;
	log_ofs << "taxo_rank: " << taxo_rank << endl;

	BacterialTaxo_Info bacterialTaxoInfo;
	bacterialTaxoInfo.initiate_bacterialTaxoFile_NCBIfullTaxoId2NameFile(
		inputSpeciesTaxoIdFile, NCBIfullTaxoId2NameFile);
	bacterialTaxoInfo.reissueTaxoIdName_all();	
	string outputDir_taxoInfo = outputDir + "taxoInfo/";
	bacterialTaxoInfo.print(outputDir_taxoInfo);

    IOHelper<keyT,valueT> *helper;
    helper = new ConstantLengthKmerHelper<keyT,valueT>(Kmer_length,splitbit);
    cout << "start to load union kmer index bit array file" << endl;
    MulOthIndex<keyT> * moth;
    moth = new MulOthIndex<keyT>(argv[1], helper);
    cout << "start to load union kmer frequency file and output repetitive Kmer" << endl;
    freqT* freqArray = new freqT[KmerNumMax];
	
    cout << "start to initiate setIdArray" << endl;
   	setIdT* setIdArray = new setIdT[KmerNumMax];

   	cout << "start to load each species Kmer file" << endl;
   	int species_num = bacterialTaxoInfo.return_species_num();
    cout << "species_num: " << species_num << endl;
    for(uint16_t tmp_reissued_species_id = 0; tmp_reissued_species_id < (uint16_t)species_num; tmp_reissued_species_id ++)
    {
    	int tmp_species_id_ori = bacterialTaxoInfo.return_species_id((int)tmp_reissued_species_id);// tmp_species_index + 1;
    	string tmp_species_Kmer_file = inputKmerDir + "/" + int_to_str(tmp_species_id_ori) + ".Kmer";
        int tmp_taxo_reissued_id = bacterialTaxoInfo.return_taxo_reissued_id(taxo_rank, tmp_reissued_species_id);
        cout << "tmp_reissued_species_id: " << tmp_reissued_species_id << endl;
        cout << "tmp_taxo_reissued_id: " << tmp_taxo_reissued_id << endl;
    	cout << "tmpFileIndex: " << tmp_reissued_species_id << endl;
    	cout << "tmp_species_Kmer_file: " << tmp_species_Kmer_file << endl;
        log_ofs << "tmp_reissued_species_id: " << tmp_reissued_species_id << endl;
        log_ofs << "tmp_taxo_reissued_id: " << tmp_taxo_reissued_id << endl;
    	log_ofs << "tmpFileIndex: " << tmp_reissued_species_id << endl;
    	log_ofs << "tmp_species_Kmer_file: " << tmp_species_Kmer_file << endl;    	
    	ifstream tmpKmerSet_ifs(tmp_species_Kmer_file.c_str());
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
			if(freqArray[tmpValueT] == 0)
			{
				freqArray[tmpValueT] = 1;
				setIdArray[tmpValueT] = tmp_taxo_reissued_id + 1;
			}
			else if(freqArray[tmpValueT] == 1)
			{
                if(setIdArray[tmpValueT] == tmp_taxo_reissued_id + 1)
                {}
                else
                    freqArray[tmpValueT] = 255;
            }
			else
			{}
		}
    	tmpKmerSet_ifs.close();
    }

    int total_taxo_num = bacterialTaxoInfo.return_taxo_num(taxo_rank);
    int repetitive_class_id = total_taxo_num + 1;
    cout << "total_taxo_num: " << total_taxo_num << endl;
    cout << "repetitive_class_id: " << repetitive_class_id << endl;
    log_ofs << "total_taxo_num: " << total_taxo_num << endl;
    log_ofs << "repetitive_class_id: " << repetitive_class_id << endl;    
    string outputFile_repetitive = outputDir + "valid.repetitive.Kmer";
    ofstream repetitive_ofs(outputFile_repetitive.c_str());
    string outputFile_invalid = outputDir + "invalid.Kmer";
   	ofstream invalid_ofs(outputFile_invalid.c_str());
   	string outputFile_valid = outputDir + "valid.Kmer";
   	ofstream valid_ofs(outputFile_valid.c_str());
    ifstream sortedKmerUnion_ifs(inputSortedKmerUnionFile.c_str());
    unsigned long long tmpLineNO = 0; 
    while(!sortedKmerUnion_ifs.eof())
    {
    	string tmpStr;
    	getline(sortedKmerUnion_ifs, tmpStr);
    	if(tmpStr == "")
    		break;
        tmpLineNO ++;
        unsigned long long tmpThousandIndex = tmpLineNO / 1000000;
        if(tmpLineNO == tmpThousandIndex * 1000000)          
            cout << "Processed Line #: " << tmpLineNO << endl; 
    	string tmpKmerStr = tmpStr.substr(0, Kmer_length);
	    char tmpKmerCharArray[Kmer_length];
		std::copy(tmpKmerStr.begin(), tmpKmerStr.end(), tmpKmerCharArray);
		keyT tmpKey = convertKmer2Key(tmpKmerCharArray, Kmer_length);
		valueT tmpValue = moth->query(tmpKey);
        freqT tmpFreq = freqArray[tmpValue];
        setIdT tmpSetId = setIdArray[tmpValue];
        if(tmpFreq == 0)
            invalid_ofs << tmpKmerStr << "\t" << tmpSetId << "\t" << (int)tmpFreq << endl;
        else if(tmpFreq == 1)
        {
            if(((int)tmpSetId > total_taxo_num)||(tmpSetId == 0))
                invalid_ofs << tmpKmerStr << "\t" << tmpSetId << "\t1" << endl;
            else
                valid_ofs << tmpKmerStr << "\t" << tmpSetId << endl; 
        }
        else
        { 
            valid_ofs << tmpKmerStr << "\t" << repetitive_class_id << endl;
            repetitive_ofs << tmpKmerStr << "\t" << repetitive_class_id << endl;
        }
    }

    invalid_ofs.close();
    valid_ofs.close();
    repetitive_ofs.close();
    sortedKmerUnion_ifs.close();// _ifs.close();
    //repetitiveKmer_ofs.close();
    //alienKmer_ofs.close();
    cout << "All jobs done!" << endl;
    log_ofs << "All jobs done!" << endl;
    delete freqArray;
    delete helper;
    delete setIdArray;
	delete moth;
	log_ofs.close();
	//taxoInfo_ofs.close();
	return 0;
}