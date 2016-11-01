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

typedef uint64_t keyT;
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

void generate_A_bit_vec(vector<keyT>& A_bit_vec)
{
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0); 
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0); 
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0);	
	A_bit_vec.push_back(0x0); 
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0); 
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0); 
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0); 
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0); 
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0); 
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0);
	A_bit_vec.push_back(0x0);
}

void generate_C_bit_vec(vector<keyT>& C_bit_vec)
{
	C_bit_vec.push_back(0x1);
	C_bit_vec.push_back(0x4);
	C_bit_vec.push_back(0x10);
	C_bit_vec.push_back(0x40);
	C_bit_vec.push_back(0x100);
	C_bit_vec.push_back(0x400);		
	C_bit_vec.push_back(0x1000);
	C_bit_vec.push_back(0x4000);
	C_bit_vec.push_back(0x10000);
	C_bit_vec.push_back(0x40000);
	C_bit_vec.push_back(0x100000);
	C_bit_vec.push_back(0x400000);	
	C_bit_vec.push_back(0x1000000);
	C_bit_vec.push_back(0x4000000);
	C_bit_vec.push_back(0x10000000);
	C_bit_vec.push_back(0x40000000);
	C_bit_vec.push_back(0x100000000);
	C_bit_vec.push_back(0x400000000);
	C_bit_vec.push_back(0x1000000000);
	C_bit_vec.push_back(0x4000000000);
	C_bit_vec.push_back(0x10000000000);
	C_bit_vec.push_back(0x40000000000);
	C_bit_vec.push_back(0x100000000000);
	C_bit_vec.push_back(0x400000000000);
	C_bit_vec.push_back(0x1000000000000);
	C_bit_vec.push_back(0x4000000000000);
	C_bit_vec.push_back(0x10000000000000);
	C_bit_vec.push_back(0x40000000000000);	
	C_bit_vec.push_back(0x100000000000000);
	C_bit_vec.push_back(0x400000000000000);	
	C_bit_vec.push_back(0x1000000000000000);					
}

void generate_G_bit_vec(vector<keyT>& G_bit_vec)
{
	G_bit_vec.push_back(0x2);
	G_bit_vec.push_back(0x8);
	G_bit_vec.push_back(0x20);
	G_bit_vec.push_back(0x80);
	G_bit_vec.push_back(0x200);
	G_bit_vec.push_back(0x800);
	G_bit_vec.push_back(0x2000);
	G_bit_vec.push_back(0x8000);			
	G_bit_vec.push_back(0x20000);
	G_bit_vec.push_back(0x80000);
	G_bit_vec.push_back(0x200000);
	G_bit_vec.push_back(0x800000);
	G_bit_vec.push_back(0x2000000);
	G_bit_vec.push_back(0x8000000);
	G_bit_vec.push_back(0x20000000);
	G_bit_vec.push_back(0x80000000);	
	G_bit_vec.push_back(0x200000000);
	G_bit_vec.push_back(0x800000000);
	G_bit_vec.push_back(0x2000000000);
	G_bit_vec.push_back(0x8000000000);
	G_bit_vec.push_back(0x20000000000);
	G_bit_vec.push_back(0x80000000000);
	G_bit_vec.push_back(0x200000000000);
	G_bit_vec.push_back(0x800000000000);
	G_bit_vec.push_back(0x2000000000000);
	G_bit_vec.push_back(0x8000000000000);
	G_bit_vec.push_back(0x20000000000000);
	G_bit_vec.push_back(0x80000000000000);
	G_bit_vec.push_back(0x200000000000000);
	G_bit_vec.push_back(0x800000000000000);
	G_bit_vec.push_back(0x2000000000000000);	
}

void generate_T_bit_vec(vector<keyT>& T_bit_vec)
{
	T_bit_vec.push_back(0x3);
	T_bit_vec.push_back(0xC);
	T_bit_vec.push_back(0x30);
	T_bit_vec.push_back(0xC0);
	T_bit_vec.push_back(0x300);
	T_bit_vec.push_back(0xC00);
	T_bit_vec.push_back(0x3000);
	T_bit_vec.push_back(0xC000);		
	T_bit_vec.push_back(0x30000);
	T_bit_vec.push_back(0xC0000);
	T_bit_vec.push_back(0x300000);
	T_bit_vec.push_back(0xC00000);
	T_bit_vec.push_back(0x3000000);
	T_bit_vec.push_back(0xC000000);
	T_bit_vec.push_back(0x30000000);
	T_bit_vec.push_back(0xC0000000);
	T_bit_vec.push_back(0x300000000);
	T_bit_vec.push_back(0xC00000000);		
	T_bit_vec.push_back(0x3000000000);
	T_bit_vec.push_back(0xC000000000);
	T_bit_vec.push_back(0x30000000000);
	T_bit_vec.push_back(0xC0000000000);
	T_bit_vec.push_back(0x300000000000);
	T_bit_vec.push_back(0xC00000000000);
	T_bit_vec.push_back(0x3000000000000);
	T_bit_vec.push_back(0xC000000000000);
	T_bit_vec.push_back(0x30000000000000);
	T_bit_vec.push_back(0xC0000000000000);		
	T_bit_vec.push_back(0x300000000000000);
	T_bit_vec.push_back(0xC00000000000000);
	T_bit_vec.push_back(0x3000000000000000);		
}

void generate_toAnd_bit_vec(vector<keyT>& toAnd_bit_vec)
{
	generate_T_bit_vec(toAnd_bit_vec);
}

string convertKmerKey2KmerSeq(keyT tmpKey, int Kmer_length, vector<keyT>& A_bit_vec,
	vector<keyT>& C_bit_vec, vector<keyT>& G_bit_vec, vector<keyT>& T_bit_vec, vector<keyT>& toAnd_bit_vec)
{
	//cout << "start to do convertKmerKey2KmerSeq" << endl;
	//cout << "tmpKey: " << tmpKey << endl;
	vector<string> tmpKmerCharVec;
	string tmpKmerSeq = "";
	for(int tmpLocInKmer = 0; tmpLocInKmer < Kmer_length; tmpLocInKmer ++)
	{
		//cout << "tmpLocInKmer: " << tmpLocInKmer << endl;
		keyT tmpKey_A = A_bit_vec[Kmer_length - tmpLocInKmer - 1];
		keyT tmpKey_C = C_bit_vec[Kmer_length - tmpLocInKmer - 1];
		keyT tmpKey_G = G_bit_vec[Kmer_length - tmpLocInKmer - 1];
		keyT tmpKey_T = T_bit_vec[Kmer_length - tmpLocInKmer - 1];
		keyT tmpKey_toAnd = toAnd_bit_vec[Kmer_length - tmpLocInKmer - 1];
		//cout << "tmpKey_A: " << tmpKey_A << endl;
		//cout << "tmpKey_C: " << tmpKey_C << endl;
		//cout << "tmpKey_G: " << tmpKey_G << endl;
		//cout << "tmpKey_T: " << tmpKey_T << endl;
		//cout << "tmpKey_toAnd" << tmpKey_toAnd << endl;
		if((tmpKey & tmpKey_toAnd) == tmpKey_A)
			tmpKmerSeq += "A";
		else if((tmpKey & tmpKey_toAnd) == tmpKey_C)
			tmpKmerSeq += "C";
		else if((tmpKey & tmpKey_toAnd) == tmpKey_G)
			tmpKmerSeq += "G";
		else
			tmpKmerSeq += "T";
	}
	return tmpKmerSeq;
}

bool KmerKeySetIdPairComparison(const pair<keyT, int>& a, const pair<keyT, int>& b)
{
	if(a.first < b.first)
		return true;
	else if(a.first == b.first)
	{
		if(a.second < b.second)
			return true;
		else
			return false;
	}
	else
		return false;
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "#0 Executable " << endl;
		cout << "#1 inputSortedMetagenomicsKmerSetIdFile" << endl;
		cout << "#2 inputSortedHumanKmerCountFile" << endl;
		cout << "#3 outputMergedKmerSetIdFile" << endl;
		cout << "#4 lastIdInMetagenomicsLothello(repetitiveSetId)" << endl;
		cout << "#5 Kmer_length" << endl;
		exit(1);
	}
	string inputSortedMetagenomicsKmerSetIdFile = argv[1];
	string inputSortedHumanKmerCountFile = argv[2];
	string outputMergedKmerSetIdFile = argv[3];
	string lastIdInMetagenomicsLothelloStr = argv[4];
	string Kmer_length_str = argv[5];
	cout << "inputSortedMetagenomicsKmerSetIdFile: " << inputSortedMetagenomicsKmerSetIdFile << endl;
	cout << "inputSortedHumanKmerCountFile: " << inputSortedHumanKmerCountFile << endl;
	cout << "outputMergedKmerSetIdFile: " << outputMergedKmerSetIdFile << endl;
	cout << "lastIdInMetagenomicsLothello: " << lastIdInMetagenomicsLothelloStr << endl;
	cout << "Kmer_length: " << Kmer_length_str << endl;

	cout << "start to load metagenomicsKmerSetId file and generate metagenomicsKmerKeySetIdPairVec" << endl;
	int Kmer_length = atoi(Kmer_length_str.c_str());
	ifstream metagenomicsKmerSetId_ifs(inputSortedMetagenomicsKmerSetIdFile.c_str());
	vector< pair<keyT, int> > metagenomicsKmerKeySetIdPairVec;
	while(!metagenomicsKmerSetId_ifs.eof())
	{
		string tmpStr;
		getline(metagenomicsKmerSetId_ifs, tmpStr);
		if(tmpStr == "")
			break;
		string tmpMetagenomicsKmer = tmpStr.substr(0, Kmer_length);
		char tmpMetagenomicsKmerCharArray[50];
		std::copy(tmpMetagenomicsKmer.begin(), tmpMetagenomicsKmer.end(), tmpMetagenomicsKmerCharArray);
		keyT tmpMetagenomicsKmerKey = convertKmer2Key(tmpMetagenomicsKmerCharArray, Kmer_length);
		string tmpIdStr = tmpStr.substr(Kmer_length + 1);
		int tmpId = atoi(tmpIdStr.c_str());
		metagenomicsKmerKeySetIdPairVec.push_back(pair<keyT,int>(tmpMetagenomicsKmerKey, tmpId));
	}
	metagenomicsKmerSetId_ifs.close();
	cout << "metagenomicsKmerKeySetIdPairVec.size(): " << metagenomicsKmerKeySetIdPairVec.size() << endl;

	cout << "start to load humanKmerCount file and generate humanKmerKeyVec" << endl;
	int lastIdInMetagenomicsLothello = atoi(lastIdInMetagenomicsLothelloStr.c_str());
	int humanSpecificSetId = lastIdInMetagenomicsLothello + 1;
	vector< pair<keyT, int> > humanKmerKeySetIdPairVec;
	ifstream humanKmerCount_ifs(inputSortedHumanKmerCountFile.c_str());
	while(!humanKmerCount_ifs.eof())
	{
		string tmpStr;
		getline(humanKmerCount_ifs, tmpStr);
		if(tmpStr == "")
			break;
		string tmpHumanKmer = tmpStr.substr(0, Kmer_length);
		char tmpHumanKmerCharArray[50];
		std::copy(tmpHumanKmer.begin(), tmpHumanKmer.end(), tmpHumanKmerCharArray);
		keyT tmpHumanKmerKey = convertKmer2Key(tmpHumanKmerCharArray, Kmer_length);
		humanKmerKeySetIdPairVec.push_back(pair<keyT,int>(tmpHumanKmerKey, humanSpecificSetId));
	}
	humanKmerCount_ifs.close();
	cout << "humanKmerKeyVec.size(): " << humanKmerKeySetIdPairVec.size() << endl;

	cout << "start to merge metagenomicsKmerKeySetIdPairVec and humanKmerKeyVec" << endl;
	cout << "start to copy metagenomicsKmerKeySetIdPairVec to mergedKmerKeySetIdPairVec" << endl;
	vector< pair<keyT, int> > mergedKmerKeySetIdPairVec;
	for(int tmp = 0; tmp < metagenomicsKmerKeySetIdPairVec.size(); tmp++)
		mergedKmerKeySetIdPairVec.push_back(metagenomicsKmerKeySetIdPairVec[tmp]);
	vector< pair<keyT, int> >().swap(metagenomicsKmerKeySetIdPairVec);

	cout << "start to copy humanKmerKeySetIdPairVec to mergedKmerKeySetIdPairVec" << endl;
	for(int tmp = 0; tmp < humanKmerKeySetIdPairVec.size(); tmp++)
		mergedKmerKeySetIdPairVec.push_back(humanKmerKeySetIdPairVec[tmp]);
	vector< pair<keyT, int> >().swap(humanKmerKeySetIdPairVec);

	cout << "start to sort mergedKmerKeySetIdPairVec based on kmer key" << endl;
	sort(mergedKmerKeySetIdPairVec.begin(), mergedKmerKeySetIdPairVec.end(), KmerKeySetIdPairComparison);

	cout << "start to initiate A_bit_vec, C_bit_vec, G_bit_vec, T_bit_vec" << endl;
	vector<keyT> A_bit_vec, C_bit_vec, G_bit_vec, T_bit_vec, toAnd_bit_vec;
	generate_A_bit_vec(A_bit_vec);
	generate_C_bit_vec(C_bit_vec);
	generate_G_bit_vec(G_bit_vec);
	generate_T_bit_vec(T_bit_vec);
	generate_toAnd_bit_vec(toAnd_bit_vec);

	cout << "start to output mergedKmerKeySetIdPairVec" << endl;
	// metagenomics(1, 2690) & specific --> remain the same
	// human & specific --> 2691
	// metagenomics & human shared --> 2692
	ofstream mergedKmerSetId_ofs(outputMergedKmerSetIdFile.c_str());
	int sharedKmerSetId = humanSpecificSetId + 1;
	keyT lastKmerKey = mergedKmerKeySetIdPairVec[0].first;
	int lastKmerSetId = mergedKmerKeySetIdPairVec[0].second;
	for(int tmp = 1; tmp < mergedKmerKeySetIdPairVec.size(); tmp++)
	{
		keyT tmpKmerKey = mergedKmerKeySetIdPairVec[tmp].first;
		//cout << "tmpKmerKey: " << tmpKmerKey << endl;
		if(tmpKmerKey != lastKmerKey)
		{
			//cout << "lastKmerKey: " << lastKmerKey << endl;
			string lastKmerSeq = convertKmerKey2KmerSeq(lastKmerKey, Kmer_length,
				A_bit_vec, C_bit_vec, G_bit_vec, T_bit_vec, toAnd_bit_vec);
			//cout << "lastKmerSeq: " << lastKmerSeq << endl;
			mergedKmerSetId_ofs << lastKmerSeq << "\t" << lastKmerSetId << endl;
			lastKmerKey = tmpKmerKey;
			lastKmerSetId = mergedKmerKeySetIdPairVec[tmp].second;
		}
		else
			lastKmerSetId = sharedKmerSetId;
	}
	cout << "lastKmerKey: " << lastKmerKey << endl;
	string lastKmerSeq = convertKmerKey2KmerSeq(lastKmerKey, Kmer_length,
		A_bit_vec, C_bit_vec, G_bit_vec, T_bit_vec, toAnd_bit_vec);
	cout << "lastKmerSeq: " << lastKmerSeq << endl;
	mergedKmerSetId_ofs << lastKmerSeq << "\t" << lastKmerSetId << endl;
	mergedKmerSetId_ofs.close();
	cout << "All jobs done!" << endl;
	return 0;
}