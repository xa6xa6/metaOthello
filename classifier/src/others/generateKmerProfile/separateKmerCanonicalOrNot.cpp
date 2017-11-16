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
#include <hash_map>
#include <map>
#include <set>
#include "../mps3Lib/read_block_test.h"
#include "../mps3Lib/otherFunc.h"
#include "../mps3Lib/index_info.h"

using namespace std;

void rcmCharArray(char* oriChar, char* targetChar, int length)
{
   	for(int tmp = 0; tmp < length; tmp++)
       	*(targetChar + tmp) = getCharRevComp(*(oriChar + length - 1 - tmp));
    targetChar[length] = '\0';
}

bool cmp_lexicoOrder(char* tmpChar_1, char* tmpChar_2, int len)
{
   	char tmp2_1[len + 1];
   	char tmp2_2[len + 1];
   	strncpy(tmp2_1, tmpChar_1, len);
   	strncpy(tmp2_2, tmpChar_2, len);
   	tmp2_1[len] = '\0';
   	tmp2_2[len] = '\0';
   	return (strcmp(tmp2_1, tmp2_2) <= 0);
}

bool KmerCanonicalOrNot(string& tmpSeq_ori)
{
	int seqLength = tmpSeq_ori.length();
	char charArray_ori[30];
	std::copy(tmpSeq_ori.begin(), tmpSeq_ori.end(), charArray_ori);
	char charArray_rcm[30];
	rcmCharArray(charArray_ori, charArray_rcm, seqLength);
	return cmp_lexicoOrder(charArray_ori, charArray_rcm, seqLength);
}


int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputRawKmerProfile outputFilePrefix" << endl;
		exit(1);
	}
	string inputRawKmerProfile = argv[1];
	string outputFilePrefix = argv[2];
	string outputFile_canonical = outputFilePrefix + ".canonical";
	string outputFile_nonCanonical = outputFilePrefix + ".noncanonical";
	ifstream raw_ifs(inputRawKmerProfile.c_str());
	ofstream canonical_ofs(outputFile_canonical.c_str());
	ofstream noncanonical_ofs(outputFile_nonCanonical.c_str());
	while(!raw_ifs.eof())
	{
		string tmpStr;
		getline(raw_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpStr_kmer = tmpStr.substr(0, tabLoc);
		if(KmerCanonicalOrNot(tmpStr_kmer))
			canonical_ofs << tmpStr << endl;
		else
			noncanonical_ofs << tmpStr << endl;
	}
	raw_ifs.close();
	canonical_ofs.close();
	noncanonical_ofs.close();
	return 0;
}