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

int main(int argc, char** argv)
{
	if((argc != 6)&&(argc != 7))
	{
		cout << "#0 Executable" << endl;
		cout << "#1 readLengthMin" << endl;
		cout << "#2 outputReadBaseNumFile" << endl;
		cout << "#3 SE_or_PE" << endl;
		cout << "#4 fa_or_fq" << endl;
		cout << "#5 seReadFile or peReadFile_1" << endl;
		cout << "#6 (peReadFile_2)" << endl;
		exit(1);
	}	
	string readLengthMinStr = argv[1];
	int readLengthMin = atoi(readLengthMinStr.c_str());
	string SE_or_PE_str = argv[3];
	string fa_or_fq_str = argv[4];
	bool SE_or_PE_bool;
	if(SE_or_PE_str == "SE")
		SE_or_PE_bool = true;
	else if(SE_or_PE_str == "PE")
		SE_or_PE_bool = false;
	else{
		cout << "invalid option for SE_or_PE_str: " << SE_or_PE_str << endl;
		exit(1);
	}

	if(SE_or_PE_bool && (argc != 6)){
		cout << "error! (SE_or_PE_bool && (argc != 6))" << endl;
		exit(1);
	}
	else if((!SE_or_PE_bool) && (argc != 7)){
		cout << "error! ((!SE_or_PE_bool) && (argc != 7))" << endl;
		exit(1);
	}
	else
	{}		

	bool fa_or_fq_bool;
	if((fa_or_fq_str == "fa")||(fa_or_fq_str == "Fa")||(fa_or_fq_str == "FA"))
		fa_or_fq_bool = true;
	else if((fa_or_fq_str == "fq")||(fa_or_fq_str == "Fq")||(fa_or_fq_str == "FQ"))
		fa_or_fq_bool = false;
	else{
		cout << "invalid option for fa_or_fq_str: " << fa_or_fq_str << endl;
		exit(1);
	}

	string seReadFile, peReadFile_1, peReadFile_2;
	if(SE_or_PE_bool)
		seReadFile = argv[5];
	else
	{
		peReadFile_1 = argv[5];
		peReadFile_2 = argv[6];
	}

	unsigned long long int total_base_num = 0;
	unsigned long long int total_read_num = 0;
	unsigned long long int total_read_num_tooShort = 0;
	if(SE_or_PE_bool)
	{
		ifstream SE_ifs(seReadFile.c_str());
		while(!SE_ifs.eof())
		{
			string tmpReadId;
			getline(SE_ifs, tmpReadId);
			if(tmpReadId == "")
				break;
			string tmpReadSeq;
			getline(SE_ifs, tmpReadSeq);
			int tmpReadSeqLength = tmpReadSeq.length();
			if(tmpReadSeqLength >= readLengthMin)
			{	
				total_base_num += tmpReadSeqLength;
				total_read_num ++;
			}
			else
				total_read_num_tooShort ++;
			if(!fa_or_fq_bool)
			{
				string tmpComment, tmpQual;
				getline(SE_ifs, tmpComment);
				getline(SE_ifs, tmpQual);
			}
		}
		SE_ifs.close();
	}
	else
	{
		ifstream PE_1_ifs(peReadFile_1.c_str());
		ifstream PE_2_ifs(peReadFile_2.c_str());
		while((!PE_1_ifs.eof())||(!PE_2_ifs.eof()))
		{
			string tmpReadId_1, tmpReadId_2;
			getline(PE_1_ifs, tmpReadId_1);
			getline(PE_2_ifs, tmpReadId_2);
			if((tmpReadId_1 == "")||(tmpReadId_2 == ""))
				break;
			string tmpReadSeq_1, tmpReadSeq_2;
			getline(PE_1_ifs, tmpReadSeq_1);
			getline(PE_2_ifs, tmpReadSeq_2);
			int tmpReadSeqLength_1 = tmpReadSeq_1.length();
			int tmpReadSeqLength_2 = tmpReadSeq_2.length();
			if(tmpReadSeqLength_1 + tmpReadSeqLength_2 >= readLengthMin)
			{	
				total_base_num += tmpReadSeqLength_1;
				total_base_num += tmpReadSeqLength_2;
				total_read_num ++;
			}
			else
				total_read_num_tooShort ++;
			if(!fa_or_fq_bool)
			{
				string tmpComment_1, tmpQual_1, tmpComment_2, tmpQual_2;
				getline(PE_1_ifs, tmpComment_1);
				getline(PE_1_ifs, tmpQual_1);
				getline(PE_2_ifs, tmpComment_2);
				getline(PE_2_ifs, tmpQual_2);				
			}
		}
		PE_1_ifs.close();
		PE_2_ifs.close();
	}
	string outputReadBaseNumFile = argv[2];
	ofstream readBaseNum_ofs(outputReadBaseNumFile.c_str());
	readBaseNum_ofs << "#1 outputReadBaseNumFile: " << outputReadBaseNumFile << endl;
	readBaseNum_ofs << "#2 SE_or_PE: " << SE_or_PE_str << endl;
	readBaseNum_ofs << "#3 fa_or_fq: " << fa_or_fq_str << endl;
	if(SE_or_PE_bool)
		readBaseNum_ofs << "#4 seReadFile: " << seReadFile << endl;
	else
	{
		readBaseNum_ofs << "#4 peReadFile_1: " << peReadFile_1 << endl;
		readBaseNum_ofs << "#5 peReadFile_2: " << peReadFile_2 << endl;
	}
	readBaseNum_ofs << endl;
	readBaseNum_ofs << "Too short read_num or pair_num: " << total_read_num_tooShort << endl << endl;
	if(SE_or_PE_bool)
		readBaseNum_ofs << "Total_read_num: " << total_read_num << endl;
	else
	{
		readBaseNum_ofs << "Total_readPair_num: " << total_read_num << endl;
		readBaseNum_ofs << "Total_read_num: " << total_read_num * 2 << endl;
	}
	readBaseNum_ofs << endl;
	readBaseNum_ofs << "Total_base_num: " << total_base_num << endl;
	readBaseNum_ofs.close();
	return 0;
}