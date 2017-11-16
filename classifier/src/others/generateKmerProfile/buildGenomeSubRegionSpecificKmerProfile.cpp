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
time_t nowtime;
struct tm *local;

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
	if(argc != 8)
	{
		cout << "Executable inputIndexFolderPath inputClassificationGroupNum outputFolder KmerLength JellyFishBin thread_num incompleteKmerCombinationOrNot" << endl;
		exit(1);
	}
	string incompleteKmerCombinationOrNot_str = argv[7];
	bool incompleteKmerCombinationOrNot_bool;
	if(incompleteKmerCombinationOrNot_str == "Y")
		incompleteKmerCombinationOrNot_bool = true;
	else if(incompleteKmerCombinationOrNot_str == "N")
		incompleteKmerCombinationOrNot_bool = false;
	else
	{
		cout << "invalid incompleteKmerCombinationOrNot option, should be Y or N" << endl;
		exit(1);
	}

	string thread_num_str = argv[6];
	string inputJellyFishBinPath = argv[5];
	inputJellyFishBinPath += "/";
	string inputIndexFolderPath = argv[1];
	string inputClassificationGroupNumStr = argv[2];
	int inputClassificationGroupNum = atoi(inputClassificationGroupNumStr.c_str());
	int genomeSubRegionNum = inputClassificationGroupNum - 2; // shared kmer group & alien kmer group
	string KmerLengthStr = argv[4];
	int KmerLength = atoi(KmerLengthStr.c_str());
	string outputDirStr = argv[3];
	outputDirStr += "/";
  	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());
   	string log_file = outputDirStr + "log.txt";
   	ofstream log_ofs(log_file.c_str());

	string outputDirStr_fa = outputDirStr + "fa/";
	string cmd_mkdir_fa = "mkdir " + outputDirStr_fa;
	system(cmd_mkdir_fa.c_str());
	string outputDirStr_jf = outputDirStr + "jf/";
	string cmd_mkdir_jf = "mkdir " + outputDirStr_jf;
	system(cmd_mkdir_jf.c_str());
	string outputDirStr_Kmer = outputDirStr + "Kmer/";
	string cmd_mkdir_Kmer = "mkdir " + outputDirStr_Kmer;
	system(cmd_mkdir_Kmer.c_str());
	string outputDirStr_Kmer_sorted = outputDirStr + "Kmer_sorted/";
	string cmd_mkdir_Kmer_sorted = "mkdir " + outputDirStr_Kmer_sorted;
	system(cmd_mkdir_Kmer_sorted.c_str());			
	string outputDirStr_Kmer_enumerate = outputDirStr + "Kmer_enumerate/";
	string cmd_mkdir_Kmer_enumerate = "mkdir " + outputDirStr_Kmer_enumerate;
	system(cmd_mkdir_Kmer_enumerate.c_str());

	log_ofs << "start to initiate indexInfo" << endl;
	cout << "start to initiate indexInfo" << endl;
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;
	log_ofs << "end of initiating indexInfo" << endl;	

	cout << "start to generate genomeSubRegion fa" << endl;
	log_ofs << "start to generate genomeSubRegion fa" << endl;	
	int faReadSeqLength = 100;
	int jumpLength = faReadSeqLength - KmerLength;
	unsigned int genomeSubRegionSize = (indexInfo->returnIndexSize())/genomeSubRegionNum;
	for(int tmp = 0; tmp < genomeSubRegionNum; tmp++)
	{
		unsigned int genomeSubRegion_startPosInGenome = tmp * genomeSubRegionSize + 1;
		unsigned int genomeSubRegion_endPosInGenome = genomeSubRegion_startPosInGenome + genomeSubRegionSize - 1 + 1;
		log_ofs << "tmpRegion: " << tmp + 1 << endl;
		log_ofs << "posInWholeGenome: " << genomeSubRegion_startPosInGenome << " ~ " << genomeSubRegion_endPosInGenome << endl;
		unsigned int genomeSubRegion_startPosInChr, genomeSubRegion_startPosChrNameInt;
		indexInfo->getChrLocation(genomeSubRegion_startPosInGenome, &genomeSubRegion_startPosChrNameInt, &genomeSubRegion_startPosInChr);
		unsigned int genomeSubRegion_endPosInChr, genomeSubRegion_endPosChrNameInt;
		indexInfo->getChrLocation(genomeSubRegion_endPosInGenome, &genomeSubRegion_endPosChrNameInt, &genomeSubRegion_endPosInChr);		
		log_ofs << "posInChr: " << indexInfo->returnChrNameStr(genomeSubRegion_startPosChrNameInt) << ":" << genomeSubRegion_startPosInChr
			<< " ~ " << indexInfo->returnChrNameStr(genomeSubRegion_endPosChrNameInt) << ":" << genomeSubRegion_endPosInChr << endl;
		string tmpGenomeSubRegion_fa = outputDirStr_fa + int_to_str(tmp + 1) + ".fa";
		string tmpGenomeSubRegion_seq = indexInfo->returnChromStringSubstr(genomeSubRegion_startPosInGenome, genomeSubRegionSize);
		ofstream tmpGenomeSubRegion_fa_ofs(tmpGenomeSubRegion_fa.c_str());
		for(int tmpBase = 0; tmpBase < genomeSubRegionSize; tmpBase +=jumpLength)
		{
			int tmpSeqEndPos = tmpBase + faReadSeqLength - 1;
			if(tmpSeqEndPos >= genomeSubRegionSize)
				break;
			string tmpReadSeq = tmpGenomeSubRegion_seq.substr(tmpBase, faReadSeqLength);
			tmpGenomeSubRegion_fa_ofs << ">" << tmpBase + 1 << "-" << tmpSeqEndPos + 1 << endl
				<< tmpReadSeq << endl;
		}
		tmpGenomeSubRegion_fa_ofs.close();
	}
	cout << "start to generate genomeSubRegion jf" << endl;
	log_ofs << "start to generate genomeSubRegion jf" << endl;	
	for(int tmp = 0; tmp < genomeSubRegionNum; tmp++)
	{
		string tmpGenomeSubRegion_fa = outputDirStr_fa + int_to_str(tmp + 1) + ".fa";
		string tmpGenomeSubRegion_jf = outputDirStr_jf + int_to_str(tmp + 1) + ".jf";		
		string cmd_jf_count_subRegion = inputJellyFishBinPath + "/jellyfish count -o " + tmpGenomeSubRegion_jf 
			+ " -m " + KmerLengthStr + " -t " + thread_num_str + " -s 3000000000 -C " + tmpGenomeSubRegion_fa;
		system(cmd_jf_count_subRegion.c_str());
	}

	cout << "start to generate genomeSubRegion Kmer" << endl;
	log_ofs << "start to generate genomeSubRegion Kmer" << endl;
	for(int tmp = 0; tmp < genomeSubRegionNum; tmp++)
	{
		string tmpGenomeSubRegion_jf = outputDirStr_jf + int_to_str(tmp + 1) + ".jf";
		string tmpGenomeSubRegion_Kmer = outputDirStr_Kmer + int_to_str(tmp + 1) + ".Kmer";				
		string cmd_jf_dump_subRegion = inputJellyFishBinPath + "/jellyfish dump -t -c -o " + tmpGenomeSubRegion_Kmer 
			+ " " + tmpGenomeSubRegion_jf;
		system(cmd_jf_dump_subRegion.c_str());
	}

	cout << "start to generate genomeSubRegion sorted Kmer" << endl;
	log_ofs << "start to generate genomeSubRegion sorted Kmer" << endl;
	string tmpSortDir = outputDirStr + "/tmpSort/";
	string cmd_mkdir_tmpSortDir = "mkdir " + tmpSortDir;
	system(cmd_mkdir_tmpSortDir.c_str());
	for(int tmp = 0; tmp < genomeSubRegionNum; tmp++)
	{
		string tmpGenomeSubRegion_Kmer = outputDirStr_Kmer + int_to_str(tmp + 1) + ".Kmer";
		string tmpGenomeSubRegion_Kmer_sorted = outputDirStr_Kmer_sorted + int_to_str(tmp + 1) + ".Kmer.sorted";						
		string cmd_sort_subRegion = "sort -k1 -T " + tmpSortDir + " " + tmpGenomeSubRegion_Kmer + " > " 
			+ tmpGenomeSubRegion_Kmer_sorted;
		system(cmd_sort_subRegion.c_str());
	}

	if(!incompleteKmerCombinationOrNot_bool) // complete Kmer file, do not allow alien kmers
	{	
		cout << "start to enumerate all Kmers" << endl;
		log_ofs << "start to enumerate all Kmers" << endl;
		string Kmer_enumerate_1_file = outputDirStr_Kmer_enumerate + "1mer.txt";
		ofstream Kmer_enumerate_1_ofs(Kmer_enumerate_1_file.c_str());
		Kmer_enumerate_1_ofs << "A" << endl << "C" << endl << "G" << endl << "T" << endl;
		Kmer_enumerate_1_ofs.close();
		for(int tmp = 2; tmp <= KmerLength - 1; tmp++)
		{
			string last_Kmer_enumerate_file = outputDirStr_Kmer_enumerate + int_to_str(tmp-1) + "mer.txt"; 
			string tmp_Kmer_enumerate_file = outputDirStr_Kmer_enumerate + int_to_str(tmp) + "mer.txt";
			ifstream last_Kmer_enumerate_ifs(last_Kmer_enumerate_file.c_str());
			ofstream tmp_Kmer_enumerate_ofs(tmp_Kmer_enumerate_file.c_str());
			while(!last_Kmer_enumerate_ifs.eof())
			{
				string tmpStr;
				getline(last_Kmer_enumerate_ifs, tmpStr);
				if(tmpStr == "")
					break;
				tmp_Kmer_enumerate_ofs << tmpStr << "A" << endl;
				tmp_Kmer_enumerate_ofs << tmpStr << "C" << endl;
				tmp_Kmer_enumerate_ofs << tmpStr << "G" << endl;
				tmp_Kmer_enumerate_ofs << tmpStr << "T" << endl;
			}
			last_Kmer_enumerate_ifs.close();
			tmp_Kmer_enumerate_ofs.close();
		}
		// tmp = KmerLength
		string penultimate_Kmer_enumerate_file = outputDirStr_Kmer_enumerate + int_to_str(KmerLength-1) + "mer.txt";
		string final_Kmer_enumerate_file = outputDirStr_Kmer_enumerate + int_to_str(KmerLength) + "mer.txt";
		ifstream penultimate_Kmer_enumerate_ifs(penultimate_Kmer_enumerate_file.c_str());
		ofstream final_Kmer_enumerate_ofs(final_Kmer_enumerate_file.c_str());		
		while(!penultimate_Kmer_enumerate_ifs.eof())
		{
			string tmpStr;
			getline(penultimate_Kmer_enumerate_ifs, tmpStr);
			if(tmpStr == "")
				break;			
			string tmpStr_A = tmpStr + "A";
			string tmpStr_C = tmpStr + "C";
			string tmpStr_G = tmpStr + "G";
			string tmpStr_T = tmpStr + "T";
			bool tmpStr_A_canonical_bool = KmerCanonicalOrNot(tmpStr_A);
			bool tmpStr_C_canonical_bool = KmerCanonicalOrNot(tmpStr_C);
			bool tmpStr_G_canonical_bool = KmerCanonicalOrNot(tmpStr_G);
			bool tmpStr_T_canonical_bool = KmerCanonicalOrNot(tmpStr_T);
			if(tmpStr_A_canonical_bool)
				final_Kmer_enumerate_ofs << tmpStr_A << endl;
			if(tmpStr_C_canonical_bool)
				final_Kmer_enumerate_ofs << tmpStr_C << endl;
			if(tmpStr_G_canonical_bool)
				final_Kmer_enumerate_ofs << tmpStr_G << endl;
			if(tmpStr_T_canonical_bool)
				final_Kmer_enumerate_ofs << tmpStr_T << endl;											
		}
		penultimate_Kmer_enumerate_ifs.close();
		final_Kmer_enumerate_ofs.close();
	}
	else // allow alien kmers
	{
		string cmd_merge_jf = inputJellyFishBinPath + "/jellyfish merge -o " + outputDirStr_jf + "merged.jf";
		for(int tmp = 0; tmp < genomeSubRegionNum; tmp++)
		{	
			string tmpJfFile = " " + outputDirStr + "jf/" + int_to_str(tmp + 1) + ".jf";
			cmd_merge_jf += tmpJfFile;
		}
		log_ofs << "cmd_merge_jf:" << endl << cmd_merge_jf << endl;
		system(cmd_merge_jf.c_str());
		string cmd_merge_dumpJf2Kmer = inputJellyFishBinPath + "/jellyfish dump -t -c -o " + outputDirStr_Kmer 
			+ "merged.Kmer " + outputDirStr_jf + "merged.jf";
		system(cmd_merge_dumpJf2Kmer.c_str());
		log_ofs << "cmd_merge_dumpJf2Kmer:" << endl << cmd_merge_dumpJf2Kmer << endl;
		string cmd_merge_sortKmer = "sort -k1 -T " + tmpSortDir + " " + outputDirStr_Kmer + "merged.Kmer" + " > "
			+ outputDirStr_Kmer_sorted + "merged.Kmer.sorted";
		system(cmd_merge_sortKmer.c_str());
		log_ofs << "cmd_merge_sortKmer:" << endl << cmd_merge_sortKmer << endl;
	}


	cout << "start to generate subRegion specific kmers" << endl;
	cout << "start to initiate localRegionKmerIfsVec" << endl;
	vector<ifstream*> KmerReadableIfsVec;
	for(int tmp = 0; tmp < genomeSubRegionNum; tmp++)
	{
		string tmpSubRegionReadFile_kmer_sorted = outputDirStr_Kmer_sorted + int_to_str(tmp + 1) + ".Kmer.sorted";
		ifstream *tmpSubRegionKmer_ifs = new ifstream(tmpSubRegionReadFile_kmer_sorted.c_str());
		KmerReadableIfsVec.push_back(tmpSubRegionKmer_ifs);
	}
	
	vector<bool> kmer_file_end_bool_vec;
	for(int tmp = 0; tmp < genomeSubRegionNum; tmp ++)
		kmer_file_end_bool_vec.push_back(false);
	vector<string> currentKmerStrVec;
	for(int tmp = 0; tmp < genomeSubRegionNum; tmp ++)
		currentKmerStrVec.push_back("");
	//cout << "head lines" << endl;
	for(int tmp = 0; tmp < genomeSubRegionNum; tmp ++)
	{
		string tmpKmerFileStr;
		getline((*KmerReadableIfsVec[tmp]), tmpKmerFileStr);
		//cout << "tmpKmerFileStr: " << tmpKmerFileStr << endl;
		int tmpTabLoc = tmpKmerFileStr.find("\t");
		currentKmerStrVec[tmp] = tmpKmerFileStr.substr(0, tmpTabLoc);
		//cout << "tmpCurrentKmerStr: " << tmp << endl << currentKmerStrVec[tmp] << endl;
	}

	string outputFile = outputDirStr + "regionSpecific.Kmer";
	string mergedKmer_file;
	if(!incompleteKmerCombinationOrNot_bool)
		mergedKmer_file = outputDirStr_Kmer_enumerate + KmerLengthStr + "mer.txt";
	else
		mergedKmer_file = outputDirStr_Kmer_sorted + "merged.Kmer.sorted";;
	ifstream mergedKmer_ifs(mergedKmer_file.c_str());
	ofstream kmerClass_ofs(outputFile.c_str());
	while(!mergedKmer_ifs.eof())
	{
		string tmpStr;
		getline(mergedKmer_ifs, tmpStr);
		if(tmpStr == "")
			break;
		string tmpMergedStr;
		if(!incompleteKmerCombinationOrNot_bool)
			tmpMergedStr = tmpStr;
		else
		{
			int tabLoc = tmpStr.find("\t");
			tmpMergedStr = tmpStr.substr(0, tabLoc);
		}
		//int tabLoc = tmpStr.find("\t");
		//string tmpMergedStr = tmpStr.substr(0, tabLoc);
		int tmpKmerExistence_fileNum = 0;
		int tmpKmerExistence_lastFileId = -1;
		vector<bool> tmp_kmer_exist_in_sample_bool_vec;
		for(int tmp = 0; tmp < genomeSubRegionNum; tmp ++)
			tmp_kmer_exist_in_sample_bool_vec.push_back(false);

		//int tmpKmerExistingClassFlag = 0;
		for(int tmp = 0; tmp < genomeSubRegionNum; tmp ++)
		{
			if(tmpMergedStr == currentKmerStrVec[tmp])
			{
				tmp_kmer_exist_in_sample_bool_vec[tmp] = true;
				tmpKmerExistence_fileNum ++;
				tmpKmerExistence_lastFileId = tmp + 1;
				//tmpKmerExistingClassFlag += pow(2, tmp);
				if(!kmer_file_end_bool_vec[tmp])
				{
					if((*KmerReadableIfsVec[tmp]).eof())
						kmer_file_end_bool_vec[tmp] = true;
					else
					{	
						string tmpKmerFileStr;
						getline((*KmerReadableIfsVec[tmp]), tmpKmerFileStr);
						if(tmpKmerFileStr == "")
							kmer_file_end_bool_vec[tmp] = true;
						else
						{
							int tmpTabLoc = tmpKmerFileStr.find("\t");
							currentKmerStrVec[tmp] = tmpKmerFileStr.substr(0, tmpTabLoc);
						}
					}
				}
			}
			else
				tmp_kmer_exist_in_sample_bool_vec[tmp] = false;
		}
		if(tmpKmerExistence_fileNum == 0) // alien Kmer
			kmerClass_ofs << tmpMergedStr << "\t0" << endl;
		else if(tmpKmerExistence_fileNum == 1) // local region specific kmers
			kmerClass_ofs << tmpMergedStr << "\t" << tmpKmerExistence_lastFileId << endl;
		else if(tmpKmerExistence_fileNum > 1) // multi local region kmers
			kmerClass_ofs << tmpMergedStr << "\t" << inputClassificationGroupNum - 1 << endl;
		else
		{
			cout << "error ! tmpKmerExistence_fileNum < 0 !" << endl;
			exit(1); 
		}
	}

	for(int tmp = 0; tmp < genomeSubRegionNum; tmp++)
	{
		(*KmerReadableIfsVec[tmp]).close();
		delete KmerReadableIfsVec[tmp];
	}
	mergedKmer_ifs.close();
	kmerClass_ofs.close();	

	parameter_ifs.close();
	chrom_bit_file_ifs.close();
	log_ofs.close();
	delete indexInfo;
	free(chrom);
	return 0;
}