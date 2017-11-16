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
#include "../Ye_implementation/othello.h"
#include "../Ye_implementation/muloth.h"
#include <chrono>
#include <inttypes.h>
#include "../Ye_implementation/io_helper.h"
#include "../mps3Lib/read_block_test.h"
#include "../mps3Lib/otherFunc.h"
#include "../mps3Lib/index_info.h"
#include "../query/general/queryConstantDef.h"
#include "../query/general/querySeq_info.h"
#include "general/regionAssignment_info.h"

using namespace std;

time_t nowtime;
struct tm *local;

bool repetitiveKmerWindowOrNot(int tmpDiscrimitiveKmerCount, int tmpRepetitiveKmerCount, double repetitiveKmerProportionMin)
{
	double tmpRepetitiveKmerProportion = (double)tmpRepetitiveKmerCount/(double)(tmpDiscrimitiveKmerCount + tmpRepetitiveKmerCount);
	if(tmpRepetitiveKmerProportion >= repetitiveKmerProportionMin)
		return true;
	else
		return false;
}

bool overlapWithLastRepetitiveKmerWindowOrNot(string& lastRepetitiveKmerWindow_chrName,
	int lastRepetitiveKmerWindow_startPos, int lastRepetitiveKmerWindow_endPos, 
	string& tmpChrName, int tmpStartPos, int tmpEndPos, int& mergedRegion_startPos, int& mergedRegion_endPos)
{
	if(lastRepetitiveKmerWindow_chrName == tmpChrName)
	{
		if(lastRepetitiveKmerWindow_endPos < tmpStartPos)
			return false;
		else
		{
			mergedRegion_startPos = lastRepetitiveKmerWindow_startPos;
			mergedRegion_endPos = tmpEndPos;
			return true;
		}
	}
	else
		return false;
}

void parseWindowKmerStr(string& tmpStr, string& tmpChrName, int& tmpStartPos, 
	int& tmpEndPos, int& discrimitiveKmerCount, int& repetitiveKmerCount)
{
	int tabLoc = tmpStr.find("\t");
	string tmpChrNamePosStr = tmpStr.substr(0, tabLoc);
	int lineLoc_1 = tmpChrNamePosStr.find("_");
	int lineLoc_2 = tmpChrNamePosStr.find("_", lineLoc_1 + 1);
	string tmpKmerCountStr = tmpStr.substr(tabLoc + 1);
	int lineLoc_3 = tmpKmerCountStr.find("_");
	int lineLoc_4 = tmpKmerCountStr.find("_", lineLoc_3 + 1);
	tmpChrName = tmpChrNamePosStr.substr(0, lineLoc_1);
	string tmpStartPosStr = tmpChrNamePosStr.substr(lineLoc_1 + 1, lineLoc_2 - lineLoc_1 - 1);
	string tmpEndPosStr = tmpChrNamePosStr.substr(lineLoc_2 + 1);
	tmpStartPos = atoi(tmpStartPosStr.c_str());
	tmpEndPos = atoi(tmpEndPosStr.c_str());
	string discrimitiveKmerCountStr = tmpKmerCountStr.substr(lineLoc_3 + 1, lineLoc_4 - lineLoc_3 - 1);
	string repetitiveKmerCountStr = tmpKmerCountStr.substr(lineLoc_4 + 1);
	discrimitiveKmerCount = atoi(discrimitiveKmerCountStr.c_str());
	repetitiveKmerCount = atoi(repetitiveKmerCountStr.c_str());
}

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputGenomeIndex inputRepeatKmerWindowFile outputMergedRepeatRegionFile repetitiveKmerProportionMin" << endl;
		exit(1);
	}
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load indexInfo" << endl;
	string indexFolderPath = argv[1];
	cout << "initiate indexInfo ..." << endl;	
	indexFolderPath += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	indexInfo->initiate_withoutLoadingSeq();
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... end of loading indexInfo" << endl;

	string inputWindowKmerFile = argv[2];
	string outputMergedRepeatRegionFile = argv[3];
	string repetitiveKmerProportionMinStr = argv[4];
	double repetitiveKmerProportionMin = atof(repetitiveKmerProportionMinStr.c_str());
	ifstream windowKmer_ifs(inputWindowKmerFile.c_str());
	ofstream mergedRepeatRegion_ofs(outputMergedRepeatRegionFile.c_str());
	string lastRepetitiveKmerWindow_chrName = "NULL";
	int lastRepetitiveKmerWindow_startPos, lastRepetitiveKmerWindow_endPos;
	unsigned int repetitiveRegionSize = 0;
	while(!windowKmer_ifs.eof())
	{
		string tmpStr;
		getline(windowKmer_ifs, tmpStr);
		if(tmpStr == "")
		{
			// printout the last region
			if(lastRepetitiveKmerWindow_chrName != "NULL")
			{	
				mergedRepeatRegion_ofs << lastRepetitiveKmerWindow_chrName << "\t" << lastRepetitiveKmerWindow_startPos << "\t" 			
					<< lastRepetitiveKmerWindow_endPos << "\t" << lastRepetitiveKmerWindow_endPos - lastRepetitiveKmerWindow_startPos + 1 << endl;
				repetitiveRegionSize += (lastRepetitiveKmerWindow_endPos - lastRepetitiveKmerWindow_startPos + 1);
			}
			break;
		}
		string tmpChrName; 
		int tmpStartPos, tmpEndPos, tmpDiscrimitiveKmerCount, tmpRepetitiveKmerCount;
		parseWindowKmerStr(tmpStr, tmpChrName, tmpStartPos, tmpEndPos, tmpDiscrimitiveKmerCount, tmpRepetitiveKmerCount);
		bool repetitiveKmerWindow_bool = repetitiveKmerWindowOrNot(tmpDiscrimitiveKmerCount, tmpRepetitiveKmerCount, repetitiveKmerProportionMin);
		if(repetitiveKmerWindow_bool)
		{
			int mergedRegion_startPos, mergedRegion_endPos;
			bool overlapWithLastRepetitiveKmerWindow_bool = overlapWithLastRepetitiveKmerWindowOrNot(
				lastRepetitiveKmerWindow_chrName, lastRepetitiveKmerWindow_startPos, lastRepetitiveKmerWindow_endPos, 
				tmpChrName, tmpStartPos, tmpEndPos, mergedRegion_startPos, mergedRegion_endPos);
			if(overlapWithLastRepetitiveKmerWindow_bool)
			{
				// update with new chrPos
				lastRepetitiveKmerWindow_startPos = mergedRegion_startPos;
				lastRepetitiveKmerWindow_endPos = mergedRegion_endPos;				
			}
			else
			{
				// print out lastRepetitiveRegion
				if(lastRepetitiveKmerWindow_chrName != "NULL")
				{	
					mergedRepeatRegion_ofs << lastRepetitiveKmerWindow_chrName << "\t" << lastRepetitiveKmerWindow_startPos << "\t" 			
						<< lastRepetitiveKmerWindow_endPos << "\t" << lastRepetitiveKmerWindow_endPos - lastRepetitiveKmerWindow_startPos + 1 << endl;
					repetitiveRegionSize += (lastRepetitiveKmerWindow_endPos - lastRepetitiveKmerWindow_startPos + 1);
				}			
				// update with new chrName+pos
				lastRepetitiveKmerWindow_chrName = tmpChrName;
				lastRepetitiveKmerWindow_startPos = tmpStartPos;
				lastRepetitiveKmerWindow_endPos = tmpEndPos;		
			} 
		}
		else
		{
			// int mergedRegion_startPos, mergedRegion_endPos;
			// bool overlapWithLastRepetitiveKmerWindow_bool = overlapWithLastRepetitiveKmerWindowOrNot(
			// 	lastRepetitiveKmerWindow_chrName, lastRepetitiveKmerWindow_startPos, lastRepetitiveKmerWindow_endPos, 
			// 	tmpChrName, tmpStartPos, tmpEndPos, mergedRegion_startPos, mergedRegion_endPos);
			// if(overlapWithLastRepetitiveKmerWindow_bool)
			// {}
			// else // print out lastRepetitiveRegion
			// {
			// 	if(lastRepetitiveKmerWindow_chrName != "NULL")
			// 	{
			// 		mergedRepeatRegion_ofs << lastRepetitiveKmerWindow_chrName << "\t" << lastRepetitiveKmerWindow_startPos << "\t" 
			// 			<< lastRepetitiveKmerWindow_endPos << "\t" << lastRepetitiveKmerWindow_endPos - lastRepetitiveKmerWindow_startPos + 1 << endl;				
			// 		repetitiveRegionSize += (lastRepetitiveKmerWindow_endPos - lastRepetitiveKmerWindow_startPos + 1);
			// 	}
			// }
		}
	}
	cout << "repetitiveRegionSize: " << endl << repetitiveRegionSize << endl;
	mergedRepeatRegion_ofs.close();
	windowKmer_ifs.close();
	return 0;
}