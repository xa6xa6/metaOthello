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

void parseStr2fieldVec_comma(vector<string>& tmpFieldVec, string& tmpStr)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find(",", startLoc);
		if(tabLoc == string::npos)
			break;
		string tmpField = tmpStr.substr(startLoc, tabLoc-startLoc);
		tmpFieldVec.push_back(tmpField);
		startLoc = tabLoc + 1;
	}
	tmpFieldVec.push_back(tmpStr.substr(startLoc));
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputReadAssignmentFile" << endl;
		cout << "#2 outputStatsFile" << endl;
		exit(1);
	}
	string inputReadAssignmentFile = argv[1];
	string outputStatsFile = argv[2];

	int total_num = 0;
	int mapped_num = 0;
	int unmapped_num = 0;	

	ifstream readAssignment_ifs(inputReadAssignmentFile.c_str());
	ofstream stats_ofs(outputStatsFile.c_str());

	unsigned long long tmpLineNO = 0;
	string tmp1stLine;
	getline(readAssignment_ifs, tmp1stLine);
	while(!readAssignment_ifs.eof())
	{
		string tmpStr;
		getline(readAssignment_ifs, tmpStr);
		if(tmpStr == "")
			break;
        tmpLineNO ++;
        unsigned long long tmpThousandIndex = tmpLineNO / 5000000;
        if(tmpLineNO == tmpThousandIndex * 5000000)
            cout << "Processed Line #: " << tmpLineNO << endl;

        total_num ++;
        vector<string> tmpFieldVec_comma;
        parseStr2fieldVec_comma(tmpFieldVec_comma, tmpStr);
        
        string tmpAssignmentStr = tmpFieldVec_comma[2];
        //
		if(tmpAssignmentStr == "NA")
			unmapped_num ++;
		else
		{
			int tmpAssignmentId = atoi(tmpAssignmentStr.c_str());
			if(tmpAssignmentId < 0)
				unmapped_num ++;
			else
				mapped_num ++;
		}
	}

	double mappingRate = ((double)mapped_num/(double)total_num) * 100;

	stats_ofs << "Total     #:\t" << total_num << endl << endl;
	stats_ofs << "taxo_rank\tmapping_rate\tmapped_read_#\tunmapped_read_#" << endl;
	stats_ofs << "\t" << mappingRate << "%\t" << mapped_num << "\t" << unmapped_num << endl;

	stats_ofs.close();
	readAssignment_ifs.close();
	return 0;
}