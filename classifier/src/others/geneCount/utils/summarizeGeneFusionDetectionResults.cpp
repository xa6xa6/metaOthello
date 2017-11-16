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
#include "general/fusionDetectionEvaluation_info.h"

//typedef map<string, int> GeneId2indexMap;

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 geneIdListFile" << endl;
		cout << "#2 fusionReadDetectionResults" << endl;
		cout << "#3 fusionDetectionSummaryFile_prefix" << endl;
		exit(1);
	}
	int pairReadLength = 200;

	string geneIdListFile = argv[1];
	string fusionReadDetectionResults = argv[2];
	string fusionDetectionSummaryFile_prefix = argv[3];
	string fusionDetectionSummaryFile_discriminative 
		= fusionDetectionSummaryFile_prefix + ".discriminative.txt";
	string fusionDetectionSummaryFile_compatible
		= fusionDetectionSummaryFile_prefix + ".compatble.txt";	

	cout << "start to initiate geneId2indexMap ..." << endl;
	FusionDetectionEvaluation_Vec_Info tmpFusionEvaluationVecInfo;
	tmpFusionEvaluationVecInfo.initiate_geneIdList(geneIdListFile);
	cout << "start to do load_fusionReadInfoFile ..." << endl;
	tmpFusionEvaluationVecInfo.load_fusionReadInfoFile(fusionReadDetectionResults, pairReadLength);
	cout << "start to do print_fusion_summary ..." << endl;
	tmpFusionEvaluationVecInfo.print_fusion_summary(fusionDetectionSummaryFile_discriminative,
		fusionDetectionSummaryFile_compatible);
	cout << "start to load fusionReadDetection results" << endl;

	return 0;
}