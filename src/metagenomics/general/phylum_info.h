#ifndef PHYLUM_INFO_H
#define PHYLUM_IFNO_H
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
#include "../../mps3Lib/read_block_test.h"
#include "../../mps3Lib/otherFunc.h"
//#include "taxonomy_info.h"
//#include "chromosomeSeq_info_vec.h"
using namespace std;

class Phylum_Info
{
private:
	vector< pair<int,string> > phylumIdNamePairVec; // phylumIdNamePairVec.size() == chromosomeSeqInfoIndexVecVec.size();
	vector< vector<int> > phylumCorrespondingChromosomeSeqInfoIndexVecVec;

public:
	Phylum_Info()
	{}

	void initiate(ChromosomeSeq_Info_Vec& tmpChromosomeSeqInfoVec)
	{
		int chromosomeSeqInfoVecSize = tmpChromosomeSeqInfoVec.return_chromosomeSeqInfoVecSize();
		for(int tmp = 0; tmp < chromosomeSeqInfoVecSize; tmp++)
		{
			int tmpChromosomeSeqInfo_phylumId = tmpChromosomeSeqInfoVec.return_chromosomeSeqInfo_phylum_id(tmp);
			string tmpChromosomeSeqInfo_phylumName = tmpChromosomeSeqInfoVec.return_chromosomeSeqInfo_phylum_name(tmp);
			int currentPhylumIdNamePairVecSize = phylumIdNamePairVec.size();
			bool alreadyExistsInCurrentPhylumIdNamePairVec_bool = false;
			for(int tmp2 = 0; tmp2 < currentPhylumIdNamePairVecSize; tmp2++)
			{
				int currentTmpPhylumId = phylumIdNamePairVec[tmp2].first;
				string currentTmpPhylumName = phylumIdNamePairVec[tmp2].second;
				if((currentTmpPhylumId == tmpChromosomeSeqInfo_phylumId)
					&&(currentTmpPhylumName == tmpChromosomeSeqInfo_phylumName)) // existing phylum
				{
					phylumCorrespondingChromosomeSeqInfoIndexVecVec[tmp2].push_back(tmp);
					alreadyExistsInCurrentPhylumIdNamePairVec_bool = true;
					break;
				}
				else if(((currentTmpPhylumId == tmpChromosomeSeqInfo_phylumId)
					&&(currentTmpPhylumName != tmpChromosomeSeqInfo_phylumName))
					||((currentTmpPhylumId != tmpChromosomeSeqInfo_phylumId)
					&&(currentTmpPhylumName == tmpChromosomeSeqInfo_phylumName)))
				{
					cout << "(((currentTmpPhylumId == tmpChromosomeSeqInfo_phylumId)"<< endl;
					cout << "&&(currentTmpPhylumName != tmpChromosomeSeqInfo_phylumName))" << endl;
					cout << "||((currentTmpPhylumId != tmpChromosomeSeqInfo_phylumId)" << endl;
					cout << "&&(currentTmpPhylumName == tmpChromosomeSeqInfo_phylumName)))" << endl;
					exit(1);
				}
				else
				{}
			}
			if(!alreadyExistsInCurrentPhylumIdNamePairVec_bool) // new phylum
			{
				phylumIdNamePairVec.push_back(pair<int,string>(
					tmpChromosomeSeqInfo_phylumId, tmpChromosomeSeqInfo_phylumName));
				vector<int> tmpPhylumCorrespondingChromosomeSeqInfoIndexVec;
				tmpPhylumCorrespondingChromosomeSeqInfoIndexVec.push_back(tmp);
				phylumCorrespondingChromosomeSeqInfoIndexVecVec.push_back(
					tmpPhylumCorrespondingChromosomeSeqInfoIndexVec);
			}
		}
	}

	void generateTaxoId2sourceChromosomeSeqIdVecPairVec(vector< pair<int, vector<int> > >& taxoId2sourceChromosomeSeqIdVecPairVec)
	{
		for(int tmp = 0; tmp < speciesIdNamePairVec.size(); tmp++)
			taxoId2sourceChromosomeSeqIdVecPairVec.push_back(pair<int, vector<int> >(phylumIdNamePairVec[tmp].first,
				phylumCorrespondingChromosomeSeqInfoIndexVecVec[tmp]));
	}

	void printPhylumInfo(string& phylumInfo_output_file)
	{
		ofstream phylumInfo_ofs(phylumInfo_output_file.c_str());
		for(int tmp = 0; tmp < phylumIdNamePairVec.size(); tmp++)
		{
			phylumInfo_ofs << phylumIdNamePairVec[tmp].first << "\t" << phylumIdNamePairVec[tmp].second << "\t";
			for(int tmp2 = 0; tmp2 < phylumCorrespondingChromosomeSeqInfoIndexVecVec[tmp].size(); tmp2++)
				phylumInfo_ofs << (phylumCorrespondingChromosomeSeqInfoIndexVecVec[tmp])[tmp2] << ",";
			phylumInfo_ofs << endl;
		}
		phylumInfo_ofs.close();
	}
};
#endif