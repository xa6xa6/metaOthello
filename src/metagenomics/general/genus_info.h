#ifndef GENUS_INFO_H
#define GENUS_IFNO_H
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

class Genus_Info
{
private:
	vector< pair<int,string> > genusIdNamePairVec; // genusIdNamePairVec.size() == chromosomeSeqInfoIndexVecVec.size();
	vector< vector<int> >  genusCorrespondingChromosomeSeqInfoIndexVecVec;

public:
	Genus_Info()
	{}

	void initiate(ChromosomeSeq_Info_Vec& tmpChromosomeSeqInfoVec)
	{
		int chromosomeSeqInfoVecSize = tmpChromosomeSeqInfoVec.return_chromosomeSeqInfoVecSize();
		for(int tmp = 0; tmp < chromosomeSeqInfoVecSize; tmp++)
		{
			int tmpChromosomeSeqInfo_genusId = tmpChromosomeSeqInfoVec.return_chromosomeSeqInfo_genus_id(tmp);
			string tmpChromosomeSeqInfo_genusName = tmpChromosomeSeqInfoVec.return_chromosomeSeqInfo_genus_name(tmp);
			int currentGenusIdNamePairVecSize = genusIdNamePairVec.size();
			bool alreadyExistsInCurrentGenusIdNamePairVec_bool = false;
			for(int tmp2 = 0; tmp2 < currentGenusIdNamePairVecSize; tmp2++)
			{
				int currentTmpGenusId = genusIdNamePairVec[tmp2].first;
				string currentTmpGenusName = genusIdNamePairVec[tmp2].second;
				if((currentTmpGenusId == tmpChromosomeSeqInfo_genusId)
					&&(currentTmpGenusName == tmpChromosomeSeqInfo_genusName)) // existing genus
				{
					genusCorrespondingChromosomeSeqInfoIndexVecVec[tmp2].push_back(tmp);
					alreadyExistsInCurrentGenusIdNamePairVec_bool = true;
					break;
				}
				else if(((currentTmpGenusId == tmpChromosomeSeqInfo_genusId)
					&&(currentTmpGenusName != tmpChromosomeSeqInfo_genusName))
					||((currentTmpGenusId != tmpChromosomeSeqInfo_genusId)
					&&(currentTmpGenusName == tmpChromosomeSeqInfo_genusName)))
				{
					cout << "(((currentTmpGenusId == tmpChromosomeSeqInfo_genusId)"<< endl;
					cout << "&&(currentTmpGenusName != tmpChromosomeSeqInfo_genusName))" << endl;
					cout << "||((currentTmpGenusId != tmpChromosomeSeqInfo_genusId)" << endl;
					cout << "&&(currentTmpGenusName == tmpChromosomeSeqInfo_genusName)))" << endl;
					exit(1);
				}
				else
				{}
			}
			if(!alreadyExistsInCurrentGenusIdNamePairVec_bool) // new genus
			{
				genusIdNamePairVec.push_back(pair<int,string>(
					tmpChromosomeSeqInfo_genusId, tmpChromosomeSeqInfo_genusName));
				vector<int> tmpGenusCorrespondingChromosomeSeqInfoIndexVec;
				tmpGenusCorrespondingChromosomeSeqInfoIndexVec.push_back(tmp);
				genusCorrespondingChromosomeSeqInfoIndexVecVec.push_back(
					tmpGenusCorrespondingChromosomeSeqInfoIndexVec);
			}
		}
	}

	void generateTaxoId2sourceChromosomeSeqIdVecPairVec(vector< pair<int, vector<int> > >& taxoId2sourceChromosomeSeqIdVecPairVec)
	{
		for(int tmp = 0; tmp < genusIdNamePairVec.size(); tmp++)
			taxoId2sourceChromosomeSeqIdVecPairVec.push_back(pair<int, vector<int> >(genusIdNamePairVec[tmp].first,
				genusCorrespondingChromosomeSeqInfoIndexVecVec[tmp]));
	}

	void printGenusInfo(string& genusInfo_output_file)
	{
		ofstream genusInfo_ofs(genusInfo_output_file.c_str());
		for(int tmp = 0; tmp < genusIdNamePairVec.size(); tmp++)
		{
			genusInfo_ofs << genusIdNamePairVec[tmp].first << "\t" << genusIdNamePairVec[tmp].second << "\t";
			for(int tmp2 = 0; tmp2 < genusCorrespondingChromosomeSeqInfoIndexVecVec[tmp].size(); tmp2++)
				genusInfo_ofs << (genusCorrespondingChromosomeSeqInfoIndexVecVec[tmp])[tmp2] << ",";
			genusInfo_ofs << endl;
		}
		genusInfo_ofs.close();
	}
};
#endif