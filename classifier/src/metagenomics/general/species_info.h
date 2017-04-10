#ifndef SPECIES_INFO_H
#define SPECIES_IFNO_H
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

class Species_Info
{
private:
	vector< pair<int,string> > speciesIdNamePairVec; // speciesIdNamePairVec.size() == chromosomeSeqInfoIndexVecVec.size();
	vector< vector<int> >  speciesCorrespondingChromosomeSeqInfoIndexVecVec;

public:
	Species_Info()
	{}

	void initiate(ChromosomeSeq_Info_Vec& tmpChromosomeSeqInfoVec)
	{
		int chromosomeSeqInfoVecSize = tmpChromosomeSeqInfoVec.return_chromosomeSeqInfoVecSize();
		for(int tmp = 0; tmp < chromosomeSeqInfoVecSize; tmp++)
		{
			int tmpChromosomeSeqInfo_speciesId = tmpChromosomeSeqInfoVec.return_chromosomeSeqInfo_species_id(tmp);
			string tmpChromosomeSeqInfo_speciesName = tmpChromosomeSeqInfoVec.return_chromosomeSeqInfo_species_name(tmp);
			int currentSpeciesIdNamePairVecSize = speciesIdNamePairVec.size();
			bool alreadyExistsInCurrentSpeciesIdNamePairVec_bool = false;
			for(int tmp2 = 0; tmp2 < currentSpeciesIdNamePairVecSize; tmp2++)
			{
				int currentTmpSpeciesId = speciesIdNamePairVec[tmp2].first;
				string currentTmpSpeciesName = speciesIdNamePairVec[tmp2].second;
				if((currentTmpSpeciesId == tmpChromosomeSeqInfo_speciesId)
					&&(currentTmpSpeciesName == tmpChromosomeSeqInfo_speciesName)) // existing species
				{
					speciesCorrespondingChromosomeSeqInfoIndexVecVec[tmp2].push_back(tmp);
					alreadyExistsInCurrentSpeciesIdNamePairVec_bool = true;
					break;
				}
				else if(((currentTmpSpeciesId == tmpChromosomeSeqInfo_speciesId)
					&&(currentTmpSpeciesName != tmpChromosomeSeqInfo_speciesName))
					||((currentTmpSpeciesId != tmpChromosomeSeqInfo_speciesId)
					&&(currentTmpSpeciesName == tmpChromosomeSeqInfo_speciesName)))
				{
					cout << "(((currentTmpSpeciesId == tmpChromosomeSeqInfo_speciesId)"<< endl;
					cout << "&&(currentTmpSpeciesName != tmpChromosomeSeqInfo_speciesName))" << endl;
					cout << "||((currentTmpSpeciesId != tmpChromosomeSeqInfo_speciesId)" << endl;
					cout << "&&(currentTmpSpeciesName == tmpChromosomeSeqInfo_speciesName)))" << endl;
					exit(1);
				}
				else
				{}
			}
			if(!alreadyExistsInCurrentSpeciesIdNamePairVec_bool) // new species
			{
				speciesIdNamePairVec.push_back(pair<int,string>(
					tmpChromosomeSeqInfo_speciesId, tmpChromosomeSeqInfo_speciesName));
				vector<int> tmpSpeciesCorrespondingChromosomeSeqInfoIndexVec;
				tmpSpeciesCorrespondingChromosomeSeqInfoIndexVec.push_back(tmp);
				speciesCorrespondingChromosomeSeqInfoIndexVecVec.push_back(
					tmpSpeciesCorrespondingChromosomeSeqInfoIndexVec);
			}
		}
	}

	void generateTaxoId2sourceChromosomeSeqIdVecPairVec(vector< pair<int, vector<int> > >& taxoId2sourceChromosomeSeqIdVecPairVec)
	{
		for(int tmp = 0; tmp < speciesIdNamePairVec.size(); tmp++)
			taxoId2sourceChromosomeSeqIdVecPairVec.push_back(pair<int, vector<int> >(speciesIdNamePairVec[tmp].first,
				speciesCorrespondingChromosomeSeqInfoIndexVecVec[tmp]));
	}

	void printSpeciesInfo(string& speciesInfo_output_file)
	{
		ofstream speciesInfo_ofs(speciesInfo_output_file.c_str());
		for(int tmp = 0; tmp < speciesIdNamePairVec.size(); tmp++)
		{
			speciesInfo_ofs << speciesIdNamePairVec[tmp].first << "\t" << speciesIdNamePairVec[tmp].second << "\t";
			for(int tmp2 = 0; tmp2 < speciesCorrespondingChromosomeSeqInfoIndexVecVec[tmp].size(); tmp2++)
				speciesInfo_ofs << (speciesCorrespondingChromosomeSeqInfoIndexVecVec[tmp])[tmp2] << ",";
			speciesInfo_ofs << endl;
		}
		speciesInfo_ofs.close();
	}
};
#endif