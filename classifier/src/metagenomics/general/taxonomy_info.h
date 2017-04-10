#ifndef TAXONOMY_INFO_H
#define TAXONOMY_INFO_H
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

using namespace std;

class Taxonomy_Info
{
private:
	vector< pair<int, string> > speciesId2namePairVec;
	vector< pair<int, string> > genusId2namePairVec;
	vector< pair<int, string> > familyId2namePairVec;
	vector< pair<int, string> > orderId2namePairVec;
	vector< pair<int, string> > classId2namePairVec;
	vector< pair<int, string> > phylumId2namePairVec;
public:
	Taxonomy_Info()
	{}

	void initiate()
	{}

	string return_name_from_id_species(int tmpId)
	{
		for(int tmp = 0; tmp < speciesId2namePairVec.size(); tmp++)
		{
			if((speciesId2namePairVec[tmp]).first == tmpId)
				return (speciesId2namePairVec[tmp]).second;
		}
		cout << "error! invalid species_id" << endl;
		exit(1);
	}

	string return_name_from_id_genus(int tmpId)
	{
		for(int tmp = 0; tmp < genusId2namePairVec.size(); tmp++)
		{
			if((genusId2namePairVec[tmp]).first == tmpId)
				return (genusId2namePairVec[tmp]).second;
		}
		cout << "error! invalid genus_id" << endl;
		exit(1);
	}

	string return_name_from_id_family(int tmpId)
	{
		for(int tmp = 0; tmp < familyId2namePairVec.size(); tmp++)
		{
			if((familyId2namePairVec[tmp]).first == tmpId)
				return (familyId2namePairVec[tmp]).second;
		}
		cout << "error! invalid family_id" << endl;
		exit(1);
	}

	string return_name_from_id_order(int tmpId)
	{
		for(int tmp = 0; tmp < orderId2namePairVec.size(); tmp++)
		{
			if((orderId2namePairVec[tmp]).first == tmpId)
				return (orderId2namePairVec[tmp]).second;
		}
		cout << "error! invalid order_id" << endl;
		exit(1);
	}

	string return_name_from_id_class(int tmpId)
	{
		for(int tmp = 0; tmp < classId2namePairVec.size(); tmp++)
		{
			if((classId2namePairVec[tmp]).first == tmpId)
				return (classId2namePairVec[tmp]).second;
		}
		cout << "error! invalid class_id" << endl;
		exit(1);
	}

	string return_name_from_id_phylum(int tmpId)
	{
		for(int tmp = 0; tmp < phylumId2namePairVec.size(); tmp++)
		{
			if((phylumId2namePairVec[tmp]).first == tmpId)
				return (phylumId2namePairVec[tmp]).second;
		}
		cout << "error! invalid phylum_id" << endl;
		exit(1);
	}				
};
#endif