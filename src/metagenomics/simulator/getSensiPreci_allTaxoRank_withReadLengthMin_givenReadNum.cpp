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

void parseStr2fieldVec_tab(vector<string>& tmpFieldVec, string& tmpStr)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find("\t", startLoc);
		if(tabLoc == string::npos)
			break;
		string tmpField = tmpStr.substr(startLoc, tabLoc-startLoc);
		tmpFieldVec.push_back(tmpField);
		startLoc = tabLoc + 1;
	}
	tmpFieldVec.push_back(tmpStr.substr(startLoc));
}

void parseStr2fieldVec_line(vector<string>& tmpFieldVec, string& tmpStr)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find("_", startLoc);
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
	if(argc != 5)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputReadAssignmentFile_withReadLength" << endl;
		cout << "#2 outputSensiPreciFile" << endl;
		cout << "#3 inputReadLengthMin" << endl;
		cout << "#4 inputTotalReadNumToComputeSensitivity" << endl;
		exit(1);
	}
	string inputReadAssignmentFile = argv[1];
	string outputSensiPreciFile = argv[2];
	//string inputSEreadFaFile = argv[3];
	string inputReadLengthMinStr = argv[3];
	int readLengthMin = atoi(inputReadLengthMinStr.c_str());// << endl;
	string inputTotalReadNumToComputeSensitivityStr = argv[4];
	int inputTotalReadNumToComputeSensitivity = atoi(inputTotalReadNumToComputeSensitivityStr.c_str());

	int total_num_tooShortRead = 0;
	int total_num_withUnidentifiedTaxoId = 0;
	int total_num = 0;

	int correct_num_species = 0;
	int incorrect_num_species = 0;
	int unmapped_num_species = 0;
	
	int correct_num_genus = 0;
	int incorrect_num_genus = 0;
	int unmapped_num_genus = 0;
	
	int correct_num_family = 0;
	int incorrect_num_family = 0;
	int unmapped_num_family = 0;
	
	int correct_num_order = 0;
	int incorrect_num_order = 0;
	int unmapped_num_order = 0;
	
	int correct_num_class = 0;
	int incorrect_num_class = 0;
	int unmapped_num_class = 0;
	
	int correct_num_phylum = 0;
	int incorrect_num_phylum = 0;
	int unmapped_num_phylum = 0;				

	//ifstream readFa_ifs(inputSEreadFaFile.c_str());
	ifstream readAssignment_ifs(inputReadAssignmentFile.c_str());
	unsigned long long tmpLineNO = 0;
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

        // string tmpReadNameStr, tmpReadSeqStr;
        // getline(readFa_ifs, tmpReadNameStr);
        // if(tmpReadNameStr == "")
        // 	break;
        // getline(readFa_ifs, tmpReadSeqStr);
        // if(tmpReadSeqStr == "")
        // 	break;
        vector<string> tmpFieldVec_tab;
        parseStr2fieldVec_tab(tmpFieldVec_tab, tmpStr);
        
        string tmpReadId_withReadLengthStr = tmpFieldVec_tab[0];
        vector<string> tmpFieldVec_comma_readId;
        parseStr2fieldVec_comma(tmpFieldVec_comma_readId, tmpReadId_withReadLengthStr);
        string tmpReadSeqLengthStr = tmpFieldVec_comma_readId[1];
        int tmpReadSeqLength = atoi(tmpReadSeqLengthStr.c_str());
        if(tmpReadSeqLength < readLengthMin)
        {
        	total_num_tooShortRead ++;
        	continue;
        }

        string tmpAssignmentStr_species = tmpFieldVec_tab[1];
        int tmpAssignmentId_species = atoi(tmpAssignmentStr_species.c_str());
        string tmpAssignmentStr_genus = tmpFieldVec_tab[2];
        int tmpAssignmentId_genus = atoi(tmpAssignmentStr_genus.c_str());
        string tmpAssignmentStr_family = tmpFieldVec_tab[3];
        int tmpAssignmentId_family = atoi(tmpAssignmentStr_family.c_str());                
        string tmpAssignmentStr_order = tmpFieldVec_tab[4];
        int tmpAssignmentId_order = atoi(tmpAssignmentStr_order.c_str());
        string tmpAssignmentStr_class = tmpFieldVec_tab[5];
        int tmpAssignmentId_class = atoi(tmpAssignmentStr_class.c_str());
        string tmpAssignmentStr_phylum = tmpFieldVec_tab[6];
        int tmpAssignmentId_phylum = atoi(tmpAssignmentStr_phylum.c_str());
        
        //cout << "tmpAssignmentStr_phylum: " << tmpAssignmentStr_phylum << endl;
        string tmpReadIdStr = tmpFieldVec_tab[0];
        string tmpReadIdStr_trimmed;
        if((tmpReadIdStr.at(0) == '>')||(tmpReadIdStr.at(0) == '@'))
        	tmpReadIdStr_trimmed = tmpReadIdStr.substr(1);
        else
        	tmpReadIdStr_trimmed = tmpReadIdStr;
        //cout << "tmpReadIdStr_trimmed: " << tmpReadIdStr_trimmed << endl;
        vector<string> tmpFieldVec_line;
        parseStr2fieldVec_line(tmpFieldVec_line, tmpReadIdStr_trimmed);

        string tmpTruthStr_species = tmpFieldVec_line[0];
        int tmpTruthId_species = atoi(tmpTruthStr_species.c_str());
        string tmpTruthStr_genus = tmpFieldVec_line[1];
        int tmpTruthId_genus = atoi(tmpTruthStr_genus.c_str());        
        string tmpTruthStr_family = tmpFieldVec_line[2];
        int tmpTruthId_family = atoi(tmpTruthStr_family.c_str());
        string tmpTruthStr_order = tmpFieldVec_line[3];
        int tmpTruthId_order = atoi(tmpTruthStr_order.c_str());
        string tmpTruthStr_class = tmpFieldVec_line[4];
        int tmpTruthId_class = atoi(tmpTruthStr_class.c_str());
        string tmpTruthStr_phylum = tmpFieldVec_line[5];
        int tmpTruthId_phylum = atoi(tmpTruthStr_phylum.c_str());
        //cout << "tmpTruthStr_phylum: " << tmpTruthStr_phylum << endl; 

        if((tmpTruthId_species < 0)||(tmpTruthId_genus < 0)||(tmpTruthId_family < 0)
        	||(tmpTruthId_order < 0)||(tmpTruthId_class < 0)||(tmpTruthId_phylum < 0))
       	{
        	total_num_withUnidentifiedTaxoId ++;
        	continue;
       	}
       	total_num ++;
        
        if(tmpAssignmentId_species < 0)
        	unmapped_num_species ++;
        else if(tmpAssignmentId_species == tmpTruthId_species)
        	correct_num_species ++;
        else
        	incorrect_num_species ++;

        if(tmpAssignmentId_genus < 0)
        	unmapped_num_genus ++;
        else if(tmpAssignmentId_genus == tmpTruthId_genus)
        	correct_num_genus ++;
        else
        	incorrect_num_genus ++;        

        if(tmpAssignmentId_family < 0)
        	unmapped_num_family ++;
        else if(tmpAssignmentId_family == tmpTruthId_family)
        	correct_num_family ++;
        else
        	incorrect_num_family ++;

        if(tmpAssignmentId_order < 0)
        	unmapped_num_order ++;
        else if(tmpAssignmentId_order == tmpTruthId_order)
        	correct_num_order ++;
        else
        	incorrect_num_order ++;   

        if(tmpAssignmentId_class < 0)
        	unmapped_num_class ++;
        else if(tmpAssignmentId_class == tmpTruthId_class)
        	correct_num_class ++;
        else
        	incorrect_num_class ++;

        if(tmpAssignmentId_phylum < 0)
        	unmapped_num_phylum ++;
        else if(tmpAssignmentId_phylum == tmpTruthId_phylum)
        	correct_num_phylum ++;
        else
        	incorrect_num_phylum ++;           
    }
	readAssignment_ifs.close();
	//readFa_ifs.close();

	int mapped_num_species = correct_num_species + incorrect_num_species;
	int mapped_num_genus = correct_num_genus + incorrect_num_genus;
	int mapped_num_family = correct_num_family + incorrect_num_family;
	int mapped_num_order = correct_num_order + incorrect_num_order;
	int mapped_num_class = correct_num_class + incorrect_num_class;
	int mapped_num_phylum = correct_num_phylum + incorrect_num_phylum;

	double sensi_species = ((double)correct_num_species/(double)inputTotalReadNumToComputeSensitivity) * 100;
	double preci_species = ((double)correct_num_species/(double)mapped_num_species) * 100;
	double sensi_genus = ((double)correct_num_genus/(double)inputTotalReadNumToComputeSensitivity) * 100;
	double preci_genus = ((double)correct_num_genus/(double)mapped_num_genus) * 100;
	double sensi_family = ((double)correct_num_family/(double)inputTotalReadNumToComputeSensitivity) * 100;
	double preci_family = ((double)correct_num_family/(double)mapped_num_family) * 100;
	double sensi_order = ((double)correct_num_order/(double)inputTotalReadNumToComputeSensitivity) * 100;
	double preci_order = ((double)correct_num_order/(double)mapped_num_order) * 100;
	double sensi_class = ((double)correct_num_class/(double)inputTotalReadNumToComputeSensitivity) * 100;
	double preci_class = ((double)correct_num_class/(double)mapped_num_class) * 100;
	double sensi_phylum = ((double)correct_num_phylum/(double)inputTotalReadNumToComputeSensitivity) * 100;
	double preci_phylum = ((double)correct_num_phylum/(double)mapped_num_phylum) * 100;

	double fscore_species = 2/(100/sensi_species + 100/preci_species);
	double fscore_genus = 2/(100/sensi_genus + 100/preci_genus);	
	double fscore_family = 2/(100/sensi_family + 100/preci_family);	
	double fscore_order = 2/(100/sensi_order + 100/preci_order);	
	double fscore_class = 2/(100/sensi_class + 100/preci_class);	
	double fscore_phylum = 2/(100/sensi_phylum + 100/preci_phylum);			

    ofstream sensiPreci_ofs(outputSensiPreciFile.c_str());
    sensiPreci_ofs << "Given total_num: " << inputTotalReadNumToComputeSensitivity << endl << endl << endl;
    sensiPreci_ofs << "In the result file: " << endl;
    sensiPreci_ofs << "total_num_withUnidentifiedTaxoId: " << total_num_withUnidentifiedTaxoId << endl;
    sensiPreci_ofs << "Total read #: " << total_num_tooShortRead + total_num << endl;
    sensiPreci_ofs << "Short read #: " << total_num_tooShortRead << endl;
    sensiPreci_ofs << "Short read perc: " 
    	<< (double)total_num_tooShortRead/(double)(total_num + total_num_tooShortRead) << endl;
    sensiPreci_ofs << "Long read #: " << total_num << endl;
    sensiPreci_ofs << "Long read perc: " 
    	<< (double)total_num/(double)(total_num + total_num_tooShortRead) << endl << endl;

    sensiPreci_ofs << endl;
	sensiPreci_ofs << "Taxo_Rank\tPrecision\tSensitivity\tF-score\tCorrect_Read_#\tIncorrect_Read_#\tUnmapped_Read_#" << endl;
	sensiPreci_ofs << "Phylum\t" << preci_phylum << "\t" << sensi_phylum << "\t" << fscore_phylum << "\t"
		<< correct_num_phylum << "\t" << incorrect_num_phylum << "\t" << unmapped_num_phylum << endl;
	sensiPreci_ofs << "Class\t" << preci_class << "\t" << sensi_class << "\t" << fscore_class << "\t"
		<< correct_num_class << "\t" << incorrect_num_class << "\t" << unmapped_num_class << endl;
	sensiPreci_ofs << "Order\t" << preci_order << "\t" << sensi_order << "\t" << fscore_order << "\t"
		<< correct_num_order << "\t" << incorrect_num_order << "\t" << unmapped_num_order << endl;				
	sensiPreci_ofs << "Family\t" << preci_family << "\t" << sensi_family << "\t" << fscore_family << "\t"
		<< correct_num_family << "\t" << incorrect_num_family << "\t" << unmapped_num_family << endl;
	sensiPreci_ofs << "Genus\t" << preci_genus << "\t" << sensi_genus << "\t" << fscore_genus << "\t"
		<< correct_num_genus << "\t" << incorrect_num_genus << "\t" << unmapped_num_genus << endl;
	sensiPreci_ofs << "Species\t" << preci_species << "\t" << sensi_species << "\t" << fscore_species << "\t"
		<< correct_num_species << "\t" << incorrect_num_species << "\t" << unmapped_num_species << endl;
    sensiPreci_ofs.close();
	return 0;
}