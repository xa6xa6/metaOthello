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
#include "../../Ye_implementation/othello.h"
#include "../../Ye_implementation/muloth.h"
#include <chrono>
#include <inttypes.h>
#include "../../Ye_implementation/io_helper.h"
#include "../../mps3Lib/read_block_test.h"
#include "../../mps3Lib/otherFunc.h"
#include "../../mps3Lib/index_info.h"
#include "../../query/general/queryConstantDef.h"
#include "../../query/general/querySeq_info.h"
//#include "../../genomeRegionAssignment/general/regionAssignment_info.h"
//#include "../general/reissuedGenomeID2TaxoID_info.h"
#include "../general/NCBIfullTaxoID2Name_info.h"
#include "../general/bacterialTaxo_info.h"
#include "../general/reissuedGenomeID2TaxoID_info.h"
#include "../general/MSKW_info.h"
#include "../general/taxoClassAssignment_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#1 inputSpeciesTaxoIdFile" << endl;
		cout << "#2 NCBIfullTaxoId2NameFile" << endl; 
		cout << "#3 outputDir" << endl;
		exit(1);
	}
	string inputSpeciesTaxoIdFile = argv[1];
	string NCBIfullTaxoId2NameFile = argv[2];
	string outputDirStr = argv[3];

	outputDirStr += "/";
  	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());
   	string log_file = outputDirStr + "log.txt";
   	ofstream log_ofs(log_file.c_str());

	BacterialTaxo_Info bacterialTaxoInfo;
	bacterialTaxoInfo.initiate_bacterialTaxoFile_NCBIfullTaxoId2NameFile(
		inputSpeciesTaxoIdFile, NCBIfullTaxoId2NameFile);
	bacterialTaxoInfo.reissueTaxoIdName_all();
	bacterialTaxoInfo.initiate_lowerRankTaxoReissuedId_to_higherRankTaxoReissuedId_array();
	//string outputDir_taxoInfo = outputDirStr + "taxoInfo/";
	bacterialTaxoInfo.print(outputDirStr);

    int taxo_num_total = bacterialTaxoInfo.return_taxo_num_total();   	
    log_ofs << "taxo_num_total: " << taxo_num_total << endl;
    log_ofs.close();
	return 0;
}