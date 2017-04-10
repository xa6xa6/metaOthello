#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "io_helper.h"
#include "muloth.h"
#include <chrono>
#include <inttypes.h>
#include <string>
using namespace std;
typedef uint64_t keyT;
typedef uint16_t valueT;


int main(int argc, char * argv[]) { 
    if (argc < 8) {
        printf("args: descriptiveFilename KmerFnamePrefix KmerFnameSuffix Kmer_length splitbits OutputFile workingDirectory\n");
        return 0;
    }    
    string Kmer_length_str = argv[4];
    string splitbitStr = argv[5];
    int Kmer_length = atoi(Kmer_length_str.c_str());
    int splitbit = atoi(splitbitStr.c_str());
    
    IOHelper<keyT,uint16_t> *helper;
    taxoTreeBuilder<uint64_t, uint16_t> builder(argv[1],argv[2],argv[3],argv[7],Kmer_length, splitbit,false);

    MulOth<keyT,uint16_t> * moth;
    moth = new MulOth<keyT,uint16_t>(VALUELENGTH,splitbit, &builder);

    if (!moth->buildsucc) return 1;
    printf("Build Succ, write to file %s\n", argv[6]);
    moth->writeToFile(argv[6]);
    delete moth;
    
}






