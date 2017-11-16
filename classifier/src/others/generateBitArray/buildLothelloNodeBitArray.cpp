// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#include <iostream>
#include <map>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include "../Ye_implementation/othello.h"
#include "../Ye_implementation/muloth.h"
#include <chrono>
#include <inttypes.h>
#include "../Ye_implementation/io_helper.h"
using namespace std;
typedef unsigned long long keyT;
typedef uint64_t valueT;

vector<keyT> keys;
vector<valueT> values;

IOHelper<keyT,valueT> *helper;

int main(int argc, char * argv[]) {
    if (argc != 5)
    {
        printf("Executable Kmer_length splitbits keyInputFile OutputFile .... \n");
        return 0;
    }
    string Kmer_length_str = argv[1];
    string splitbitStr = argv[2];
    int Kmer_length = atoi(Kmer_length_str.c_str());
    int splitbit = atoi(splitbitStr.c_str());

    helper = new ConstantLengthKmerHelper<keyT,valueT>(Kmer_length,splitbit);
    MulOth<keyT> * moth;
    moth = new MulOth<keyT>(VALUELENGTH, argv[3], splitbit, helper, true);
    if (!moth->buildsucc) return 1;
    printf("Build Succ, write to file %s\n", argv[3]);
    moth->writeToFile(argv[4]);
    delete moth;

    // MulOth<VALUELENGTH,keyT> * moth;
    // printf("Split %d groups\n",1U<< splitbit);
    // moth = new MulOth<VALUELENGTH, keyT>(argv[2],  splitbit, true);   
    // printf("Build Succ, write to file %s\n", argv[3]);
    // //moth.printall();
    // moth->writeToFile(argv[3]);
    // delete moth;
    return 0;
}