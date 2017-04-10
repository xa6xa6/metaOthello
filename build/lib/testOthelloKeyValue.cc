#include <iostream>
#include <map>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include "othello.h"
#include "muloth.h"
#include <chrono>
#include <inttypes.h>
#include <algorithm>
#include "io_helper.h"
#include "time.h"
using namespace std;
typedef unsigned long long keyT;
typedef uint16_t valueT;


vector<keyT> keys;
vector<valueT> values;

IOHelper<keyT,valueT> *helper;


int main(int argc, char * argv[]) {
    if (argc < 4) {
        printf(" splitbits keyInputFile OutputFile 'B'.... \n");
        printf(" splitbits: a number <=16, divide the keys into 2^(splitbits) sets, according to the Least (splitbit) significant bits. \n");
        printf(" splitbit <0: -splitbit keyInputfile, outputfile: skip build, just query\n");
        printf(" if B is specified, user int64+uint16 Kmer-value binary files");
        return 0;
    }
    bool usebinaryfile = false;
    if (argc >=5) {
        usebinaryfile = (argv[4][0]=='B');
        if (usebinaryfile) printf("Using binary file\n");
    }
    int splitbit;
    sscanf(argv[1],"%d",&splitbit);
    if (splitbit > 0) 
    helper = new ConstantLengthKmerHelper<keyT,valueT>(KMERLENGTH,splitbit);
    else
    helper = new ConstantLengthKmerHelper<keyT,valueT>(KMERLENGTH,-splitbit);

    MulOth<keyT,valueT> * moth;
    if (splitbit >=0) {
        printf("Split %d groups\n",1U<< splitbit);
        if (usebinaryfile) 
            moth = new MulOth<keyT>(VALUELENGTH,splitbit, new compressFileReader<keyT, uint16_t>(argv[2], helper, 8, 2, true));
        else 
            moth = new MulOth<keyT>(VALUELENGTH,splitbit, new KmerFileReader<keyT,uint16_t>(argv[2],helper,true) );
        if (!moth->buildsucc) return 1;

        printf("Build Succ, write to file %s\n", argv[3]);
        //moth.printall();
        moth->writeToFile(argv[3]);



        delete moth;
    }
    moth = new MulOth<keyT>( argv[3],helper);

        printf("Testing using keys from file %s\n", argv[2]);
        FileReader<keyT, valueT> *reader;
        if (usebinaryfile)
            reader = new compressFileReader<keyT,valueT>(argv[2], helper, 8, 2, true);
        else 
            reader = new KmerFileReader<keyT,valueT>(argv[2], helper, true); 
        uint64_t cnt = 0;
        while (true) {
            char buf[1024];
            keyT k;
            valueT v;
            if (!reader->getNext(&k,&v)) break;
            valueT qv = moth->query(k);
            if ((qv & ((1<<VALUELENGTH)-1))!=(v & ((1<<VALUELENGTH)-1) )) {
                printf("Err %llx -> %x : %x", k,
                       v & ((1<<VALUELENGTH)-1),
                       qv & ((1<<VALUELENGTH)-1));
                if (find(moth->removedKeys.begin(),moth->removedKeys.end(), k)!=moth->removedKeys.end()) {
                    printf(" But it is removed during construction\n");
                    continue;
                } else {
                    continue;
                    reader -> finish();
                    return 0;
                }
            }
            //    printf("%" PRIx64"--> %d\n", k,v);

        }
        reader -> finish();
        printf("Test Succ\n");
    delete moth;
    return 0;
}
