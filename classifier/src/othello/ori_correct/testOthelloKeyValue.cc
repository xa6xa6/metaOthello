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
typedef uint64_t valueT;


vector<keyT> keys;
vector<valueT> values;

IOHelper<keyT,valueT> *helper;


int main(int argc, char * argv[]) {
    if (argc < 4) {
        printf(" splitbits keyInputFile OutputFile .... \n");
        printf(" splitbits: a number <=16, divide the keys into 2^(splitbits) sets, according to the Least (splitbit) significant bits. \n");
        printf(" splitbit == -1: -1, keyInputfile, outputfile: skip build, just query\n");
        return 0;
    }
    int splitbit;
    sscanf(argv[1],"%d",&splitbit);
    helper = new ConstantLengthKmerHelper<keyT,valueT>(KMERLENGTH,splitbit);

    MulOth<keyT> * moth;
    if (splitbit >=0) {
        printf("Split %d groups\n",1U<< splitbit);
        moth = new MulOth<keyT>(VALUELENGTH,argv[2],  splitbit, helper, true);
        if (!moth->buildsucc) return 1;

        printf("Build Succ, write to file %s\n", argv[3]);
        //moth.printall();
        moth->writeToFile(argv[3]);



        delete moth;
    }
    moth = new MulOth<keyT>( argv[3],helper);

    for (int i = 2; i< argc; i+=2) {
        printf("Testing using keys from file %s\n", argv[i]);
        FILE *pFile;
        pFile = fopen (argv[i],"r");
        uint64_t cnt = 0;
        while (true) {
            char buf[1024];
            if (fgets(buf,1024,pFile)==NULL) break;
            keyT k;
            valueT v;
            if (!helper->convert(buf, &k, &v)) break;
            valueT qv = moth->query(k);
            if ((qv & ((1<<VALUELENGTH)-1))!=(v & ((1<<VALUELENGTH)-1) )) {
                printf("Err %llx -> %x : %x", k,
                       v & ((1<<VALUELENGTH)-1),
                       qv & ((1<<VALUELENGTH)-1));
                if (find(moth->removedKeys.begin(),moth->removedKeys.end(), k)!=moth->removedKeys.end()) {
                    printf(" But it is removed during construction\n");
                    continue;
                } else {
                    fclose(pFile);
                    return 0;
                }
            }
            //    printf("%" PRIx64"--> %d\n", k,v);

        }
        fclose(pFile);
        printf("Test Succ\n");
    }
    delete moth;
    return 0;
}
