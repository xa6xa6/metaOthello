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
        printf(" splitbits keyFile MulOthFile 'B/N' kmerlength valuelength \n");
    }
    bool usebinaryfile = false;
    if (argc >=5) {
        usebinaryfile = (argv[4][0]=='B');
        if (usebinaryfile) printf("Using binary file\n");
    }
    int splitbit;
    sscanf(argv[1],"%d",&splitbit);
    int kmerlength;
    sscanf(argv[5],"%d",&kmerlength);
    int VALUELENGTH;
    sscanf(argv[6],"%d",&VALUELENGTH);
    helper = new ConstantLengthKmerHelper<keyT,valueT>(kmerlength,splitbit);
    MulOth<keyT> * moth = new MulOth<keyT>( argv[3],helper);

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
                printf("%llx %d\n", k,qv);

        }
        reader -> finish();
    delete moth;
    return 0;
}
