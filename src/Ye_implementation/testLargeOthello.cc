#include <iostream>
#include <map>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include "othello.h"
#include "muloth.h"
#include <chrono>
#include <inttypes.h>
#include "io_helper.h"
using namespace std;
typedef unsigned long long keyT;
typedef uint64_t valueT;


vector<keyT> keys;
vector<valueT> values;

IOHelper<keyT,valueT> *helper;


int main(int argc, char * argv[]) {
    int splitbit;
    MulOth<keyT> *moth;
    const int NN = 204857600;
    moth = new MulOth<keyT>( 6, NN);
    for (int i = 0 ; i < NN; i++) {
        uint64_t vv =  moth->query((((uint64_t) i)<<32)+i+1);
        vv ^= i;
        if (vv & 0x3FULL) {
            printf("%llx, %x %x\n", i,vv,vv);
        }
    }
    return 0;
}
