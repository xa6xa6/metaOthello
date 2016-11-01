#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <random>
#include <inttypes.h>
#include "othelloindex.h"
#include <array>
#include <set>
using namespace std;
typedef uint64_t keyType;
typedef uint16_t valueType;
#define VALUELEN (16)
std::random_device rd;
std::uniform_int_distribution<keyType> disKey;
std::uniform_int_distribution<valueType> disValue;
std::mt19937_64 *g;
int main() {
    const int NN = 596;
    vector<keyType> kV;
    vector<valueType> vV;
    g = new std::mt19937_64 (rd());
    for (int i = 0; i < NN; i++) {
        kV.push_back(disKey(*g));
        valueType v = disValue(*g);
        v &= ((1<<VALUELEN)-1);
        vV.push_back(v);
    }
    set<valueType> mset;
    OthelloIndex <keyType> othidx(&kV[0],NN);
    for (int i = 0 ; i < NN; i++) {
        valueType t;
        t = othidx.query(kV[i]);
        mset.insert(t);
    }
    printf("%d\n", mset.size());
    printf("Write to file\n");
    FILE *pF;
    pF = fopen("testtmp","wb");
    unsigned char buf0x20[0x20];
    memset(&(buf0x20[0]),0,sizeof(buf0x20));
    othidx.exportInfo(buf0x20);
    fwrite(buf0x20, sizeof(buf0x20),1,pF);
    othidx.writeDataToBinaryFile(pF);
    fclose(pF);

    printf("Load from file\n");
    FILE *pFr;
    pFr = fopen("testtmp","rb");
    fread(buf0x20, sizeof(buf0x20),1, pFr);
    OthelloIndex<keyType> *newoth;
    newoth = new OthelloIndex<keyType>((unsigned char *) buf0x20);
    newoth->loadDataFromBinaryFile(pFr);

    for (int i = 0; i< NN; i++) {
        if (othidx.query(kV[i])!=newoth->query(kV[i])) {
            printf("Err: %llx %x %x\n", kV[i], othidx.query(kV[i]), newoth->query(kV[i]));
            return 0;
        }
    }
    if (mset.size()==NN && (*mset.begin())==0 && (*mset.rbegin())==NN-1 ) 
        printf("Test Succ!\n");
    else 
        printf("fail\n");

}
