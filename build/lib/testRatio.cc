#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <random>
#include <inttypes.h>
#include "othello.h"
#include <array>
#include <vector>
using namespace std;
typedef uint64_t keyType;
typedef uint16_t valueType;
#define VALUELEN (16)
std::random_device rd;
std::uniform_int_distribution<keyType> disKey;
std::uniform_int_distribution<valueType> disValue;
std::mt19937_64 *g;
int main() {
    const int NN = 96;
    vector<keyType> kV;
    vector<valueType> vV;
    g = new std::mt19937_64 (rd());
    for (int i = 0; i < NN; i++) {
        kV.push_back(disKey(*g));
        valueType v = disValue(*g);
        v &= ((1<<VALUELEN)-1);
        vV.push_back(v);
    }
    Othello<keyType> oth(VALUELEN, kV, vV,false);

    vector< uint32_t> cnt;
    vector<uint32_t> ans;
    vector< double> rat;

    cnt  = oth.getCnt();
    rat  = oth.getRatio();
    for (int i = 0; i < VALUELEN; i++) {  printf("%d ", cnt[i]);   }    printf("\n");
    for (int i = 0; i < VALUELEN; i++) {        printf("%d ", cnt[VALUELEN+i]);    }    printf("\n");
    for (int i = 0; i < VALUELEN; i++) {        printf("%.3lf ", rat[i]);   }    printf("\n");
    const int CCNT = 10000000;
    ans.resize(VALUELEN);
    for (int i = 0; i < CCNT; i++) {
        keyType T = disKey(*g);
        valueType V = oth.queryInt(T);
        for (int j = 0; j < VALUELEN; j++)
            ans[j] += ((V >> j) & 1);
    }
    printf("Query Rates:");
    for (int i = 0; i < VALUELEN; i++)
        printf("%.3lf ",ans[i]*1.0/CCNT);
    printf("\n");


    for (int i = 0; i < kV.size(); i++) {
        valueType v = oth.queryInt(kV[i]);
        if (v!=vV[i]) {            printf("Err!!!!");            return 0; }
    }
    printf("Verify Succ\n");


    printf("AlienPreference");
    oth.setAlienPreference(1.0);
    cnt  = oth.getCnt();
    rat  = oth.getRatio();
    for (int j = 0; j < VALUELEN; j++) ans[j] =0;
    for (int i = 0; i < VALUELEN; i++) {  printf("%d ", cnt[i]);   }    printf("\n");
    for (int i = 0; i < VALUELEN; i++) {        printf("%d ", cnt[VALUELEN+i]);    }    printf("\n");
    for (int i = 0; i < VALUELEN; i++) {        printf("%.3lf ", rat[i]);   }    printf("\n");
    for (int i = 0; i < CCNT; i++) {
        keyType T = disKey(*g);
        valueType V = oth.queryInt(T);
        for (int j = 0; j < VALUELEN; j++)
            ans[j] += ((V >> j) & 1);
    }
    printf("Query Rates:");
    for (int i = 0; i < VALUELEN; i++)
        printf("%.3lf ",ans[i]*1.0/CCNT);
    printf("\n");

    for (int i = 0; i < kV.size(); i++) {
        valueType v = oth.queryInt(kV[i]);
        if (v!=vV[i]) {            printf("Err!!!!");            return 0; }
    }
    printf("Verify Succ\n");
    return 0;
}
