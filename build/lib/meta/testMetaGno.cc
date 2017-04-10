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
#include <queue>
using namespace std;
typedef unsigned long long keyT;
typedef uint64_t valueT;


vector<keyT> keys;
vector<valueT> values;

IOHelper<keyT,valueT> *helper;


int main(int argc, char * argv[]) {
    if (argc < 4) {
        printf(" splitbits InputfnameList OutputFile .... \n");
        printf(" splitbits: a number <=16, divide the keys into 2^(splitbits) sets, according to the Least (splitbit) significant bits. \n");
        printf(" splitbit == -1: -1, keyInputfile, outputfile: skip build, just query\n");
        return 0;
    }
    int splitbit;
    sscanf(argv[1],"%d",&splitbit);
    helper = new ConstantLengthKmerHelper<keyT,valueT>(KMERLENGTH,splitbit);
    FILE *fin; fin = fopen(argv[2],"r");
    char buf[1024];
    vector<string> fnames;
    int kfid = 0;
    while (true) {
        if(fgets(buf,sizeof(buf),fin)==NULL) break;
        if (strlen(buf)<1) break;
        printf("%s\n", buf);
        string ss(buf);
        fnames.push_back(ss);
        kfid ++;
    }
    vector< KmerFileReader<keyT,valueT> > readers;
    for (int i = 0; i < fnames.size(); i++) {
     KmerFileReader<keyT,valueT> reader(fnames[i].c_str(),helper);
     readers.push_back( reader );
    }
    vector<int32_t> tempV;
    priority_queue< pair<keyT, uint32_t>,vector< pair<keyT,uint32_t> > , greater< pair< keyT, uint32_t> > > PQ;
    for (int i = 0 ; i < fnames.size(); i++) {
        PQ.push(make_pair(0ULL,(uint32_t) i));
    }
    tempV.resize(fnames.size(),-1);
    keyT cur = 0;
    int cnt = 0;
    map<uint32_t, uint64_t> freq;
    vector<uint16_t> vSampleIDs;
    while (!PQ.empty()) {
        pair<keyT,uint32_t> PP = PQ.top(); PQ.pop();
        //printf("From file %d key %llx\n", PP.second, PP.first);
        keyT k; valueT v;
        if (PP.first != cur) {
                //new key discovered , push old key ....
                freq[cnt] ++;
                cur = PP.first;
                cnt = 1;
                vSampleIDs.clear();
        }
        else 
            if (tempV[PP.second]>=0) cnt ++;
        vSampleIDs.push_back(PP.second);
        if (readers[PP.second].getNext(&k,&v)) {
            tempV[PP.second] = v;
            PQ.push(make_pair(k,PP.second));
        }
        
    }
    for (auto a = freq.begin(); a!=freq.end(); a++)
        printf("%lld %lld\n",a->first,a->second);

}
