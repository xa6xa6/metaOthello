#include <cstdio>
#include "io_helper.h"
#include <vector>
#include <algorithm>
using namespace std;
int main(int argc, char * argv[]) {
    ConstantLengthKmerHelper<uint64_t, uint16_t> iohelper(20,0);
    KmerFileReader< uint64_t,uint16_t > freader(argv[1], &iohelper);
    KVpair<uint64_t, uint16_t> pp;
    vector<uint64_t> VKmer;
    while (freader.getNext(&pp.k, &pp.v)) {
        VKmer.push_back(pp.k);
    }
    printf("Sorting\n");
    sort(VKmer.begin(),VKmer.end());
    BinaryKmerWriter< uint64_t > fwriter(argv[2]);
    for (auto a: VKmer) {
        fwriter.write(&a);
    }
    fwriter.finish();
   return 0; 

}
