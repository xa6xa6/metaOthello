#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <inttypes.h>
const char * ac = "ATGC";
int main(int argc, char * argv[]) {
    int n,l,mx;
    sscanf(argv[1],"%d",&n);
    sscanf(argv[2],"%d",&l);
    sscanf(argv[3],"%d",&mx);
    for (int i = 0 ;i < n; i++) {
        uint64_t LL = rand();
        LL <<= 28;
        LL ^= rand();
        LL<<=20;
        LL |= (uint64_t) i;
        for ( int j = 0 ; j < l; j++) {
            printf("%c",  ac[LL%4]);
            LL/=4;
        }
        printf(" %d\n",rand()%mx);
    }

}
