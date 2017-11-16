// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#include <cstdio>
#include "io_helper.h"
#include <vector>
#include <algorithm>
using namespace std;
int main(int argc, char * argv[]) {
    BinaryKmerReader< uint64_t > breader(argv[1]);
    uint64_t k;
    while (breader.getNext(&k)) {
        printf("%llx\n", k);
    }

   return 0; 

}
