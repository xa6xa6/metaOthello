// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#include <iostream>
#include <fstream>
#include <map>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include "../Ye_implementation/othello.h"
#include "../Ye_implementation/muloth.h"
#include <chrono>
#include <inttypes.h>
#include "../Ye_implementation/io_helper.h"

using namespace std;

time_t nowtime;
struct tm *local;

typedef unsigned long long keyT;
typedef uint64_t valueT;

IOHelper<keyT,valueT> *helper;

keyT convertKmer2Key(char*s, int kmerLength)
{
    int tmpLength = 0;
    switch (*s) 
    {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            keyT ret = 0;
            while (*s == 'A' || *s == 'C' || *s =='T' || *s =='G')
            {
                ret <<=2;
                switch (*s) 
                {
                    case 'T':
                        ret ++;
                    case 'G':
                        ret ++;
                    case 'C':
                        ret ++;
                }
                s++;
                tmpLength ++;
                if(tmpLength == kmerLength)
                    return ret;
            }
    }
}

int main(int argc, char* argv[])
{
    if(argc != 5)
    {
        cout << "Executable indexFile querySeq kmer_length splitBit" << endl;
        exit(1);
    }
    string kmer_length_str = argv[3];
    string splitbit_str = argv[4];
    int kmer_length = atoi(kmer_length_str.c_str());
    int splitbit = atoi(splitbit_str.c_str());
    helper = new ConstantLengthKmerHelper<keyT,valueT>(KMERLENGTH,splitbit);

    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "start to load index" << endl;
    MulOth<keyT> * moth;
    moth = new MulOth<keyT>(argv[1], helper);

    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "end of loading index" << endl;
    cout << endl << "[" << asctime(local) << "start to do query" << endl;
    if(argc == 5) // query sequence provided in command line
    {
        keyT k = convertKmer2Key(argv[2], kmer_length);
        cout << "K: " << k << endl;
        // if(!helper->convertKmer2Key(argv[2], &k))
        // {
        //     cout << "error in helper->convert(argv[2], &k)" << endl;
        //     exit(1);
        // }
        valueT qv = moth->query(k);
        cout << "query seq: " << argv[2] << endl;
        //cout << "query key: " << k << endl;
        cout << "query results: " << qv << endl;
    }
    else
    {
        cout << "argc != 4" << endl;
        cout << "error! currently it is not supported " << endl;
    }
    nowtime = time(NULL);
    local = localtime(&nowtime);    
    cout << endl << "[" << asctime(local) << "end of doing query" << endl;
    delete moth;
    return 0;  
}