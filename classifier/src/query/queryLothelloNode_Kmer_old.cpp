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
typedef uint8_t valueT;

int main(int argc, char** argv )
{
    if((argc != 4)&&(argc != 5))
    {
        cout << "Executable indexFile querySeq kmer_length" << endl;
        cout << "or" << endl;
        cout << "Executable indexFile inputQuerySeqFile outputResultsFile kmer_length" << endl;
        exit(1);
    }
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "start to load index" << endl;
    MulOth<VALUELENGTH, keyT> * moth = new MulOth<VALUELENGTH, keyT>(argv[1]);
    
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "end of loading index" << endl;
    cout << endl << "[" << asctime(local) << "start to do query" << endl;
    if(argc == 4) // MODE 1: query sequence provided in command line
    {
        string kmer_length_str = argv[3];
        int kmer_length = atoi(kmer_length_str.c_str());
        keyT k;
        valueT qv;        
        if(!lineToKVpair<keyT,valueT>(argv[2], &k, kmer_length))
            cout << "invalid querySeq: " << argv[2] << endl;
        else
        {
            qv = moth->query(k);
            printf("%d\n", qv);
        }
    }
    else // MODE 2: query sequence provided in file
    {
        string kmer_length_str = argv[4];
        int kmer_length = atoi(kmer_length_str.c_str());   
        string inputQuerySeqFile = argv[2];
        string outputResultsFile = argv[3];
        FILE *pQuerySeqFile = fopen(inputQuerySeqFile.c_str(), "r");
        //ifstream querySeq_ifs(inputQuerySeqFile.c_str());
        FILE *pResultsFile = fopen(outputResultsFile.c_str(), "w");
        while (1)//!querySeq_ifs.eof())
        {
            // string tmpStr;
            // getline(querySeq_ifs, tmpStr);
            // if(tmpStr == "")
            //     break;
            // string tmpKmerStr = tmpStr.substr(0,20);
            char buf[50];
            if(fgets(buf, 50, pQuerySeqFile)==NULL)
                break;
            keyT k;
            if(!lineToKVpair<keyT,valueT>(buf, &k, kmer_length)) 
                break;
            valueT qv = moth->query(k);
            //results_ofs << tmpKmerStr << "\t" << (int)qv << endl;
            fprintf(pResultsFile, "%s\t%d\n", buf, qv);
        }
        //results_ofs.close();
        fclose(pQuerySeqFile);
        fclose(pResultsFile); 
    }
    nowtime = time(NULL);
    local = localtime(&nowtime);    
    cout << endl << "[" << asctime(local) << "end of doing query" << endl;
    delete moth;
    return 0;  
}