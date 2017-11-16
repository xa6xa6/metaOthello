// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#include <iostream>
#include <map>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include "Ye_implementation/othello.h"
#include "Ye_implementation/muloth.h"
#include <chrono>
#include <inttypes.h>
#include "Ye_implementation/io_helper.h"
using namespace std;
typedef unsigned long long keyT;
typedef uint64_t valueT;

vector<keyT> keys;
vector<valueT> values;


void add(char * fname) {
	printf("Reading file %s \n", fname);
	FILE *pFile;
	pFile = fopen (fname,"r");
	uint64_t cnt = 0;
	while (true) {
		char buf[1024];
		if (fgets(buf,1024,pFile)==NULL) break;
		keyT k;
		valueT v;
		if (!lineToKVpair(buf, &k, &v)) break;
		keys.push_back(k);
		values.push_back(v);
		cnt++;
	}
	fclose(pFile);
	cout << "Keycnt" << human(cnt) <<"  totKey# "<< human(keys.size())<<endl;
}


int main(int argc, char * argv[]) {
	if (argc < 4) {
		printf(" splitbits keyInputFile OutputFile .... \n");
		printf(" splitbits: a number <=16, divide the keys into 2^(splitbits) sets, according to the Least (splitbit) significant bits. \n");
		printf(" splitbit == -1: -1, keyInputfile, outputfile: skip build, just query\n");
		return 0;
	}
	int splitbit;
	sscanf(argv[1],"%d",&splitbit);
	MulOth<VALUELENGTH,keyT> * moth;
	if (splitbit >=0) {
		printf("Split %d groups\n",1U<< splitbit);
		moth = new MulOth<VALUELENGTH, keyT>(argv[2],  splitbit, true) ;
        if (!moth->buildsucc) return 1;
          
		printf("Build Succ, write to file %s\n", argv[3]);
		//moth.printall();
		moth->writeToFile(argv[3]);
		delete moth;
	}
	moth = new MulOth<VALUELENGTH, keyT>( argv[3]);
	printf("start to do querying\n");
	for (int i = 2; i< argc; i+=2) 
	{
		printf("Testing using keys from file %s\n", argv[i]);
		FILE *pFile;
		pFile = fopen (argv[i],"r");
		uint64_t cnt = 0;
		while (true)
		{
			char buf[1024];
			if(fgets(buf,1024,pFile)==NULL) 
			{
				printf("(fgets(buf,1024,pFile)==NULL)\n");
				break;
			}
			keyT k;
			valueT v;
			if(!lineToKVpair<keyT,valueT>(buf, &k, &v))
			{
				printf("(!lineToKVpair<keyT,valueT>(buf, &k, &v))\n");
				break;
			}
			valueT qv = moth->query(k);
			printf("%s\t%d\n", buf, qv);
			if ((qv & ((1<<VALUELENGTH)-1))!=(v & ((1<<VALUELENGTH)-1) )) 
			{
                printf("Err %llx -> %x : %x\n", k,v,qv);
				fclose(pFile);
				return 0;
			}
		}
		fclose(pFile);
		printf("Test Succ\n");
	}
	return 0;
}
