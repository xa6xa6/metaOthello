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
#include "../mps3Lib/read_block_test.h"
#include "../mps3Lib/otherFunc.h"
#include "../mps3Lib/index_info.h"
#include "general/queryConstantDef.h"
#include "general/querySeq_info.h"
using namespace std;

time_t nowtime;
struct tm *local;

typedef unsigned long long keyT;

#if VALUELENGTH==1 || VALUELENGTH == 2 || VALUELENGTH ==4 || VALUELENGTH==8
typedef uint8_t valueT;
#elif VALUELENGTH==16
typedef uint16_t valueT;
#elif VALUELENGTH==32
typedef uint32_t valueT;
#endif

void parseStr2vec_comma(vector<int>& intVec, string& tmpStr)
{
    int startLoc = 0;
    for(int tmp = 0; ; tmp++)
    {
        int tabLoc = tmpStr.find(",", startLoc);
        if(tabLoc == string::npos)
            break;
        string tmpField = tmpStr.substr(startLoc, tabLoc - startLoc);
        startLoc = tabLoc + 1;
        int tmpInt = atoi(tmpField.c_str());
        intVec.push_back(tmpInt);
    }
}

void rcmCharArray(char* oriChar, char* targetChar, int length)
{
    for(int tmp = 0; tmp < length; tmp++)
        *(targetChar + tmp) = getCharRevComp(*(oriChar + length - 1 - tmp));
    targetChar[length] = '\0';
}

void printCharArray(char* toPrintCharArray, int len)
{
    char tmp[len + 1];
    strncpy(tmp, toPrintCharArray, len);
    tmp[len] = '\0';
    printf("%s", tmp);
}

bool cmp_lexicoOrder(char* tmpChar_1, char* tmpChar_2, int len)
{
    char tmp2_1[len + 1];
    char tmp2_2[len + 1];
    strncpy(tmp2_1, tmpChar_1, len);
    strncpy(tmp2_2, tmpChar_2, len);
    tmp2_1[len] = '\0';
    tmp2_2[len] = '\0';
    return (strcmp(tmp2_1, tmp2_2) <= 0);
}

void queryKmer2KeyT(keyT &k, char* oriChar, char* rcmChar, int kmer_length, int tmpBaseIndex)
{
    //int length_divided_by2 = length / 2;
    //bool length_even_or_odd_bool = (length_divided_by2 * 2 == length);
    bool lexiOrderedKmer_ori_or_rcm_bool = (strcmp(oriChar, rcmChar) <= 0);
    // if(length_even_or_odd_bool)
    // {
    //     for(int tmp = 0; tmp < length_divided_by2; tmp++)
    //     {

    //     }
    // }
    // else
    // {
    //     for(int tmp = 0; tmp < length_divided_by2 + 1; tmp++)
    //     {
            
    //     }
    // }
    if(lexiOrderedKmer_ori_or_rcm_bool)
        lineToKVpair<keyT,valueT>(oriChar + tmpBaseIndex, &k, kmer_length);
}

int main(int argc, char** argv )
{
    if(argc != 6)
    {
        cout << "Executable indexFile PATTERN kmer_length inputQuerySeq jump_length" << endl;
        cout << "or" << endl;
        cout << "Executable indexFile PATTERN kmer_length inputQuerySeqFile outputResultsFile" << endl;
        exit(1);
    }
    string PATTERN_str = argv[2];
    bool cmd_or_file_mode_bool;
    if(PATTERN_str == "CMD")
        cmd_or_file_mode_bool = true;
    else if(PATTERN_str == "FILE")
        cmd_or_file_mode_bool = false;
    else
    {
        cout << "invalid PATTERN !" << endl;
        cout << "set as CMD or FILE" << endl;
        exit(1);
    }
    string kmer_length_str = argv[3];
    int kmer_length = atoi(kmer_length_str.c_str());

    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "start to load index" << endl;
    MulOth<VALUELENGTH, keyT> * moth = new MulOth<VALUELENGTH, keyT>(argv[1]);
    
    nowtime = time(NULL);
    local = localtime(&nowtime);        
    cout << endl << "[" << asctime(local) << "end of loading index" << endl;
    cout << endl << "[" << asctime(local) << "start to do query" << endl;

    // uint8_t fullUint8 = 0xFF;
    // cout << "uint8_t full: " << (int)fullUint8 << endl;
    // fullUint8 >>=1;
    // cout << "uint8_t left move: " << (int)fullUint8 << endl;

    if(cmd_or_file_mode_bool) // MODE 1: query sequence provided in command line
    {
        string inputQuerySeq = argv[4];
        string jump_length_str = argv[5];
        int jump_length = atoi(jump_length_str.c_str());
        string tmpId = "tmpQueryID";
        
        QuerySeq_Info tmpQuerySeqInfo;
        tmpQuerySeqInfo.initiate(tmpId, inputQuerySeq);
        
        vector< pair<valueT,int> > tmpQuerySeq_classIdCountPairVec;
        tmpQuerySeqInfo.querySeq_returnClassIdCountPairVec(tmpQuerySeq_classIdCountPairVec, moth, kmer_length, jump_length);
        
        valueT tmpQuerySeq_bestClassId;
        int tmpQuerySeq_bestClassKmerCount;
        tmpQuerySeqInfo.select_best_class(tmpQuerySeq_bestClassId, tmpQuerySeq_bestClassKmerCount, tmpQuerySeq_classIdCountPairVec);
        
        valueT tmpQuerySeq_secondBestClassId;
        int tmpQuerySeq_secondBestClassKmerCount;
        tmpQuerySeqInfo.select_secondBest_class(tmpQuerySeq_secondBestClassId, tmpQuerySeq_secondBestClassKmerCount, tmpQuerySeq_bestClassId, tmpQuerySeq_classIdCountPairVec);

        // string baseIndexStr = argv[5];
        // cout << "baseIndexStr: " << baseIndexStr << endl;
        
        // vector<int> baseIndexVec;
        // parseStr2vec_comma(baseIndexVec, baseIndexStr);
        // cout << "baseIndexVec.size(): " << baseIndexVec.size() << endl;
        // for(int tmp = 0; tmp < baseIndexVec.size(); tmp++)
        //     cout << "tmpBase: " << baseIndexVec[tmp] << endl;
        
        // char charArray_ori[1024];
        // strcpy(charArray_ori, argv[4]);
        // int seqLength = strlen(charArray_ori);
        // cout << "seqLength: " << seqLength << endl;
        // char charArray_rcm[1024];// = new char[seqLength];
        // rcmCharArray(charArray_ori, charArray_rcm, seqLength);
        // printf("tmpSeq_ori: %s\n", charArray_ori);
        // printf("tmpSeq_rcm: %s\n", charArray_rcm);

        // for(int tmp = 0; tmp < baseIndexVec.size(); tmp++)
        // {
        //     int tmpBaseIndex = baseIndexVec[tmp];
        //     printf("\ntmpBaseIndex: %d\n", tmpBaseIndex);
        //     bool cmp_ori_or_rcm_bool = cmp_lexicoOrder(charArray_ori + tmpBaseIndex, charArray_rcm + seqLength - (tmpBaseIndex + kmer_length), kmer_length);
        //     printf("ori_char_array:\n");
        //     printCharArray(charArray_ori + tmpBaseIndex, kmer_length);
        //     printf("\nrcm_char_array:\n");
        //     printCharArray(charArray_rcm + seqLength - (tmpBaseIndex + kmer_length), kmer_length);
        //     cout << "\ncmp_ori_or_rcm_bool: " << cmp_ori_or_rcm_bool << endl;
        //     keyT tmpK;
        //     valueT tmpQv;
        //     if(cmp_ori_or_rcm_bool)
        //         lineToKVpair<keyT,valueT>(charArray_ori + tmpBaseIndex, &tmpK, kmer_length);
        //     else
        //         lineToKVpair<keyT,valueT>(charArray_rcm + seqLength - (tmpBaseIndex + kmer_length), &tmpK, kmer_length);
        //     tmpQv = moth->query(tmpK);
        //     printf("query results: %d\n", tmpQv);
        // }
    }
    else // MODE 2: query sequence provided in file
    {
        // cout << "MODE 2 starts ..." << endl;
        // cout << "start to query sequence provided in file" << endl;
        // vector<int> baseIndexVec;
        // baseIndexVec.push_back(0);
        // baseIndexVec.push_back(20);
        // baseIndexVec.push_back(40);
        // baseIndexVec.push_back(60);
        // baseIndexVec.push_back(80);

        // string inputQuerySeqFile = argv[4];
        // string outputResultsFile = argv[5];
        // //FILE *pQuerySeqFile = fopen(inputQuerySeqFile.c_str(), "r");
        // FILE *pResultsFile = fopen(outputResultsFile.c_str(), "w");
        // ifstream querySeq_ifs(inputQuerySeqFile.c_str());
        // while (!querySeq_ifs.eof())
        // {
        //     string tmpId, tmpSeq_ori;
        //     getline(querySeq_ifs, tmpId);
        //     getline(querySeq_ifs, tmpSeq_ori);
        //     if((tmpId == "")||(tmpSeq_ori == ""))
        //         break;
        //     int seqLength = tmpSeq_ori.length();
        //     //int seqLength = strlen(charArray_ori);
        //     cout << "seqLength: " << seqLength << endl;
        //     char charArray_ori[1024];
        //     std::copy(tmpSeq_ori.begin(), tmpSeq_ori.end(), charArray_ori);
            
        //     char charArray_rcm[1024];// = new char[seqLength];
        //     rcmCharArray(charArray_ori, charArray_rcm, seqLength);
        //     //printf("tmpSeq_id: %s\n", charArray_id);
        //     printf("tmpSeq_ori: %s\n", charArray_ori);
        //     printf("tmpSeq_rcm: %s\n", charArray_rcm);
        //     fprintf(pResultsFile, "%s\n%s\n", charArray_ori, charArray_rcm);
        //     vector<valueT> tmpSeqQvVec;
        //     for(int tmp = 0; tmp < baseIndexVec.size(); tmp++)
        //     {
        //         int tmpBaseIndex = baseIndexVec[tmp];
        //         cout << "tmpBaseIndex: " << tmpBaseIndex << endl;
        //         bool cmp_ori_or_rcm_bool = cmp_lexicoOrder(charArray_ori + tmpBaseIndex, 
        //             charArray_rcm + seqLength - (tmpBaseIndex + kmer_length), kmer_length);

        //         keyT tmpK;
        //         if(cmp_ori_or_rcm_bool)
        //         {    
        //             if(!(lineToKVpair<keyT,valueT>(charArray_ori + tmpBaseIndex, &tmpK, kmer_length)))
        //                 break;
        //         }
        //         else
        //         {
        //             if(!(lineToKVpair<keyT,valueT>(charArray_rcm + seqLength - (tmpBaseIndex + kmer_length), &tmpK, kmer_length)))
        //                 break;
        //         }
        //         valueT tmpQv = moth->query(tmpK);
        //         tmpSeqQvVec.push_back(tmpQv);
        //         fprintf(pResultsFile, "%d,", tmpQv);
        //     }
        //     fprintf(pResultsFile, "\n");
        // }
        // querySeq_ifs.close();
        // fclose(pResultsFile); 
    }
    nowtime = time(NULL);
    local = localtime(&nowtime);    
    cout << endl << "[" << asctime(local) << "end of doing query" << endl;
    delete moth;
    return 0;  
}