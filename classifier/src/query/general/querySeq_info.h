// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#ifndef QUERYSEQ_INFO_H
#define QUERYSEQ_INFO_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>
#include "queryConstantDef.h"

using namespace std;

const keyT char2keyArray[26] = {0x0, 0x0, 0x1, 0x0, 0x0, 0x0, 0x2, //ABCDEFG
                        0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, //HIJKLMN
                        0x0, 0x0, 0x0, 0x0, 0x0, 0x3,      //OPQRST
                        0x0, 0x0, 0x0, 0x0, 0x0, 0x3};//HIJKLMN

const keyT twoChar2keyArray[676] = 
    {0x0,0x0,0x1,0x0,0x0,0x0,0x2, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x3, 0x0,0x0,0x0,0x0,0x0,0x0,//AA,AB,AC,AD,AE,AF,AG,AH,AI,AJ,AK,AL,AM,AN,AO,AP,AQ,AR,AS,AT,AU,AV,AW,AX,AY,AZ
     0x4,0x0,0x5,0x0,0x0,0x0,0x6, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x7, 0x0,0x0,0x0,0x0,0x0,0x0,//BA,BB,BC,BD,BE,BF,BG,BH,BI,BJ,BK,BL,BM,BN,BO,BP,BQ,BR,BS,BT,BU,BV,BW,BX,BY,BZ
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//CA,CB,CC,CD,CE,CF,CG,CH,CI,CJ,CK,CL,CM,CN,CO,CP,CQ,CR,CS,CT,CU,CV,CW,CX,CY,CZ
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//D
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//E
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//F
     0x8,0x0,0x9,0x0,0x0,0x0,0xA, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0xB, 0x0,0x0,0x0,0x0,0x0,0x0,//G
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//H
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//I
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//J
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//K
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//L
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//M
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//N
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//O
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//P
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//Q
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//R
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//S
     0xC,0x0,0xD,0x0,0x0,0x0,0xE, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0xF, 0x0,0x0,0x0,0x0,0x0,0x0,//T
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//U
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//V
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//W
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//X
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,//Y
     0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0, 0x0,0x0,0x0,0x0,0x0,0x0};//Z

class QuerySeq_Info
{
private:
	string queryId;
	char charArray_ori[QUERY_SEQ_LENGTH_MAX];
	char charArray_rcm[QUERY_SEQ_LENGTH_MAX];
	int seqLength;
public:
	QuerySeq_Info()
	{}

	void initiate(string& tmpId, string& tmpSeq_ori)
	{
		queryId = tmpId;
		int tmpSeqLength = tmpSeq_ori.length();
		seqLength = tmpSeqLength;
		std::copy(tmpSeq_ori.begin(), tmpSeq_ori.end(), charArray_ori);
		this->rcmCharArray(charArray_ori, charArray_rcm, tmpSeqLength);
        //cout << "after initiating ..." << endl;
        //cout << "charArray_ori: " << charArray_ori << endl;
        //cout << "charArray_rcm: " << charArray_rcm << endl;
	}

	void rcmCharArray(char* oriChar, char* targetChar, int length)
	{
    	for(int tmp = 0; tmp < length; tmp++)
        	*(targetChar + tmp) = getCharRevComp(*(oriChar + length - 1 - tmp));
    	targetChar[length] = '\0';
	}

	bool return_kmerInSeq_canonical_ori_or_rcm(int tmpKmerStartBaseIndexInSeq, int kmer_length)
	{
		return this->cmp_lexicoOrder(charArray_ori + tmpKmerStartBaseIndexInSeq, 
    		charArray_rcm + seqLength - (tmpKmerStartBaseIndexInSeq + kmer_length), kmer_length);
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

    keyT returnKmerKey_baseIndex(int tmpKmerStartBaseIndexInSeq, int kmer_length)
    {
    	//bool cmp_ori_or_rcm_bool = this->return_kmerInSeq_canonical_ori_or_rcm(tmpKmerStartBaseIndexInSeq, kmer_length);
        //bool cmp_ori_or_rcm_bool = true;
        //keyT tmpK;
        //if(cmp_ori_or_rcm_bool)
		keyT tmpK_for, tmpK_rcm;
        tmpK_for = this->convertKmer2Key(charArray_ori + tmpKmerStartBaseIndexInSeq, kmer_length);
        tmpK_rcm = this->convertKmer2Key(charArray_rcm + seqLength - (tmpKmerStartBaseIndexInSeq + kmer_length), kmer_length);
        	//lineToKVpair<keyT,valueT>(charArray_ori + tmpKmerStartBaseIndexInSeq, &tmpK_for, kmer_length);                
        //else
            //lineToKVpair<keyT,valueT>(charArray_rcm + seqLength - (tmpKmerStartBaseIndexInSeq + kmer_length), &tmpK_rcm, kmer_length);
    	if(tmpK_for <= tmpK_rcm)  
            return tmpK_for;
       else
            return tmpK_rcm;
    }

    keyT returnKmerFullKey(int kmer_length)
    {
        switch(kmer_length)
        {
            case 16:
                return 0xFFFFFFFF;
            case 18:
                return 0xFFFFFFFFF;
            case 20:
                return 0xFFFFFFFFFF;
            case 22:
                return 0xFFFFFFFFFFF;
            case 24:
                return 0xFFFFFFFFFFFF;
            case 25:
                return 0x3FFFFFFFFFFFF;
            case 28:
                return 0xFFFFFFFFFFFFFF;
            case 29:
                return 0x3FFFFFFFFFFFFFF;
            case 30:
                return 0xFFFFFFFFFFFFFFF;
            case 31:
                return 0x3FFFFFFFFFFFFFFF;                                                    
            default:
                cout << "error! invalid kmer_length: " << kmer_length << endl;
                exit(1); 
        }
    }

    void generateKmerKeyVec_overlapped_oriOnly(vector<keyT>& tmpKmerKeyVec_oriOnly, 
        vector<int> tmpKmerStartBaseIndexVecInSeq, int kmer_length, int jump_length)
    {
        keyT tmpKmerFullKey = this->returnKmerFullKey(kmer_length);
        //keyT tmpKmerFullKey = 0xFFFFFFFFFF;
        keyT tmpK_ori_1st = this->convertKmer2Key(charArray_ori, kmer_length);
        // initiate 1st key
        tmpKmerKeyVec_oriOnly.push_back(tmpK_ori_1st);
        // others
        keyT tmpK_ori_last = tmpK_ori_1st;
        for(int tmp = 1; tmp < tmpKmerStartBaseIndexVecInSeq.size(); tmp++)
        {
            keyT tmpK_ori_current = tmpK_ori_last;
            tmpK_ori_current<<=(2*jump_length);
            keyT tmpK_ori_current_newAddedBase;
            if(jump_length > 1) 
                tmpK_ori_current_newAddedBase = this->convertKmer2Key(
                    charArray_ori + tmpKmerStartBaseIndexVecInSeq[tmp-1] + kmer_length, jump_length);
            else
                tmpK_ori_current_newAddedBase = char2keyArray[*(charArray_ori + tmpKmerStartBaseIndexVecInSeq[tmp-1] + kmer_length)-'A'];
            tmpK_ori_current = (tmpK_ori_current|tmpK_ori_current_newAddedBase)&tmpKmerFullKey;
            tmpKmerKeyVec_oriOnly.push_back(tmpK_ori_current);
            tmpK_ori_last = tmpK_ori_current;
        }
    }

    void generateKmerKeyVec_overlapped_rcmOnly(vector<keyT>& tmpKmerKeyVec_rcmOnly, 
        vector<int>& tmpKmerStartBaseIndexVecInSeq, int kmer_length, int jump_length)
    {
        keyT tmpKmerFullKey = this->returnKmerFullKey(kmer_length);
        // get tmpKmerStartBaseIndexVecInSeq_rcm
        vector<int> tmpKmerStartBaseIndexVecInSeq_rcm;
        for(int tmp = tmpKmerStartBaseIndexVecInSeq.size() - 1; tmp >= 0; tmp--)
        {
            int tmpKmerStartBaseIndexInSeq = tmpKmerStartBaseIndexVecInSeq[tmp];
            int tmpKmerEndBaseIndexInSeq = tmpKmerStartBaseIndexInSeq + kmer_length - 1;
            tmpKmerStartBaseIndexVecInSeq_rcm.push_back(seqLength - 1 - tmpKmerEndBaseIndexInSeq);
        }

        keyT tmpK_rcm_1st = this->convertKmer2Key(charArray_rcm + tmpKmerStartBaseIndexVecInSeq_rcm[0], kmer_length);
        // initiate 1st key
        //lineToKVpair<keyT,valueT>(charArray_rcm + tmpKmerStartBaseIndexVecInSeq_rcm[0], &tmpK_rcm_1st, kmer_length);
        tmpKmerKeyVec_rcmOnly.push_back(tmpK_rcm_1st);
        // others
        keyT tmpK_rcm_last = tmpK_rcm_1st;
        for(int tmp = 1; tmp < tmpKmerStartBaseIndexVecInSeq_rcm.size(); tmp++)
        {
            keyT tmpK_rcm_current = tmpK_rcm_last;
            tmpK_rcm_current<<=(2*jump_length);
            keyT tmpK_rcm_current_newAddedBase;
            if(jump_length > 1)
                tmpK_rcm_current_newAddedBase = this->convertKmer2Key(
                    charArray_rcm + tmpKmerStartBaseIndexVecInSeq_rcm[tmp-1] + kmer_length, jump_length);
            else
                tmpK_rcm_current_newAddedBase = char2keyArray[*(charArray_rcm + tmpKmerStartBaseIndexVecInSeq_rcm[tmp-1] + kmer_length)-'A'];
            tmpK_rcm_current = (tmpK_rcm_current|tmpK_rcm_current_newAddedBase)&tmpKmerFullKey;
            tmpKmerKeyVec_rcmOnly.push_back(tmpK_rcm_current);
            tmpK_rcm_last = tmpK_rcm_current;
        }
    }

    void generateKmerKeyVec_overlapped(vector<keyT>& tmpKmerKeyVec, int kmer_length, int jump_length)
    {
        //cout << "generateKmerKeyVec_overlapped starts ..." << endl;
        if(jump_length >= kmer_length)
        {
            cout << "error! (jump_length >= kmer_length) in generateKmerKeyVec_overlapped()" << endl;
            exit(1);
        }

        vector<int> tmpKmerStartBaseIndexVecInSeq;
        this->generateKmerStartBaseIndexVec(tmpKmerStartBaseIndexVecInSeq, kmer_length, jump_length);
        //cout << "tmpKmerStartBaseIndexInSeq.size(): " << tmpKmerStartBaseIndexVecInSeq.size() << endl;
        // Alter strategy 1:
        vector<keyT> tmpKmerKeyVec_ori;
        this->generateKmerKeyVec_overlapped_oriOnly(tmpKmerKeyVec_ori, tmpKmerStartBaseIndexVecInSeq, kmer_length, jump_length);
        vector<keyT> tmpKmerKeyVec_rcm;
        this->generateKmerKeyVec_overlapped_rcmOnly(tmpKmerKeyVec_rcm, tmpKmerStartBaseIndexVecInSeq, kmer_length, jump_length);

        int kmer_num = tmpKmerStartBaseIndexVecInSeq.size();
        for(int tmp = 0; tmp < kmer_num; tmp++)
        {
            //cout << "tmp: " << tmp << endl;
            //cout << "tmpKmerKeyVec_ori[tmp]: " << tmpKmerKeyVec_ori[tmp] << endl;
            //cout << "tmpKmerKeyVec_rcm[kmer_num - 1 - tmp]: " << tmpKmerKeyVec_rcm[kmer_num - 1 - tmp] << endl;
            if(tmpKmerKeyVec_ori[tmp] < tmpKmerKeyVec_rcm[kmer_num - 1 - tmp])
                tmpKmerKeyVec.push_back(tmpKmerKeyVec_ori[tmp]);
            else
                tmpKmerKeyVec.push_back(tmpKmerKeyVec_rcm[kmer_num - 1 - tmp]);
        }        

        // Alter strategy 2:
        // for(int tmp = 0; tmp < tmpKmerStartBaseIndexVecInSeq.size(); tmp++)
        // {
        //     int tmpKmerStartBaseIndexInSeq = tmpKmerStartBaseIndexVecInSeq[tmp];
        //     keyT tmpKey = this->returnKmerKey_baseIndex(tmpKmerStartBaseIndexInSeq, kmer_length);
        //     tmpKmerKeyVec.push_back(tmpKey);
        // }
    }

    void generateKmerKeyVec_nonOverlapped(vector<keyT>& tmpKmerKeyVec, int kmer_length, int jump_length)
    {
        if(jump_length < kmer_length)
        {
            cout << "error! (jump_length < kmer_length) in generateKmerKeyVec_nonOverlapped()" << endl;
            exit(1);
        }
        vector<int> tmpKmerStartBaseIndexVecInSeq;
        this->generateKmerStartBaseIndexVec(tmpKmerStartBaseIndexVecInSeq, kmer_length, jump_length);
        for(int tmp = 0; tmp < tmpKmerStartBaseIndexVecInSeq.size(); tmp++)
        {
            int tmpKmerStartBaseIndexInSeq = tmpKmerStartBaseIndexVecInSeq[tmp];
            keyT tmpKey = this->returnKmerKey_baseIndex(tmpKmerStartBaseIndexInSeq, kmer_length);
            tmpKmerKeyVec.push_back(tmpKey);
        }
    }

    valueT queryKmer_baseIndex(MulOth<keyT>* moth, int tmpKmerStartBaseIndexInSeq, int kmer_length)
    {
    	keyT tmpKmerKey = this->returnKmerKey_baseIndex(tmpKmerStartBaseIndexInSeq, kmer_length);
    	//keyT tmpKmerKey = 1000;
        return moth->query(tmpKmerKey);
    }

    void queryKmerVec_baseIndexVec(MulOth<keyT>* moth, 
    	vector<valueT>& tmpValueVec, vector<int>& tmpKmerStartBaseIndexVecInSeq, int kmer_length)
    {
    	for(int tmp = 0; tmp < tmpKmerStartBaseIndexVecInSeq.size(); tmp++)
    	{
    		valueT tmpValue = this->queryKmer_baseIndex(moth, tmpKmerStartBaseIndexVecInSeq[tmp], kmer_length);
    		tmpValueVec.push_back(tmpValue);
    	}
    }

    void generateKmerStartBaseIndexVec(vector<int>& tmpKmerStartBaseIndexVecInSeq, int kmer_length, int jump_length)
    {
    	for(int tmp = 0; ; tmp++)
    	{
    		int tmpKmerStartBaseIndex = tmp * jump_length;
    		int tmpKmerEndBaseIndex = tmpKmerStartBaseIndex + kmer_length - 1;
    		if(tmpKmerEndBaseIndex < seqLength)
        		tmpKmerStartBaseIndexVecInSeq.push_back(tmpKmerStartBaseIndex);
    	    else
                return;
        }
    }

    void generateKmerStartBaseIndexVec(vector<int>& tmpKmerStartBaseIndexVecInSeq, int kmer_length)
    {
        int jump_length = kmer_length;
        this->generateKmerStartBaseIndexVec(tmpKmerStartBaseIndexVecInSeq, kmer_length, jump_length);
    }

    void querySeq_returnClassIdCountPairVec_old(vector< pair<valueT,int> >& tmpClassIdCountPairVec,
    	MulOth<keyT>* moth, int kmer_length, int jump_length)
    {
    	vector<int> tmpKmerStartBaseIndexVecInSeq;
        //cout << "start to generateKmerStartBaseIndexVec" << endl;
    	this->generateKmerStartBaseIndexVec(tmpKmerStartBaseIndexVecInSeq, kmer_length, jump_length);
    	vector<valueT> tmpKmerValueVec;
        //cout << "start to do queryKmerVec_baseIndexVec" << endl;
    	this->queryKmerVec_baseIndexVec(moth, tmpKmerValueVec, tmpKmerStartBaseIndexVecInSeq, kmer_length);
        //cout << "start to generate tmpClassIdCountPairVec" << endl;
    	this->kmerValueVec2classCountIdPairVec(tmpClassIdCountPairVec, tmpKmerValueVec);
    }

    void querySeq_returnClassIdCountPairVec(vector< pair<valueT,int> >& tmpClassIdCountPairVec,
        MulOth<keyT>* moth, int kmer_length, int jump_length)
    {
        vector<keyT> tmpKmerKeyVec;
        //cout << "start to generateKmerStartBaseIndexVec" << endl;
        if(jump_length < kmer_length) // overlapped
            this->generateKmerKeyVec_overlapped(tmpKmerKeyVec, kmer_length, jump_length);
        else
            this->generateKmerKeyVec_nonOverlapped(tmpKmerKeyVec, kmer_length, jump_length);
        vector<valueT> tmpKmerValueVec;
        //cout << "start to do queryKmerVec_baseIndexVec" << endl;
        this->queryKmerKeyVec(moth, tmpKmerValueVec, tmpKmerKeyVec);
        //cout << "start to generate tmpClassIdCountPairVec" << endl;
        this->kmerValueVec2classCountIdPairVec(tmpClassIdCountPairVec, tmpKmerValueVec);
    }

    void querySeq_returnKmerValueVec(vector<valueT>& tmpKmerValueVec,
        MulOth<keyT>* moth, int kmer_length, int jump_length)
    {
        vector<keyT> tmpKmerKeyVec;
        //cout << "start to generateKmerStartBaseIndexVec" << endl;
        if(jump_length < kmer_length) // overlapped
            this->generateKmerKeyVec_overlapped(tmpKmerKeyVec, kmer_length, jump_length);
        else
            this->generateKmerKeyVec_nonOverlapped(tmpKmerKeyVec, kmer_length, jump_length);
        ;
        //cout << "start to do queryKmerVec_baseIndexVec" << endl;
        this->queryKmerKeyVec(moth, tmpKmerValueVec, tmpKmerKeyVec);
    }

    void queryKmerKeyVec(MulOth<keyT>* moth, vector<valueT>& tmpKmerValueVec, vector<keyT>& tmpKmerKeyVec)
    {
        for(int tmp = 0; tmp < tmpKmerKeyVec.size(); tmp++)
        {
            //cout << "tmpKmerKeyVec[tmp]: " << tmpKmerKeyVec[tmp] << endl;
            valueT tmpKmerValue = moth->query(tmpKmerKeyVec[tmp]);
            //cout << "tmpKmerValue: " << tmpKmerValue << endl;
            tmpKmerValueVec.push_back(tmpKmerValue);
        }
    }

    void querySeq_returnClassIdCountPairVec(vector< pair<valueT,int> >& tmpClassIdCountPairVec,
        MulOth<keyT>* moth, int kmer_length)
    {
        int jump_length = kmer_length;
        this->querySeq_returnClassIdCountPairVec(tmpClassIdCountPairVec, moth, kmer_length, jump_length);
    }

    void kmerValueVec2classCountIdPairVec(vector< pair<valueT,int> >& tmpClassCountPairVec, 
    	vector<valueT>& tmpValueVec)
    {
    	for(int tmp = 0; tmp < tmpValueVec.size(); tmp++)
    	{
    		valueT tmpValue = tmpValueVec[tmp];
    		bool tmpAlreadyExists_bool = false;
    		int currentTmpClassCountPairVecSize = tmpClassCountPairVec.size();
    		for(int tmp2 = 0; tmp2 < currentTmpClassCountPairVecSize; tmp2++)
    		{
    			if(tmpClassCountPairVec[tmp2].first == tmpValue)
    			{
    				(tmpClassCountPairVec[tmp2].second)++;
    				tmpAlreadyExists_bool = true;
    				break;
    			}
    		}
    		if(!tmpAlreadyExists_bool)
    			tmpClassCountPairVec.push_back(pair<int,int>(tmpValue, 1));
    	}
    }

    void select_best_class(valueT& bestClassId, int& bestClassKmerCount,
    	vector< pair<valueT,int> >& tmpClassCountPairVec)
    {	
    	int tmpClassCount_best = -1;
    	valueT tmpClassId_best = 0;
    	for(int tmp = 0; tmp < tmpClassCountPairVec.size(); tmp++)
    	{
    		int tmpClassCount = tmpClassCountPairVec[tmp].second;
    		if(tmpClassCount > tmpClassCount_best)
    		{
    			tmpClassCount_best = tmpClassCount;
    			tmpClassId_best = tmpClassCountPairVec[tmp].first;
    		}
    	}
    	bestClassId = tmpClassId_best;
    	bestClassKmerCount = tmpClassCount_best;
    }

    void select_secondBest_class(valueT& secondBestClassId, int& secondBestClassKmerCount, 
    	valueT bestClassId, vector< pair<valueT,int> >& tmpClassCountPairVec)
    {
    	int tmpClassCount_secondBest = -1;
    	valueT tmpClassId_secondBest = 0;
    	for(int tmp = 0; tmp < tmpClassCountPairVec.size(); tmp++)
    	{
    		valueT tmpClassId = tmpClassCountPairVec[tmp].first;
    		if(tmpClassId != bestClassId)
    		{	
    			int tmpClassCount = tmpClassCountPairVec[tmp].second;
    			if(tmpClassCount > tmpClassCount_secondBest)
    			{
    				tmpClassCount_secondBest = tmpClassCount;
    				tmpClassId_secondBest = tmpClassId;
    			}
    		}
    	}
    	secondBestClassId = tmpClassId_secondBest;
    	secondBestClassKmerCount = tmpClassCount_secondBest;
    }

    void kmerValueVec2branchIndexCountPairVec(vector< pair<int,int> >& tmpHitBranchIndexCountPairVec, 
    	vector<valueT>& tmpValueVec, int totalBranchNum)
    {
    	for(int tmp = 0; tmp < tmpValueVec.size(); tmp++)
    	{
    		valueT tmpValue = tmpValueVec[tmp];
    		vector<int> tmpHitBranchIndexVec;
    		this->getBranchHitVec(totalBranchNum, tmpHitBranchIndexVec, tmpValue);
    		for(int tmp2 = 0; tmp2 < tmpHitBranchIndexVec.size(); tmp2++)
    		{
    			int tmpHitBranchIndex = tmpHitBranchIndexVec[tmp2]; // tmp hit branch
    			int currentTmpHitBranchIndexCountPairVecSize = tmpHitBranchIndexCountPairVec.size();
    			bool tmpBranchAlreadyExists_bool = false;
    			for(int tmp3 = 0; tmp3 < currentTmpHitBranchIndexCountPairVecSize; tmp3++)
    			{
    				if(tmpHitBranchIndexCountPairVec[tmp3].first == tmpHitBranchIndex)
    				{
    					(tmpHitBranchIndexCountPairVec[tmp3].second)++;
    					tmpBranchAlreadyExists_bool = true;
    					break;
    				}
    			}
    			if(!tmpBranchAlreadyExists_bool)
    				tmpHitBranchIndexCountPairVec.push_back(pair<int,int>(tmpHitBranchIndex, 1));
    		}
    	}
    }

    bool branchHitOrNot(valueT tmpValue, int tmpBranchIndex)
    {
    	int toAndInt = 1;
    	toAndInt<<=tmpBranchIndex;
    	if(tmpValue & toAndInt)
    		return true;
    	else
    		return false;
    }

    void getBranchHitVec(int totalBranchNum, vector<int>& hitBranchIndexVec, valueT tmpValue)
    {
    	for(int tmp = 0; tmp < totalBranchNum; tmp++)
    	{
    		bool tmpBranchHitOrNotBool = this->branchHitOrNot(tmpValue, tmp);
    		if(tmpBranchHitOrNotBool)
    			hitBranchIndexVec.push_back(tmp);
    	}
    }

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
                    {
                        return ret;
                    }
                }
                //cout << "error in convertKmer2Key" << endl;
                //exit(1);
        }
        //cout << "error in convertKmer2Key" << endl;
        //exit(1);        
    }

    // bool kmerCharArray2key(char* tmpChar, keyT* tmpK, int kmer_length)
    // {

    // }
};
#endif