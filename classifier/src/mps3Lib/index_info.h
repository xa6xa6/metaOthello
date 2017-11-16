// This file is provided as a part of MetaOthello. Please refer to LICENSE.TXT for the 'License'
#ifndef INDEX_INFO_H
#define INDEX_INFO_H

#include <stdlib.h>
#include <string>
#include <string.h>
#include <map>

#include "read_block_test.h"
#include "otherFunc.h"

using namespace std;

int INDEX_KMER_LENGTH = 14;

int baseChar2intArray[26] = {0, 100, 1, 100, 100, 100, 2,
			100, 100, 100, 100, 100, 100, 100,
			100, 100, 100, 100, 100, 3, 
			100, 100, 100, 100, 100, 100};

int baseCharCount2intArray[14][26] = {0};


class Index_Info
{
private:
	string chromString;
	unsigned int genomeLength;

	int chromNum;

	vector<string> chrNameStr; // size = chromNum
	vector<int> chromLength; // size = chromNum
	vector<string> chromStr;
	vector<unsigned int> chrEndPosInGenome;

	map<string, int> chrNameMap;
	//map<string, int>::iterator chrNameMapIter;

	int secondLevelIndexNormalSize;// = 3000000;
	vector<int> secondLevelIndexPartsNum;
	int secondLevelIndexPartsNumSum;

	set<int> invalidSecondLevelIndexNOset;
	//vector<int> secondLevelIndexLengthVec;

	unsigned int null_num; // 2654911540 for mm9_noRandom genome
	unsigned int indexSize; //2654911539  //sequence length + 1, the length of sa-lcp-down-next 

	vector<int> chrNameIndexArray;
	int chrNameIndexIntervalSize;
	map<int, set<int> > chrNameIndexArrayMap; // index_chrNameIndexArray, set<chrNameInt>
public:
	
	void initiateChrNameIndexArray(int intervalSize)
	{
		chrNameIndexIntervalSize = intervalSize;
		int chrNameIndexArraySize = (chrEndPosInGenome[chromNum-1])/chrNameIndexIntervalSize + 1;
		for(int tmp = 0; tmp < chrNameIndexArraySize; tmp++)
		{
			chrNameIndexArray.push_back(-2);
		}

		int indexArray_start_1 = 1/chrNameIndexIntervalSize;
		int indexArray_end_1 = chrEndPosInGenome[0]/chrNameIndexIntervalSize;
		for(int tmp = indexArray_start_1; tmp <= indexArray_end_1; tmp++)
		{
			chrNameIndexArray[tmp] = 0;
		}

		for(int tmp = 1; tmp < chromNum; tmp++)
		{
			int indexArray_start_tmp = (chrEndPosInGenome[tmp-1]+2)/chrNameIndexIntervalSize;
			int indexArray_end_tmp = chrEndPosInGenome[tmp]/chrNameIndexIntervalSize;
			for(int tmp2 = indexArray_start_tmp; tmp2 <= indexArray_end_tmp; tmp2++)
			{
				if(chrNameIndexArray[tmp2] == -2)
				{
					chrNameIndexArray[tmp2] = tmp;
				}
				else if(chrNameIndexArray[tmp2] == -1)
				{
					map<int, set<int> >::iterator mapIter_found = chrNameIndexArrayMap.find(tmp2);
					(mapIter_found->second).insert(tmp);
				}
				else
				{
					set<int> tmpIntSet;
					tmpIntSet.insert(tmp-1);
					tmpIntSet.insert(tmp);
					chrNameIndexArrayMap.insert(pair<int, set<int> >(tmp2, tmpIntSet));
					chrNameIndexArray[tmp2] = -1;
				}
			}
		}
	}

	/*
	void getChrLocation(unsigned int locationInWholeGenome, unsigned int *chr_name_int, unsigned int *chr_local_location)
	{
		#ifdef CAL_TIME
		getChrLocation_begin = clock();
		#endif
		if(locationInWholeGenome <= chrEndPosInGenome[0])
		{
			(*chr_name_int) = 0;
			*chr_local_location = locationInWholeGenome;
		}
		else
		{
			for(int tmp = 1; tmp < chrEndPosInGenome.size(); tmp++)
			{
				if( (locationInWholeGenome >= chrEndPosInGenome[tmp-1] + 2) 
					&& (locationInWholeGenome <= chrEndPosInGenome[tmp]) )
				{
					*chr_name_int = tmp;
					*chr_local_location = locationInWholeGenome - chrEndPosInGenome[tmp-1] - 2;
				}
				else
				{
					continue;
				}
			}
		}
		#ifdef CAL_TIME
		getChrLocation_end = clock();
		getChrLocation_cost = getChrLocation_cost + getChrLocation_end - getChrLocation_begin;
		#endif		
	}*/

	void getChrLocation(unsigned int locationInWholeGenome, unsigned int& chr_name_int, unsigned int& chr_local_location)
	{
		int chrNameInt = this->getChr(locationInWholeGenome);
		if(chrNameInt == 0)
		{
			(chr_name_int) = 0;
			(chr_local_location) = locationInWholeGenome;
		}
		else
		{
			(chr_name_int) = chrNameInt;
			(chr_local_location) = locationInWholeGenome - chrEndPosInGenome[chrNameInt - 1] -2;
		}
	}


	void getChrLocation(unsigned int locationInWholeGenome, unsigned int *chr_name_int, unsigned int *chr_local_location)
	{
		int chrNameInt = this->getChr(locationInWholeGenome);
		if(chrNameInt == 0)
		{
			(*chr_name_int) = 0;
			(*chr_local_location) = locationInWholeGenome;
		}
		else
		{
			(*chr_name_int) = chrNameInt;
			(*chr_local_location) = locationInWholeGenome - chrEndPosInGenome[chrNameInt - 1] -2;
		}
	}

	int getChrLocation_withChrNameInt(unsigned int locationInWholeGenome,
		int chr_name_int)
	{
		if(chr_name_int == 0)
		{
			return locationInWholeGenome;
		}
		else
		{
			return (locationInWholeGenome - chrEndPosInGenome[chr_name_int-1] - 2);
		}
	}

	/*int getChr(unsigned int locationInWholeGenome)
	{
		#ifdef CAL_TIME
		getChrLocation_begin = clock();
		#endif		
		int chrInt;
		if(locationInWholeGenome <= chrEndPosInGenome[0])
		{
			chrInt = 0;
		}
		else
		{
			for(int tmp = 1; tmp < chromNum; tmp++)
			{
				if( (locationInWholeGenome >= chrEndPosInGenome[tmp-1] + 2) 
					&& (locationInWholeGenome <= chrEndPosInGenome[tmp]) )
				{
					chrInt = tmp;
					break;
				}
				else
				{
					continue;
				}				
			}
		}
		#ifdef CAL_TIME
		getChrLocation_end = clock();
		getChrLocation_cost = getChrLocation_cost + getChrLocation_end - getChrLocation_begin;
		#endif				
		return chrInt;
	}*/

	int getChr(unsigned int locationInWholeGenome)
	{
		//cout << endl << "start to run getChr " << endl;
		int chrNameIndexArray_index = locationInWholeGenome/chrNameIndexIntervalSize;
		int index_chrNameVec = chrNameIndexArray[chrNameIndexArray_index];
		
		int getChrInt = -3;
		if(index_chrNameVec == -1)
		{
			map<int, set<int> >::iterator mapIter_found = 
				chrNameIndexArrayMap.find(chrNameIndexArray_index);

			if(mapIter_found == chrNameIndexArrayMap.end())
			{
				//cout << "error in get chr -- search in chrNameIndexArrayMap" << endl;
				exit(1);
			}
			else
			{
				//cout << "locationInWholeGenome: " << locationInWholeGenome << endl;
				for(set<int>::iterator setIter = (mapIter_found->second).begin(); 
					setIter != (mapIter_found->second).end(); setIter ++) 
				{
					int tmpIndex_chrNameVec = (*setIter);
					//cout << "tmpIndex_chrNameVec: " << tmpIndex_chrNameVec << endl;
					int tmpChr_startPos, tmpChr_endPos;
					if(tmpIndex_chrNameVec == 0)
					{
						tmpChr_startPos = 1;
						tmpChr_endPos = chrEndPosInGenome[0];
					}
					else
					{
						tmpChr_startPos = chrEndPosInGenome[tmpIndex_chrNameVec-1]+1;// fix me: should be set as +2 ?
						tmpChr_endPos = chrEndPosInGenome[tmpIndex_chrNameVec];
					}
					//cout << "tmpChr_pos: " << tmpChr_startPos << " ~ " << tmpChr_endPos << endl;
					if( (locationInWholeGenome >= tmpChr_startPos) && (locationInWholeGenome <= tmpChr_endPos) )
					{
						getChrInt = tmpIndex_chrNameVec;
						break;
					}
				}
				//cout << "some error in getChr -- 1" << endl;
			}
		}
		else if(index_chrNameVec == -2)
		{
			cout << "error in getChr ..." << endl;
			exit(1);
		}
		else
		{
			getChrInt = index_chrNameVec;
		}

		if(getChrInt == -3)
		{
			cout << "getChrInt == -3" << endl;
			cout << "locationInWholeGenome: " << locationInWholeGenome << endl;
			cout << "chrNameIndexArray_index: " << chrNameIndexArray_index << endl;
			cout << "index_chrNameVec: " << index_chrNameVec << endl;
			exit(1);
		}
		return getChrInt;
	}

	char getCharInChromosome(int chrNameInt, int chrMapPos)
	{
		return chromStr[chrNameInt].at(chrMapPos-1);
	}

	char getCharInWholeGenome(unsigned int wholeGenomeMapPos)
	{
		return chromString.at(wholeGenomeMapPos-1);
	}

	bool returnInvalidSecondLevelIndexNOset_find_bool(int value_to_find)
	{
		return ((invalidSecondLevelIndexNOset).find(value_to_find)
							!= (invalidSecondLevelIndexNOset).end()
							);
	}
	void insert2invalidSecondLevelIndexNOset(int value_to_insert)
	{
		invalidSecondLevelIndexNOset.insert(value_to_insert);
	} 

	int returnSecondLevelIndexNormalSize()
	{
		return secondLevelIndexNormalSize;
	}
	int returnSecondLevelIndexPartsNum(int secondLevelIndexPartsNum_index)
	{
		return secondLevelIndexPartsNum[secondLevelIndexPartsNum_index];
	}
	int returnSecondLevelIndexPartsNumSum()
	{
		return secondLevelIndexPartsNumSum;
	}
	unsigned int returnNull_num()
	{
		return null_num;
	}
	unsigned int returnIndexSize()
	{
		return indexSize;
	}
	int returnChrNameStrSize()
	{
		return chrNameStr.size();
	}
	const string& returnChrNameStr(int chrName_index)
	{
		return chrNameStr[chrName_index];
	}
	int returnChromLength(int chromLength_index)
	{
		return chromLength[chromLength_index];
	}
	const string& returnChromStr(int chromStr_index)
	{		
		return chromStr[chromStr_index];
	}
	//const string& 
	string returnChromStrSubstr(int chromStr_index, int start_pos, int seq_length)
	{
		int tmpChromLength = this->returnChromLength(chromStr_index);
		if((start_pos <= 0)||(start_pos + seq_length - 1 > tmpChromLength))
		{
			string tmpStr = "";
			for(int tmp = 0; tmp < seq_length; tmp++)
				tmpStr += "N";
			return tmpStr;
		}

		return chromStr[chromStr_index].substr(start_pos - 1, seq_length);
	}
	unsigned int returnChrEndPosInGenome(int chrEndPosInGenome_index)
	{
		return chrEndPosInGenome[chrEndPosInGenome_index];
	}
	void initiate()
	{
		//(indexInfo->chromStr).push_back((indexInfo->returnChromString()).substr(0, (indexInfo->chrEndPosInGenome)[0]+1));
		(chromStr).push_back(this->returnChromStringSubstr(1, (chrEndPosInGenome)[0]+1));
		(chromLength).push_back(((chrEndPosInGenome)[0]+1));
		for(int tmp = 1; tmp < this->returnChromNum(); tmp++)
		{
			//chromStr[tmp] = 
			//(indexInfo->chromStr).push_back((indexInfo->returnChromString()).substr((indexInfo->chrEndPosInGenome)[tmp-1]+2, 
			//	(indexInfo->chrEndPosInGenome)[tmp]-(indexInfo->chrEndPosInGenome)[tmp-1]-1));	
			(chromStr).push_back(
				this->returnChromStringSubstr((chrEndPosInGenome)[tmp-1]+3,
				(chrEndPosInGenome)[tmp]-(chrEndPosInGenome)[tmp-1]-1));// (indexInfo->chrEndPosInGenome)[tmp]-(indexInfo->chrEndPosInGenome)[tmp-1]-1);

			(chromLength).push_back(((chrEndPosInGenome)[tmp]-(chrEndPosInGenome)[tmp-1]-1));
		}		
	}

	void initiate_withoutLoadingSeq()
	{
		(chromLength).push_back(((chrEndPosInGenome)[0]+1));
		for(int tmp = 1; tmp < this->returnChromNum(); tmp++)
			(chromLength).push_back(((chrEndPosInGenome)[tmp]-(chrEndPosInGenome)[tmp-1]-1));	
		this->initiateChrNameIndexArray(1000);
	}

	void assignTranscriptIDarray(int* transcriptIDarray)
	{
		//for(int tmp = 1; tmp <= (chrEndPosInGenome)[0])

	}

	void assignGeneIDarray(int* geneIDarray)
	{}

	unsigned int returnChromStringLength()
	{
		return chromString.length();
	}
	void readGenome(char* chrom)
	{
		chromString = chrom;
	}
	const string& returnChromString()
	{
		return chromString;
	}

	string returnChromStringSubstr(unsigned int start_pos, unsigned int string_length)
	{
		return chromString.substr(start_pos - 1, string_length);
	}

	unsigned int returnGenomeLength()
	{
		return genomeLength;
	}
	int returnChromNum()
	{
		return chromNum;
	}


	int returnFlankStringCase(int chrNameInt, int chromPos_doner, int chromPos_acceptor)
	{
		string tmpFlankString = this->returnFlankString(chrNameInt, chromPos_doner, chromPos_acceptor);
		int tmpCase = this->returnFlankStringCaseFromFlankString(tmpFlankString);
		return tmpCase;
	}

	int returnFlankStringCaseFromFlankString(const string& flankString)
	{
		if(flankString == "GTAG")
		{
			return 5;
		}
		else if(flankString == "CTAC")
		{
			return 6;
		}
		else if(flankString == "ATAC")
		{
			return 1;
		}
		else if(flankString == "GTAT")
		{
			return 2;
		}
		else if(flankString == "CTGC")
		{
			return 3;
		}
		else if(flankString == "GCAG")
		{
			return 4;
		}
		else
		{
			return 0;
		}
	}

	string returnFlankStringStrand(int chrNameInt, int chromPos_doner, int chromPos_acceptor)
	{
		string tmpFlankString = this->returnFlankString(chrNameInt, chromPos_doner, chromPos_acceptor);
		return this->returnStrand_withFlankString(tmpFlankString);
	}

	string returnStrand_withFlankString(string& tmpFlkStr)
	{
		if((tmpFlkStr == "GTAG")||(tmpFlkStr == "GCAG")||(tmpFlkStr == "ATAC"))
			return "+";
		else if((tmpFlkStr == "CTAC")||(tmpFlkStr == "CTGC")||(tmpFlkStr == "GTAT"))
			return "-";
		else
			return "N";
	}

	string returnFlankString(int chrNameInt, int chromPos_doner, int chromPos_acceptor)
	{
		// if(chromPos_doner > chromPos_acceptor) // normal splice junction
		// {

		// }
		// else // backsplice junction
		// {

		// }
		return (chromStr[chrNameInt]).substr(chromPos_doner, 2) 
			+ (chromStr[chrNameInt]).substr(chromPos_acceptor-3, 2);
	}

	string returnBackSpliceFlankString(int chrNameInt, int chromPos_doner, int chromPos_acceptor)
	{
		if((chromPos_doner + 2 > chromLength[chrNameInt])||(chromPos_acceptor - 2 < 1))
		{
			return "XXXX";
		}
		else
			return returnFlankString(chrNameInt, chromPos_doner, chromPos_acceptor);
	}

	string returnFusionAnchorSeq_oneSide(int chrNameInt, int breakPointPos,
		bool gene1_or_gene2_bool, string& strand, int anchorSize)
	{
		int anchorSeq_startPos;//, anchorSeq_endPos;
		if(gene1_or_gene2_bool)
		{
			if(strand == "+")
				anchorSeq_startPos = breakPointPos - anchorSize + 1;
			else if(strand == "-")
				anchorSeq_startPos = breakPointPos;
			else
			{
				cout << "incorrectStrand in returnFusionAnchorSeq_oneSide: " << strand << endl;
				exit(1);
			}
		}
		else
		{
			if(strand == "+")
				anchorSeq_startPos = breakPointPos;
			else if(strand == "-")
				anchorSeq_startPos = breakPointPos - anchorSize + 1;
			else
			{
				cout << "incorrectStrand in returnFusionAnchorSeq_oneSide: " << strand << endl;
				exit(1);
			}
		}
		string tmpFusionAnchorSeq = this->returnChromStrSubstr(chrNameInt, anchorSeq_startPos, anchorSize);
		if(strand == "+")
			return tmpFusionAnchorSeq;
		else if(strand == "-")
			return covertStringToReverseComplement(tmpFusionAnchorSeq);
		else
		{
			cout << "incorrectStrand in returnFusionAnchorSeq_oneSide: " << strand << endl;
			exit(1);
		}
	}

	string returnFusionAnchorSeq(int chrNameInt_1, int breakPointPos_1, string& strand_1,
		int chrNameInt_2, int breakPointPos_2, string& strand_2, int anchorSize_inEachSide)
	{
		string anchorSeq_gene1 = this->returnFusionAnchorSeq_oneSide(
			chrNameInt_1, breakPointPos_1, true, strand_1, anchorSize_inEachSide);
		string anchorSeq_gene2 = this->returnFusionAnchorSeq_oneSide(
			chrNameInt_2, breakPointPos_2, false, strand_2, anchorSize_inEachSide);
		return (anchorSeq_gene1 + anchorSeq_gene2);
	}

	string returnFusionAnchorSeq(string& chrNameStr_1, int breakPointPos_1, string& strand_1,
		string& chrNameStr_2, int breakPointPos_2, string& strand_2, int anchorSize_inEachSide)
	{
		int chrNameInt_1 = this->convertStringToInt(chrNameStr_1);
		int chrNameInt_2 = this->convertStringToInt(chrNameStr_2);
		return this->returnFusionAnchorSeq(chrNameInt_1, breakPointPos_1, strand_1,
				chrNameInt_2, breakPointPos_2, strand_2, anchorSize_inEachSide);
	}

	bool try_returnRawFusionJuncFlankString(
		int leftBreakPoint_chrNameInt, int rightBreakPoint_chrNameInt,
		int leftBreakPoint_chrPos, int rightBreakPoint_chrPos,
		bool NorOrRcm_bool_1, bool NorOrRcm_bool_2, bool breakPointInEnd1OrEnd2_bool,
		string& generatedRawFusionJuncFlankString)
	{
		int chromLength_left = this->returnChromLength(leftBreakPoint_chrNameInt);
		int chromLength_right = this->returnChromLength(rightBreakPoint_chrNameInt);

		if(NorOrRcm_bool_1 && NorOrRcm_bool_2) // For For
		{
			if((leftBreakPoint_chrPos < 0)||(leftBreakPoint_chrPos+2 > chromLength_left)
				||(rightBreakPoint_chrPos-3 < 0)||(rightBreakPoint_chrPos-3+2 > chromLength_right))
				return false;
			generatedRawFusionJuncFlankString 
				= (chromStr[leftBreakPoint_chrNameInt]).substr(leftBreakPoint_chrPos, 2)
					+ (chromStr[rightBreakPoint_chrNameInt]).substr(rightBreakPoint_chrPos-3, 2);
			return true;
		}
		else if((!NorOrRcm_bool_1) && (!NorOrRcm_bool_2)) // Rev Rev
		{
			if((leftBreakPoint_chrPos < 0)||(leftBreakPoint_chrPos+2 > chromLength_left)
				||(rightBreakPoint_chrPos-3 < 0)||(rightBreakPoint_chrPos-3+2 > chromLength_right))
				return false;
			generatedRawFusionJuncFlankString 
				= (chromStr[leftBreakPoint_chrNameInt]).substr(leftBreakPoint_chrPos, 2)
					+ (chromStr[rightBreakPoint_chrNameInt]).substr(rightBreakPoint_chrPos-3, 2);
			return true;
		}
		else if((!NorOrRcm_bool_1) && NorOrRcm_bool_2) // Rev for
		{
			if(breakPointInEnd1OrEnd2_bool)
			{
				if((leftBreakPoint_chrPos-3 < 0)||(leftBreakPoint_chrPos-3+2 > chromLength_left)
					||(rightBreakPoint_chrPos-3 < 0)||(rightBreakPoint_chrPos-3+2 > chromLength_right))
					return false;
				generatedRawFusionJuncFlankString 
					= (chromStr[leftBreakPoint_chrNameInt]).substr(leftBreakPoint_chrPos-3, 2)
						+ (chromStr[rightBreakPoint_chrNameInt]).substr(rightBreakPoint_chrPos-3, 2);
				return true;
			}
			else
			{	
				if((leftBreakPoint_chrPos < 0)||(leftBreakPoint_chrPos+2 > chromLength_left)
					||(rightBreakPoint_chrPos < 0)||(rightBreakPoint_chrPos+2 > chromLength_right))
					return false;
				generatedRawFusionJuncFlankString 
					= (chromStr[leftBreakPoint_chrNameInt]).substr(leftBreakPoint_chrPos, 2)
						+ (chromStr[rightBreakPoint_chrNameInt]).substr(rightBreakPoint_chrPos, 2);
				return true;
			}
		}
		else
		{
			cout << "error in returnRawFusionJuncFlankString .. in index_info.h" << endl;
			exit(1);
		}
		/*else if(NorOrRcm_bool_1 && (!NorOrRcm_bool_2)) // For Rev
		{
			return (chromStr[leftBreakPoint_chrNameInt]).substr(leftBreakPoint_chrPos, 2)
				+ (chromStr[rightBreakPoint_chrNameInt]).substr(rightBreakPoint_chrPos, 2);
		}
		else // Rev For
		{
			return (chromStr[leftBreakPoint_chrNameInt]).substr(leftBreakPoint_chrPos-3, 2)
				+ (chromStr[rightBreakPoint_chrNameInt]).substr(rightBreakPoint_chrPos-3, 2);
		}*/
	}

	string returnFormattedFusionJuncFlankString(int leftBreakPoint_chrNameInt, int rightBreakPoint_chrNameInt, int leftBreakPoint_chrPos, 
		int rightBreakPoint_chrPos, bool fusionAtUpstreamEndOrDownstreamEnd_bool, bool clippedSeg_forMapOrRcmMap_bool)
	{
		int chromLength_left = this->returnChromLength(leftBreakPoint_chrNameInt);
		int chromLength_right = this->returnChromLength(rightBreakPoint_chrNameInt);
		if(fusionAtUpstreamEndOrDownstreamEnd_bool) // case 2,5,10,11
		{
			if(clippedSeg_forMapOrRcmMap_bool) // case 2,5
			{
				if((leftBreakPoint_chrPos + 1 < 1)||(leftBreakPoint_chrPos + 2 > chromLength_left)
					||(rightBreakPoint_chrPos - 1 < 1)||(rightBreakPoint_chrPos - 2 > chromLength_right))
					return "XXXX";				
				string halfFlankString_left = this->returnChromStrSubstr(leftBreakPoint_chrNameInt, leftBreakPoint_chrPos + 1, 2);
				string halfFlankString_right = this->returnChromStrSubstr(rightBreakPoint_chrNameInt, rightBreakPoint_chrPos - 2, 2);
				return (halfFlankString_left + halfFlankString_right);
			}
			else // case 10,11
			{
				if((leftBreakPoint_chrPos - 2 < 1)||(leftBreakPoint_chrPos - 1 > chromLength_left)
					||(rightBreakPoint_chrPos - 2 < 1)||(rightBreakPoint_chrPos - 1 > chromLength_right))
					return "XXXX";				
				string halfFlankString_left = convertStringToReverseComplement(
					this->returnChromStrSubstr(leftBreakPoint_chrNameInt, leftBreakPoint_chrPos - 2, 2));
				string halfFlankString_right = this->returnChromStrSubstr(rightBreakPoint_chrNameInt, rightBreakPoint_chrPos - 2, 2);
				return (halfFlankString_left + halfFlankString_right);
			}
		}
		else // case 1,4,7,8
		{
			if(clippedSeg_forMapOrRcmMap_bool) // case 1,4
			{
				if((leftBreakPoint_chrPos + 1 < 1)||(leftBreakPoint_chrPos + 2 > chromLength_left)
					||(rightBreakPoint_chrPos - 2 < 1)||(rightBreakPoint_chrPos - 1 > chromLength_right))
					return "XXXX";				
				string halfFlankString_left = this->returnChromStrSubstr(leftBreakPoint_chrNameInt, leftBreakPoint_chrPos + 1, 2);
				string halfFlankString_right = this->returnChromStrSubstr(rightBreakPoint_chrNameInt, rightBreakPoint_chrPos - 2, 2);
				return (halfFlankString_left + halfFlankString_right);				
			}
			else // case 7,8
			{
				if((leftBreakPoint_chrPos + 1 < 1)||(leftBreakPoint_chrPos + 2 > chromLength_left)
					||(rightBreakPoint_chrPos + 1 < 1)||(rightBreakPoint_chrPos + 2 > chromLength_right))
					return "XXXX";				
				string halfFlankString_left = this->returnChromStrSubstr(leftBreakPoint_chrNameInt, leftBreakPoint_chrPos + 1, 2);
				string halfFlankString_right = convertStringToReverseComplement(
					this->returnChromStrSubstr(rightBreakPoint_chrNameInt, rightBreakPoint_chrPos + 1, 2));
				return (halfFlankString_left + halfFlankString_right);
			}
		}
	}

	string returnRawFusionJuncFlankString(
		int leftBreakPoint_chrNameInt, int rightBreakPoint_chrNameInt,
		int leftBreakPoint_chrPos, int rightBreakPoint_chrPos,
		bool NorOrRcm_bool_1, bool NorOrRcm_bool_2, bool breakPointInEnd1OrEnd2_bool)
	{
		//cout << "start to do returnRawFusionJuncFlankString ... " << endl; 
		int chromLength_left = this->returnChromLength(leftBreakPoint_chrNameInt);
		int chromLength_right = this->returnChromLength(rightBreakPoint_chrNameInt);

		// cout << "NorOrRcm_bool_1: " << NorOrRcm_bool_1 << endl;
		// cout << "NorOrRcm_bool_2: " << NorOrRcm_bool_2 << endl;
		// cout << "breakPointInEnd1OrEnd2_bool: " << breakPointInEnd1OrEnd2_bool << endl;
		// cout << "chromLength_left: " << chromLength_left << endl;
		// cout << "chromLength_right: " << chromLength_right << endl;
		// cout << "leftBreakPoint_chrNameInt: " << leftBreakPoint_chrNameInt << endl;
		// cout << "rightBreakPoint_chrNameInt: " << rightBreakPoint_chrNameInt << endl;
		// cout << "leftBreakPoint_chrPos: " << leftBreakPoint_chrPos << endl;
		// cout << "rightBreakPoint_chrPos: " << rightBreakPoint_chrPos << endl;


		if(NorOrRcm_bool_1 && NorOrRcm_bool_2) // For For
		{
			if((leftBreakPoint_chrPos < 0)||(leftBreakPoint_chrPos+2 > chromLength_left)
				||(rightBreakPoint_chrPos-3 < 0)||(rightBreakPoint_chrPos-3+2 > chromLength_right))
				return "XXXX";
			else
				return (chromStr[leftBreakPoint_chrNameInt]).substr(leftBreakPoint_chrPos, 2)
					+ (chromStr[rightBreakPoint_chrNameInt]).substr(rightBreakPoint_chrPos-3, 2);
		}
		else if((!NorOrRcm_bool_1) && (!NorOrRcm_bool_2)) // Rev Rev
		{
			if((leftBreakPoint_chrPos < 0)||(leftBreakPoint_chrPos+2 > chromLength_left)
				||(rightBreakPoint_chrPos-3 < 0)||(rightBreakPoint_chrPos-3+2 > chromLength_right))
				return "XXXX";
			else
				return (chromStr[leftBreakPoint_chrNameInt]).substr(leftBreakPoint_chrPos, 2)
					+ (chromStr[rightBreakPoint_chrNameInt]).substr(rightBreakPoint_chrPos-3, 2);
		}
		else if((!NorOrRcm_bool_1) && NorOrRcm_bool_2) // Rev for
		{
			if(breakPointInEnd1OrEnd2_bool)
			{
				if((leftBreakPoint_chrPos-3 < 0)||(leftBreakPoint_chrPos-3+2 > chromLength_left)
					||(rightBreakPoint_chrPos-3 < 0)||(rightBreakPoint_chrPos-3+2 > chromLength_right))
					return "XXXX";
				else
					return (chromStr[leftBreakPoint_chrNameInt]).substr(leftBreakPoint_chrPos-3, 2)
						+ (chromStr[rightBreakPoint_chrNameInt]).substr(rightBreakPoint_chrPos-3, 2);
			}
			else
			{
				if((leftBreakPoint_chrPos < 0)||(leftBreakPoint_chrPos+2 > chromLength_left)
					||(rightBreakPoint_chrPos < 0)||(rightBreakPoint_chrPos+2 > chromLength_right))
					return "XXXX";			
				else
					return (chromStr[leftBreakPoint_chrNameInt]).substr(leftBreakPoint_chrPos, 2)
						+ (chromStr[rightBreakPoint_chrNameInt]).substr(rightBreakPoint_chrPos, 2);
			}
		}
		else
		{
			cout << "error in returnRawFusionJuncFlankString .. in index_info.h" << endl;
			exit(1);
		}
		/*else if(NorOrRcm_bool_1 && (!NorOrRcm_bool_2)) // For Rev
		{
			return (chromStr[leftBreakPoint_chrNameInt]).substr(leftBreakPoint_chrPos, 2)
				+ (chromStr[rightBreakPoint_chrNameInt]).substr(rightBreakPoint_chrPos, 2);
		}
		else // Rev For
		{
			return (chromStr[leftBreakPoint_chrNameInt]).substr(leftBreakPoint_chrPos-3, 2)
				+ (chromStr[rightBreakPoint_chrNameInt]).substr(rightBreakPoint_chrPos-3, 2);
		}*/
	}

	string returnRawFusionJuncFlankString(
		string& leftBreakPoint_chrNameStr, string& rightBreakPoint_chrNameStr,
		int leftBreakPoint_chrPos, int rightBreakPoint_chrPos,
		bool NorOrRcm_bool_1, bool NorOrRcm_bool_2, Index_Info* indexInfo,
		bool breakPointInEnd1OrEnd2_bool)
	{
		int leftBreakPoint_chrNameInt = indexInfo->convertStringToInt(leftBreakPoint_chrNameStr);
		int rightBreakPoint_chrNameInt = indexInfo->convertStringToInt(rightBreakPoint_chrNameStr);
		return this->returnRawFusionJuncFlankString(
						leftBreakPoint_chrNameInt, rightBreakPoint_chrNameInt,
						leftBreakPoint_chrPos, rightBreakPoint_chrPos,
						NorOrRcm_bool_1, NorOrRcm_bool_2, breakPointInEnd1OrEnd2_bool);
	}

	/*string returnFusionJuncFlankString(
		int chrNameInt_1, int chromPos_1, int chrNameInt_2, int chromPos_2,
		string& strand_1, string& strand_2)
	{
		if((strand_1 == "+")&&(strand_2 == "+"))
			return (chromStr[chrNameInt_1]).substr(chromPos_1, 2) 
				+ (chromStr[chrNameInt_2]).substr(chromPos_2-3, 2);	
		else if((strand_1 == "-")&&(strand_2 == "-"))
		{
			return convertStringToReverseComplement(chromStr[chrNameInt_1].substr(chromPos_1-3, 2))
				+ covertStringToReverseComplement(chromStr[chrNameInt_2].substr(chromPos_2, 2));
		}
		else if((strand_1 == "+")&&(strand_2 == "-"))
		{
			return chromStr[chrNameInt_1].substr(chromPos_1, 2) 
				+ convertStringToReverseComplement(chromStr[chrNameInt_2].substr(chromPos_2, 2));
		}
		else // ((strand_1 == "-")&&(strand_2 == "+"))
		{
			return convertStringToReverseComplement(chromStr[chrNameInt_1].substr(chromPos_1-3, 2))
				+ chromStr[chrNameInt_2].substr(chromPos_2-3, 2);
		}
	}

	string returnFusionJuncFlankString(
		string& chrNameStr_1, int chrPos_1, string& chrNameStr_2, int chrPos_2,
		string& strand_1, string& strand_2)
	{
		int chrNameInt_1 = this->convertStringToInt(chrNameStr_1);
		int chrNameInt_2 = this->convertStringToInt(chrNameStr_2);
		string tmpFusionJuncFlankString = returnFusionJuncFlankString(
			chrNameInt_1, chrPos_1, chrNameInt_2, chrPos_2,
			strand_1, strand_2);
		return tmpFusionJuncFlankString;
	}*/

	string returnTwoBasesString(int chrNameInt, int firstChromPos)
	{
		return (chromStr[chrNameInt]).substr(firstChromPos-1, 2); 
			//+ (chromStr[chrNameInt]).substr(chromPos_acceptor-3, 2);
	}

	char returnOneBaseCharInGenome(int chrNameInt, int chromPos)
	{
		return (chromStr[chrNameInt]).at(chromPos-1);
	}

	string returnOneBaseStrInGenome(int chrNameInt, int chromPos)
	{
		return (chromStr[chrNameInt]).substr(chromPos-1, 1);
	}

	string returnGenomeSubstr(int chrNameInt, int startPosInChrom, int substrLength)
	{
		return (chromStr[chrNameInt]).substr(startPosInChrom - 1, substrLength);
	}

	string getReferenceGenomeSubstr(const string& chrName, int startPos, int endPos) // can be used in debugging on some special case
	{
		int chrNameInt = this->convertStringToInt(chrName);
		return (this->chromStr)[chrNameInt].substr(startPos-1, endPos - startPos + 1);
	}

	string getInvalidSecondLevelIndexNOstr()
	{
		string tmpStr = "invalidSecondLevelIndexNO: \n";
		for(set<int>::iterator setIter = invalidSecondLevelIndexNOset.begin(); setIter != invalidSecondLevelIndexNOset.end(); setIter ++)
		{
			tmpStr += int_to_str(*setIter);
			tmpStr += ",";
		} 
		tmpStr += "\n";
		return tmpStr;
	}

	Index_Info()
	{}

	Index_Info(ifstream& inputIndexInfoFile, ofstream& outputIndexInfoFile)
	{

		for(int tmp1 = 0; tmp1 < INDEX_KMER_LENGTH; tmp1++)
		{
			for(int tmp2 = 0; tmp2 < 26/*# of letters alphabet*/; tmp2++)
			{
				baseCharCount2intArray[tmp1][tmp2] = 100;
			}
		}
		int tmpBaseCount = 1;
		for(int tmp3 = 0; tmp3 < INDEX_KMER_LENGTH; tmp3++)
		{
			baseCharCount2intArray[tmp3][0] = 0*tmpBaseCount;
			baseCharCount2intArray[tmp3][2] = 1*tmpBaseCount;
			baseCharCount2intArray[tmp3][6] = 2*tmpBaseCount;
			baseCharCount2intArray[tmp3][19] = 3*tmpBaseCount;
			tmpBaseCount = 4*tmpBaseCount;
		}


		string s;
		string chromNumLine;
		string chromNameLine;
		string chromEndPosInGenomeLine;
		string secondLevelIndexSizeLine;
		string chrom2ndLevelIndexNumLine;
		outputIndexInfoFile << endl << " #####  index information:  #####" << endl;
		getline(inputIndexInfoFile, s);
		outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromNumLine);
		outputIndexInfoFile << chromNumLine << endl;
		getline(inputIndexInfoFile, s);
		outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromNameLine);
		outputIndexInfoFile << chromNameLine << endl;
		getline(inputIndexInfoFile, s);
		outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromEndPosInGenomeLine);
		outputIndexInfoFile << chromEndPosInGenomeLine << endl;
		getline(inputIndexInfoFile, s);
		outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, secondLevelIndexSizeLine);
		outputIndexInfoFile << secondLevelIndexSizeLine << endl;
		getline(inputIndexInfoFile, s);
		outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chrom2ndLevelIndexNumLine);
		outputIndexInfoFile << chrom2ndLevelIndexNumLine << endl;
		outputIndexInfoFile << "***************************************" << endl << endl;
		chromNum = atoi( (chromNumLine.substr(0, chromNumLine.length())).c_str() );
		//cout << "chromNum: " << chromNum << endl;

		int startSearchPos = 0;
		int foundSearchPos;
		string tmpChromNameStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromNameLine.find(",", startSearchPos);
			tmpChromNameStr = chromNameLine.substr(startSearchPos+1, foundSearchPos - 2 - startSearchPos - 1 + 1);
			chrNameStr.push_back(tmpChromNameStr);
			//cout << tmp+1 << " tmpChromNameStr: " << tmpChromNameStr << " strLen: " << tmpChromNameStr.length() << endl;
			startSearchPos = foundSearchPos + 1;
		}

		startSearchPos = 0;
		string tmpChromEndPosStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromEndPosInGenomeLine.find(",", startSearchPos);
			tmpChromEndPosStr = chromEndPosInGenomeLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			unsigned int tmpChromEndPos = strtoul(tmpChromEndPosStr.c_str(), NULL, 10);
			chrEndPosInGenome.push_back(tmpChromEndPos);
			//cout << tmp+1 << " tmpChromEndPos: " << tmpChromEndPos << endl;
			startSearchPos = foundSearchPos + 1;
		}		

		secondLevelIndexNormalSize = atoi( (secondLevelIndexSizeLine.substr(0, secondLevelIndexSizeLine.length())).c_str() );
		//cout << "secondLevelIndexNormalSize: " << secondLevelIndexNormalSize << endl;

		startSearchPos = 0;
		string tmpChrom2ndLevelIndexNumStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chrom2ndLevelIndexNumLine.find(",", startSearchPos);
			tmpChrom2ndLevelIndexNumStr = chrom2ndLevelIndexNumLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			int tmpChrom2ndLevelIndexNum = atoi(tmpChrom2ndLevelIndexNumStr.c_str());
			secondLevelIndexPartsNum.push_back(tmpChrom2ndLevelIndexNum);
			//cout << tmp+1 << "tmp2ndLevelIndexNum: " << tmpChrom2ndLevelIndexNum << endl;
			startSearchPos = foundSearchPos + 1;
		}

		indexSize = chrEndPosInGenome[chrEndPosInGenome.size()-1] + 2;
		//cout << "MAX: " << indexSize << endl;
		null_num = indexSize + 1;
		//cout << "NULL_NUM: " << null_num << endl;
		this->buildChrNameMap();
	}

	void initiate(ifstream& inputIndexInfoFile, ofstream& outputIndexInfoFile)
	{

		for(int tmp1 = 0; tmp1 < INDEX_KMER_LENGTH; tmp1++)
		{
			for(int tmp2 = 0; tmp2 < 26/*# of letters alphabet*/; tmp2++)
			{
				baseCharCount2intArray[tmp1][tmp2] = 100;
			}
		}
		int tmpBaseCount = 1;
		for(int tmp3 = 0; tmp3 < INDEX_KMER_LENGTH; tmp3++)
		{
			baseCharCount2intArray[tmp3][0] = 0*tmpBaseCount;
			baseCharCount2intArray[tmp3][2] = 1*tmpBaseCount;
			baseCharCount2intArray[tmp3][6] = 2*tmpBaseCount;
			baseCharCount2intArray[tmp3][19] = 3*tmpBaseCount;
			tmpBaseCount = 4*tmpBaseCount;
		}


		string s;
		string chromNumLine;
		string chromNameLine;
		string chromEndPosInGenomeLine;
		string secondLevelIndexSizeLine;
		string chrom2ndLevelIndexNumLine;
		outputIndexInfoFile << endl << " #####  index information:  #####" << endl;
		getline(inputIndexInfoFile, s);
		outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromNumLine);
		outputIndexInfoFile << chromNumLine << endl;
		getline(inputIndexInfoFile, s);
		outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromNameLine);
		outputIndexInfoFile << chromNameLine << endl;
		getline(inputIndexInfoFile, s);
		outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromEndPosInGenomeLine);
		outputIndexInfoFile << chromEndPosInGenomeLine << endl;
		getline(inputIndexInfoFile, s);
		outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, secondLevelIndexSizeLine);
		outputIndexInfoFile << secondLevelIndexSizeLine << endl;
		getline(inputIndexInfoFile, s);
		outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chrom2ndLevelIndexNumLine);
		outputIndexInfoFile << chrom2ndLevelIndexNumLine << endl;
		outputIndexInfoFile << "***************************************" << endl << endl;
		chromNum = atoi( (chromNumLine.substr(0, chromNumLine.length())).c_str() );
		//cout << "chromNum: " << chromNum << endl;

		int startSearchPos = 0;
		int foundSearchPos;
		string tmpChromNameStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromNameLine.find(",", startSearchPos);
			tmpChromNameStr = chromNameLine.substr(startSearchPos+1, foundSearchPos - 2 - startSearchPos - 1 + 1);
			chrNameStr.push_back(tmpChromNameStr);
			//cout << tmp+1 << " tmpChromNameStr: " << tmpChromNameStr << " strLen: " << tmpChromNameStr.length() << endl;
			startSearchPos = foundSearchPos + 1;
		}

		startSearchPos = 0;
		string tmpChromEndPosStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromEndPosInGenomeLine.find(",", startSearchPos);
			tmpChromEndPosStr = chromEndPosInGenomeLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			unsigned int tmpChromEndPos = strtoul(tmpChromEndPosStr.c_str(), NULL, 10);
			chrEndPosInGenome.push_back(tmpChromEndPos);
			//cout << tmp+1 << " tmpChromEndPos: " << tmpChromEndPos << endl;
			startSearchPos = foundSearchPos + 1;
		}		

		secondLevelIndexNormalSize = atoi( (secondLevelIndexSizeLine.substr(0, secondLevelIndexSizeLine.length())).c_str() );
		//cout << "secondLevelIndexNormalSize: " << secondLevelIndexNormalSize << endl;

		startSearchPos = 0;
		string tmpChrom2ndLevelIndexNumStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chrom2ndLevelIndexNumLine.find(",", startSearchPos);
			tmpChrom2ndLevelIndexNumStr = chrom2ndLevelIndexNumLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			int tmpChrom2ndLevelIndexNum = atoi(tmpChrom2ndLevelIndexNumStr.c_str());
			secondLevelIndexPartsNum.push_back(tmpChrom2ndLevelIndexNum);
			//cout << tmp+1 << "tmp2ndLevelIndexNum: " << tmpChrom2ndLevelIndexNum << endl;
			startSearchPos = foundSearchPos + 1;
		}

		indexSize = chrEndPosInGenome[chrEndPosInGenome.size()-1] + 2;
		//cout << "MAX: " << indexSize << endl;
		null_num = indexSize + 1;
		//cout << "NULL_NUM: " << null_num << endl;
		this->buildChrNameMap();
	}
	Index_Info(ifstream& inputIndexInfoFile)
	{

		for(int tmp1 = 0; tmp1 < INDEX_KMER_LENGTH; tmp1++)
		{
			for(int tmp2 = 0; tmp2 < 26/*# of letters alphabet*/; tmp2++)
			{
				baseCharCount2intArray[tmp1][tmp2] = 100;
			}
		}
		int tmpBaseCount = 1;
		for(int tmp3 = 0; tmp3 < INDEX_KMER_LENGTH; tmp3++)
		{
			baseCharCount2intArray[tmp3][0] = 0*tmpBaseCount;
			baseCharCount2intArray[tmp3][2] = 1*tmpBaseCount;
			baseCharCount2intArray[tmp3][6] = 2*tmpBaseCount;
			baseCharCount2intArray[tmp3][19] = 3*tmpBaseCount;
			tmpBaseCount = 4*tmpBaseCount;
		}


		string s;
		string chromNumLine;
		string chromNameLine;
		string chromEndPosInGenomeLine;
		string secondLevelIndexSizeLine;
		string chrom2ndLevelIndexNumLine;
		//outputIndexInfoFile << endl << " #####  index information:  #####" << endl;
		getline(inputIndexInfoFile, s);
		//outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromNumLine);
		//outputIndexInfoFile << chromNumLine << endl;
		getline(inputIndexInfoFile, s);
		//outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromNameLine);
		//outputIndexInfoFile << chromNameLine << endl;
		getline(inputIndexInfoFile, s);
		//outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromEndPosInGenomeLine);
		//outputIndexInfoFile << chromEndPosInGenomeLine << endl;
		getline(inputIndexInfoFile, s);
		//outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, secondLevelIndexSizeLine);
		//outputIndexInfoFile << secondLevelIndexSizeLine << endl;
		getline(inputIndexInfoFile, s);
		//outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chrom2ndLevelIndexNumLine);
		//outputIndexInfoFile << chrom2ndLevelIndexNumLine << endl;
		//outputIndexInfoFile << "***************************************" << endl << endl;
		chromNum = atoi( (chromNumLine.substr(0, chromNumLine.length())).c_str() );
		//cout << "chromNum: " << chromNum << endl;

		int startSearchPos = 0;
		int foundSearchPos;
		string tmpChromNameStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromNameLine.find(",", startSearchPos);
			tmpChromNameStr = chromNameLine.substr(startSearchPos+1, foundSearchPos - 2 - startSearchPos - 1 + 1);
			chrNameStr.push_back(tmpChromNameStr);
			//cout << tmp+1 << " tmpChromNameStr: " << tmpChromNameStr << " strLen: " << tmpChromNameStr.length() << endl;
			startSearchPos = foundSearchPos + 1;
		}

		startSearchPos = 0;
		string tmpChromEndPosStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromEndPosInGenomeLine.find(",", startSearchPos);
			tmpChromEndPosStr = chromEndPosInGenomeLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			unsigned int tmpChromEndPos = strtoul(tmpChromEndPosStr.c_str(), NULL, 10);
			chrEndPosInGenome.push_back(tmpChromEndPos);
			//cout << tmp+1 << " tmpChromEndPos: " << tmpChromEndPos << endl;
			startSearchPos = foundSearchPos + 1;
		}		

		secondLevelIndexNormalSize = atoi( (secondLevelIndexSizeLine.substr(0, secondLevelIndexSizeLine.length())).c_str() );
		//cout << "secondLevelIndexNormalSize: " << secondLevelIndexNormalSize << endl;

		startSearchPos = 0;
		string tmpChrom2ndLevelIndexNumStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chrom2ndLevelIndexNumLine.find(",", startSearchPos);
			tmpChrom2ndLevelIndexNumStr = chrom2ndLevelIndexNumLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			int tmpChrom2ndLevelIndexNum = atoi(tmpChrom2ndLevelIndexNumStr.c_str());
			secondLevelIndexPartsNum.push_back(tmpChrom2ndLevelIndexNum);
			//cout << tmp+1 << "tmp2ndLevelIndexNum: " << tmpChrom2ndLevelIndexNum << endl;
			startSearchPos = foundSearchPos + 1;
		}

		indexSize = chrEndPosInGenome[chrEndPosInGenome.size()-1] + 2;
		//cout << "MAX: " << indexSize << endl;
		null_num = indexSize + 1;
		//cout << "NULL_NUM: " << null_num << endl;
		this->buildChrNameMap();
	}

	Index_Info(string inputIndexInfoFilePath)
	{
		ifstream inputIndexInfoFile(inputIndexInfoFilePath.c_str());
		for(int tmp1 = 0; tmp1 < INDEX_KMER_LENGTH; tmp1++)
		{
			for(int tmp2 = 0; tmp2 < 26/*# of letters alphabet*/; tmp2++)
			{
				baseCharCount2intArray[tmp1][tmp2] = 100;
			}
		}
		int tmpBaseCount = 1;
		for(int tmp3 = 0; tmp3 < INDEX_KMER_LENGTH; tmp3++)
		{
			baseCharCount2intArray[tmp3][0] = 0*tmpBaseCount;
			baseCharCount2intArray[tmp3][2] = 1*tmpBaseCount;
			baseCharCount2intArray[tmp3][6] = 2*tmpBaseCount;
			baseCharCount2intArray[tmp3][19] = 3*tmpBaseCount;
			tmpBaseCount = 4*tmpBaseCount;
		}
		string s;
		string chromNumLine;
		string chromNameLine;
		string chromEndPosInGenomeLine;
		string secondLevelIndexSizeLine;
		string chrom2ndLevelIndexNumLine;
		//outputIndexInfoFile << endl << " #####  index information:  #####" << endl;
		getline(inputIndexInfoFile, s);
		//outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromNumLine);
		//outputIndexInfoFile << chromNumLine << endl;
		getline(inputIndexInfoFile, s);
		//outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromNameLine);
		//outputIndexInfoFile << chromNameLine << endl;
		getline(inputIndexInfoFile, s);
		//outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chromEndPosInGenomeLine);
		//outputIndexInfoFile << chromEndPosInGenomeLine << endl;
		getline(inputIndexInfoFile, s);
		//outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, secondLevelIndexSizeLine);
		//outputIndexInfoFile << secondLevelIndexSizeLine << endl;
		getline(inputIndexInfoFile, s);
		//outputIndexInfoFile << s << endl;
		getline(inputIndexInfoFile, chrom2ndLevelIndexNumLine);
		//outputIndexInfoFile << chrom2ndLevelIndexNumLine << endl;
		//outputIndexInfoFile << "***************************************" << endl << endl;
		chromNum = atoi( (chromNumLine.substr(0, chromNumLine.length())).c_str() );
		//cout << "chromNum: " << chromNum << endl;

		int startSearchPos = 0;
		int foundSearchPos;
		string tmpChromNameStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromNameLine.find(",", startSearchPos);
			tmpChromNameStr = chromNameLine.substr(startSearchPos+1, foundSearchPos - 2 - startSearchPos - 1 + 1);
			chrNameStr.push_back(tmpChromNameStr);
			//cout << tmp+1 << " tmpChromNameStr: " << tmpChromNameStr << " strLen: " << tmpChromNameStr.length() << endl;
			startSearchPos = foundSearchPos + 1;
		}

		startSearchPos = 0;
		string tmpChromEndPosStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chromEndPosInGenomeLine.find(",", startSearchPos);
			tmpChromEndPosStr = chromEndPosInGenomeLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			unsigned int tmpChromEndPos = strtoul(tmpChromEndPosStr.c_str(), NULL, 10);
			chrEndPosInGenome.push_back(tmpChromEndPos);
			//cout << tmp+1 << " tmpChromEndPos: " << tmpChromEndPos << endl;
			startSearchPos = foundSearchPos + 1;
		}		

		secondLevelIndexNormalSize = atoi( (secondLevelIndexSizeLine.substr(0, secondLevelIndexSizeLine.length())).c_str() );
		//cout << "secondLevelIndexNormalSize: " << secondLevelIndexNormalSize << endl;

		startSearchPos = 0;
		string tmpChrom2ndLevelIndexNumStr;
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			foundSearchPos = chrom2ndLevelIndexNumLine.find(",", startSearchPos);
			tmpChrom2ndLevelIndexNumStr = chrom2ndLevelIndexNumLine.substr(startSearchPos, foundSearchPos - 1 - startSearchPos + 1);
			int tmpChrom2ndLevelIndexNum = atoi(tmpChrom2ndLevelIndexNumStr.c_str());
			secondLevelIndexPartsNum.push_back(tmpChrom2ndLevelIndexNum);
			//cout << tmp+1 << "tmp2ndLevelIndexNum: " << tmpChrom2ndLevelIndexNum << endl;
			startSearchPos = foundSearchPos + 1;
		}

		indexSize = chrEndPosInGenome[chrEndPosInGenome.size()-1] + 2;
		//cout << "MAX: " << indexSize << endl;
		null_num = indexSize + 1;
		//cout << "NULL_NUM: " << null_num << endl;
		this->buildChrNameMap();
		inputIndexInfoFile.close();
	}	

	void buildChrNameMap()
	{
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			chrNameMap.insert(pair <string, int> (chrNameStr[tmp], tmp));
		}

	}

	int getSecondLevelIndexFromChrAndPos(int chrNameInt, int chrMapPos)
	{
		int tmpTimes = chrMapPos/secondLevelIndexNormalSize;
		int partsTimeBase = 0;
		for(int tmp = 0; tmp < chrNameInt; tmp++)
		{
			partsTimeBase += secondLevelIndexPartsNum[tmp];
		}
		return	(partsTimeBase + tmpTimes + 1); 
	}

	int getSecondLevelIndexFromChrStrAndPos(string chrNameStr, int chrMapPos)
	{
		int chrNameInt = this->convertStringToInt(chrNameStr);
		int tmpTimes = chrMapPos/secondLevelIndexNormalSize;
		int partsTimeBase = 0;
		for(int tmp = 0; tmp < chrNameInt; tmp++)
		{
			partsTimeBase += secondLevelIndexPartsNum[tmp];
		}
		return	(partsTimeBase + tmpTimes + 1); 
	} 

	int getChrPosFromSecondLevelIndexPos(int chrNameInt, int secondLevelIndexNum, int secondLevelIndexPos)
	{
		int partsTimeBase = 0;
		for(int tmp = 0; tmp < chrNameInt; tmp++)
		{
			partsTimeBase += secondLevelIndexPartsNum[tmp];
		}

		int tmpSecondLevelIndexNO = secondLevelIndexNum - partsTimeBase;
		return ( (tmpSecondLevelIndexNO-1) * secondLevelIndexNormalSize + secondLevelIndexPos);	
	}

	int getChrNameIntFromSecondLevelIndexNO(int secondLevelIndexNum)
	{
		int tmpIndexNOsum = 0;
		int tmp;
		for(tmp = 0; tmp < chrNameStr.size(); tmp++)
		{
			tmpIndexNOsum += secondLevelIndexPartsNum[tmp];
			if(tmpIndexNOsum >= secondLevelIndexNum)
				return tmp;
		}
	}


	unsigned int getWholeGenomeLocation(unsigned int chr_name_int, unsigned int locationInWholeGenome)
	{
		//cout << "in function chr_name_int: " << chr_name_int << endl;
		unsigned int chr_local_location;
		if(chr_name_int == 0)
		{
			chr_local_location = locationInWholeGenome;
		}
		else if(chr_name_int < chromNum)
		{
			//cout << "< chromNum chr_name_int: " << chr_name_int << endl;
			chr_local_location = locationInWholeGenome 
				+ chrEndPosInGenome[chr_name_int-1] + 2; 
			//cout << "chr_local_location: " << chr_local_location << endl;
		}
		else
		{
			cout << "chr_name_int error: " << chr_name_int << endl;
		}
		return chr_local_location;
	}

	unsigned int getWholeGenomeLocation(const string& chromNameStr, unsigned int locationInWholeGenome)
	{
		//cout << "in function chr_name_int: " << chr_name_int << endl;
		
		int chr_name_int = this->convertStringToInt(chromNameStr);
		
		unsigned int chr_local_location;
		if(chr_name_int == 0)
		{
			chr_local_location = locationInWholeGenome;
		}
		else if(chr_name_int < chromNum)
		{
			//cout << "< chromNum chr_name_int: " << chr_name_int << endl;
			chr_local_location = locationInWholeGenome 
				+ chrEndPosInGenome[chr_name_int-1] + 2; 
			//cout << "chr_local_location: " << chr_local_location << endl;
		}
		else
		{
			cout << "chr_name_int error: " << chr_name_int << endl;
		}
		return chr_local_location;
	}

	int convertStringToInt(const string& chrName)
	{
		map<string, int>::iterator chrNameMapIter;
		int chrNameInt = 1000;
		chrNameMapIter = chrNameMap.find(chrName);
		if(chrNameMapIter != chrNameMap.end())
		{
			chrNameInt = chrNameMapIter->second;
		}
		else
		{
			chrNameInt = -1;
			//cout << "...... chrom name error! ...... " << endl;
		}
		return chrNameInt;
	}

	void insertSNP2chromStr(string& SNPfilePath, ofstream& log_ofs)
	{
		cout << "start to extract SNPs from formatted SNPfile ! " << endl;
		log_ofs << "start to extract SNPs from formatted SNPfile ! " << endl;
		vector<int> SNP_chrNameIntVec;
		vector<int> SNP_chrPosVec;
		vector<string> SNP_refBaseVec;
		vector<string> SNP_alterBaseVec;
		ifstream SNP_ifs(SNPfilePath.c_str());
		int tmpSNP_num_valid = 0;
		int tmpSNP_num_invalid = 0;
		while(!SNP_ifs.eof())
		{
			string tmpSNPstr;
			getline(SNP_ifs, tmpSNPstr);
			if(tmpSNPstr == "")
				break;
			int tabLoc_1 = tmpSNPstr.find("\t");
			int tabLoc_2 = tmpSNPstr.find("\t", tabLoc_1 + 1);
			int tabLoc_3 = tmpSNPstr.find("\t", tabLoc_2 + 1);			
			string tmpChrNameStr = tmpSNPstr.substr(0, tabLoc_1);
			int tmpChrNameInt = this->convertStringToInt(tmpChrNameStr);
			string tmpPosStr = tmpSNPstr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
			int tmpPos = atoi(tmpPosStr.c_str());			
			string tmpRefBase = tmpSNPstr.substr(tabLoc_2 + 1, 1);
			string tmpAlterBase = tmpSNPstr.substr(tabLoc_3 + 1, 1);
			if(tmpChrNameInt < 0)
				tmpSNP_num_invalid ++;
			else
			{
				tmpSNP_num_valid ++;
				SNP_chrNameIntVec.push_back(tmpChrNameInt);
				SNP_chrPosVec.push_back(tmpPos);
				SNP_refBaseVec.push_back(tmpRefBase);
				SNP_alterBaseVec.push_back(tmpAlterBase);
			}
		}
		SNP_ifs.close();
		log_ofs << "tmpSNP_num_valid: " << tmpSNP_num_valid << endl;
		log_ofs << "tmpSNP_num_invalid: " << tmpSNP_num_invalid << endl;
		cout << "tmpSNP_num_valid: " << tmpSNP_num_valid << endl;
		cout << "tmpSNP_num_invalid: " << tmpSNP_num_invalid << endl;

		cout << "end of extract SNPs from formatted SNPfile !" << endl;
		cout << "start to insert formatted SNPs 2 chromStr" << endl;
		log_ofs << "end of extract SNPs from formatted SNPfile !" << endl;
		log_ofs << "start to insert formatted SNPs 2 chromStr" << endl;
		for(int tmp = 0; tmp < tmpSNP_num_valid; tmp++)
		{
			int tmpSNP_chrNameInt = SNP_chrNameIntVec[tmp];
			int tmpSNP_chrPos = SNP_chrPosVec[tmp];
			string tmpSNP_alterBase = SNP_alterBaseVec[tmp];
			
			if((SNP_refBaseVec[tmp] != (chromStr[tmpSNP_chrNameInt]).substr(tmpSNP_chrPos - 1, 1))
				&&(SNP_alterBaseVec[tmp] != (chromStr[tmpSNP_chrNameInt]).substr(tmpSNP_chrPos - 1, 1)))
			{
				cout << "error, SNP_refBaseVec[tmp] != (chromStr[tmpSNP_chrNameInt]).substr(tmpSNP_chrPos - 1, 1)" << endl;
				cout << "tmpSNP_chrNameInt: " << tmpSNP_chrNameInt << endl;
				cout << "tmpSNP_chrPos: " << tmpSNP_chrPos << endl;
				cout << "tmpSNP_refBase: " << SNP_refBaseVec[tmp] << endl;
				cout << "base in original CHROMOSOME: " << (chromStr[tmpSNP_chrNameInt]).substr(tmpSNP_chrPos - 1, 1) << endl;
				cout << "tmpSNP_alterBase: " << tmpSNP_alterBase << endl;
				//exit(1);
			}
			chromStr[tmpSNP_chrNameInt].replace(tmpSNP_chrPos - 1, 1, tmpSNP_alterBase);		
		}
		cout << "end of inserting formatted SNPs 2 chromStr" << endl;
		log_ofs << "end of inserting formatted SNPs 2 chromStr" << endl;
	}

	void insertSNP2chromStr(string& SNPfilePath)
	{
		cout << "start to extract SNPs from formatted SNPfile ! " << endl;
		vector<int> SNP_chrNameIntVec;
		vector<int> SNP_chrPosVec;
		vector<string> SNP_refBaseVec;
		vector<string> SNP_alterBaseVec;
		ifstream SNP_ifs(SNPfilePath.c_str());
		int tmpSNP_num_valid = 0;
		int tmpSNP_num_invalid = 0;
		while(!SNP_ifs.eof())
		{
			string tmpSNPstr;
			getline(SNP_ifs, tmpSNPstr);
			if(tmpSNPstr == "")
				break;
			int tabLoc_1 = tmpSNPstr.find("\t");
			int tabLoc_2 = tmpSNPstr.find("\t", tabLoc_1 + 1);
			int tabLoc_3 = tmpSNPstr.find("\t", tabLoc_2 + 1);			
			string tmpChrNameStr = tmpSNPstr.substr(0, tabLoc_1);
			int tmpChrNameInt = this->convertStringToInt(tmpChrNameStr);
			string tmpPosStr = tmpSNPstr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
			int tmpPos = atoi(tmpPosStr.c_str());			
			string tmpRefBase = tmpSNPstr.substr(tabLoc_2 + 1, 1);
			string tmpAlterBase = tmpSNPstr.substr(tabLoc_3 + 1, 1);
			if(tmpChrNameInt < 0)
				tmpSNP_num_invalid ++;
			else
			{
				tmpSNP_num_valid ++;
				SNP_chrNameIntVec.push_back(tmpChrNameInt);
				SNP_chrPosVec.push_back(tmpPos);
				SNP_refBaseVec.push_back(tmpRefBase);
				SNP_alterBaseVec.push_back(tmpAlterBase);
			}
		}
		SNP_ifs.close();
		cout << "tmpSNP_num_valid: " << tmpSNP_num_valid << endl;
		cout << "tmpSNP_num_invalid: " << tmpSNP_num_invalid << endl;

		cout << "end of extract SNPs from formatted SNPfile !" << endl;
		cout << "start to insert formatted SNPs 2 chromStr" << endl;
		for(int tmp = 0; tmp < tmpSNP_num_valid; tmp++)
		{
			int tmpSNP_chrNameInt = SNP_chrNameIntVec[tmp];
			int tmpSNP_chrPos = SNP_chrPosVec[tmp];
			string tmpSNP_alterBase = SNP_alterBaseVec[tmp];
			
			if((SNP_refBaseVec[tmp] != (chromStr[tmpSNP_chrNameInt]).substr(tmpSNP_chrPos - 1, 1))
				&&(SNP_alterBaseVec[tmp] != (chromStr[tmpSNP_chrNameInt]).substr(tmpSNP_chrPos - 1, 1)))
			{
				cout << "error, SNP_refBaseVec[tmp] != (chromStr[tmpSNP_chrNameInt]).substr(tmpSNP_chrPos - 1, 1)" << endl;
				cout << "tmpSNP_chrNameInt: " << tmpSNP_chrNameInt << endl;
				cout << "tmpSNP_chrPos: " << tmpSNP_chrPos << endl;
				cout << "tmpSNP_refBase: " << SNP_refBaseVec[tmp] << endl;
				cout << "base in original CHROMOSOME: " << (chromStr[tmpSNP_chrNameInt]).substr(tmpSNP_chrPos - 1, 1) << endl;
				cout << "tmpSNP_alterBase: " << tmpSNP_alterBase << endl;
				//exit(1);
			}
			chromStr[tmpSNP_chrNameInt].replace(tmpSNP_chrPos - 1, 1, tmpSNP_alterBase);		
		}
		cout << "end of inserting formatted SNPs 2 chromStr" << endl;
	}

	void convertSNPbaseBackToReferenceBase(string& SNPfilePath, ofstream& log_ofs)
	{
		cout << "start to extract SNPs from formatted SNPfile ! " << endl;
		log_ofs << "start to extract SNPs from formatted SNPfile ! " << endl;
		vector<int> SNP_chrNameIntVec;
		vector<int> SNP_chrPosVec;
		vector<string> SNP_refBaseVec;
		vector<string> SNP_alterBaseVec;
		ifstream SNP_ifs(SNPfilePath.c_str());
		int tmpSNP_num_valid = 0;
		int tmpSNP_num_invalid = 0;
		while(!SNP_ifs.eof())
		{
			string tmpSNPstr;
			getline(SNP_ifs, tmpSNPstr);
			if(tmpSNPstr == "")
				break;
			int tabLoc_1 = tmpSNPstr.find("\t");
			int tabLoc_2 = tmpSNPstr.find("\t", tabLoc_1 + 1);
			int tabLoc_3 = tmpSNPstr.find("\t", tabLoc_2 + 1);			
			string tmpChrNameStr = tmpSNPstr.substr(0, tabLoc_1);
			int tmpChrNameInt = this->convertStringToInt(tmpChrNameStr);
			string tmpPosStr = tmpSNPstr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
			int tmpPos = atoi(tmpPosStr.c_str());			
			string tmpRefBase = tmpSNPstr.substr(tabLoc_2 + 1, 1);
			string tmpAlterBase = tmpSNPstr.substr(tabLoc_3 + 1, 1);
			if(tmpChrNameInt < 0)
				tmpSNP_num_invalid ++;
			else
			{
				tmpSNP_num_valid ++;
				SNP_chrNameIntVec.push_back(tmpChrNameInt);
				SNP_chrPosVec.push_back(tmpPos);
				SNP_refBaseVec.push_back(tmpRefBase);
				SNP_alterBaseVec.push_back(tmpAlterBase);
			}
		}
		SNP_ifs.close();
		log_ofs << "tmpSNP_num_valid: " << tmpSNP_num_valid << endl;
		log_ofs << "tmpSNP_num_invalid: " << tmpSNP_num_invalid << endl;
		cout << "tmpSNP_num_valid: " << tmpSNP_num_valid << endl;
		cout << "tmpSNP_num_invalid: " << tmpSNP_num_invalid << endl;

		cout << "end of extract SNPs from formatted SNPfile !" << endl;
		cout << "start to convert SNPbases back to reference bases" << endl;
		log_ofs << "end of extract SNPs from formatted SNPfile !" << endl;
		log_ofs << "start to convert SNPbases back to reference bases" << endl;
		for(int tmp = 0; tmp < tmpSNP_num_valid; tmp++)
		{
			int tmpSNP_chrNameInt = SNP_chrNameIntVec[tmp];
			int tmpSNP_chrPos = SNP_chrPosVec[tmp];
			//string tmpSNP_alterBase = SNP_alterBaseVec[tmp];
			string tmpSNP_refBase = SNP_refBaseVec[tmp];
			// if((SNP_refBaseVec[tmp] != (chromStr[tmpSNP_chrNameInt]).substr(tmpSNP_chrPos - 1, 1))
			// 	&&(SNP_alterBaseVec[tmp] != (chromStr[tmpSNP_chrNameInt]).substr(tmpSNP_chrPos - 1, 1)))
			// {
			// 	cout << "error, SNP_refBaseVec[tmp] != (chromStr[tmpSNP_chrNameInt]).substr(tmpSNP_chrPos - 1, 1)" << endl;
			// 	cout << "tmpSNP_chrNameInt: " << tmpSNP_chrNameInt << endl;
			// 	cout << "tmpSNP_chrPos: " << tmpSNP_chrPos << endl;
			// 	cout << "tmpSNP_refBase: " << SNP_refBaseVec[tmp] << endl;
			// 	cout << "base in original CHROMOSOME: " << (chromStr[tmpSNP_chrNameInt]).substr(tmpSNP_chrPos - 1, 1) << endl;
			// 	cout << "tmpSNP_alterBase: " << tmpSNP_alterBase << endl;
			// 	//exit(1);
			// }
			chromStr[tmpSNP_chrNameInt].replace(tmpSNP_chrPos - 1, 1, tmpSNP_refBase);		
		}
		cout << "end of converting SNPbases back to reference bases" << endl;
		log_ofs << "end of converting SNPbases back to reference bases" << endl;
	}

	void insertSNP2chromStr_outputSNPmer(string& SNPfilePath, string& SNPmerFilePath, ofstream& log_ofs, int tmpSNPmerLength)
	{
		cout << "start to extract SNPs from formatted SNPfile ! " << endl;
		log_ofs << "start to extract SNPs from formatted SNPfile ! " << endl;
		// loading SNPs from SNPfile
		vector<int> SNP_chrNameIntVec;
		vector<int> SNP_chrPosVec;
		vector<string> SNP_refBaseVec;
		vector<string> SNP_alterBaseVec;
		ifstream SNP_ifs(SNPfilePath.c_str());
		int tmpSNP_num_valid = 0;
		int tmpSNP_num_invalid = 0;
		while(!SNP_ifs.eof())
		{
			string tmpSNPstr;
			getline(SNP_ifs, tmpSNPstr);
			if(tmpSNPstr == "")
				break;
			int tabLoc_1 = tmpSNPstr.find("\t");
			int tabLoc_2 = tmpSNPstr.find("\t", tabLoc_1 + 1);
			int tabLoc_3 = tmpSNPstr.find("\t", tabLoc_2 + 1);			
			string tmpChrNameStr = tmpSNPstr.substr(0, tabLoc_1);
			int tmpChrNameInt = this->convertStringToInt(tmpChrNameStr);
			string tmpPosStr = tmpSNPstr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
			int tmpPos = atoi(tmpPosStr.c_str());			
			string tmpRefBase = tmpSNPstr.substr(tabLoc_2 + 1, 1);
			string tmpAlterBase = tmpSNPstr.substr(tabLoc_3 + 1, 1);
			if(tmpChrNameInt < 0)
				tmpSNP_num_invalid ++;
			else
			{
				tmpSNP_num_valid ++;
				SNP_chrNameIntVec.push_back(tmpChrNameInt);
				SNP_chrPosVec.push_back(tmpPos);
				SNP_refBaseVec.push_back(tmpRefBase);
				SNP_alterBaseVec.push_back(tmpAlterBase);
			}
		}
		SNP_ifs.close();
		log_ofs << "tmpSNP_num_valid: " << tmpSNP_num_valid << endl;
		log_ofs << "tmpSNP_num_invalid: " << tmpSNP_num_invalid << endl;
		cout << "tmpSNP_num_valid: " << tmpSNP_num_valid << endl;
		cout << "tmpSNP_num_invalid: " << tmpSNP_num_invalid << endl;
		// replace refBase with altBase in SNP file
		cout << "end of extract SNPs from formatted SNPfile !" << endl;
		cout << "start to insert formatted SNPs 2 chromStr" << endl;
		log_ofs << "end of extract SNPs from formatted SNPfile !" << endl;
		log_ofs << "start to insert formatted SNPs 2 chromStr" << endl;
		for(int tmp = 0; tmp < tmpSNP_num_valid; tmp++)
		{
			int tmpSNP_chrNameInt = SNP_chrNameIntVec[tmp];
			int tmpSNP_chrPos = SNP_chrPosVec[tmp];
			string tmpSNP_alterBase = SNP_alterBaseVec[tmp];
			
			if((SNP_refBaseVec[tmp] != (chromStr[tmpSNP_chrNameInt]).substr(tmpSNP_chrPos - 1, 1))
				&&(SNP_alterBaseVec[tmp] != (chromStr[tmpSNP_chrNameInt]).substr(tmpSNP_chrPos - 1, 1)))
			{
				cout << "error, SNP_refBaseVec[tmp] != (chromStr[tmpSNP_chrNameInt]).substr(tmpSNP_chrPos - 1, 1)" << endl;
				cout << "tmpSNP_chrNameInt: " << tmpSNP_chrNameInt << endl;
				cout << "tmpSNP_chrPos: " << tmpSNP_chrPos << endl;
				cout << "tmpSNP_refBase: " << SNP_refBaseVec[tmp] << endl;
				cout << "base in original CHROMOSOME: " << (chromStr[tmpSNP_chrNameInt]).substr(tmpSNP_chrPos - 1, 1) << endl;
				cout << "tmpSNP_alterBase: " << tmpSNP_alterBase << endl;
				//exit(1);
			}
			chromStr[tmpSNP_chrNameInt].replace(tmpSNP_chrPos - 1, 1, tmpSNP_alterBase);		
		}
		cout << "end of inserting formatted SNPs 2 chromStr" << endl;
		log_ofs << "end of inserting formatted SNPs 2 chromStr" << endl;
		cout << "start to output SNPmers" << endl;
		log_ofs << "start to output SNPmers" << endl;		
		// output SNPmers
		int tmpSNPlocInSyntheticSNPseq = (tmpSNPmerLength - 1)/2;
		ofstream SNPmer_ofs(SNPmerFilePath.c_str());
		for(int tmp = 0; tmp < tmpSNP_num_valid; tmp++)
		{
			int tmpSNP_chrNameInt = SNP_chrNameIntVec[tmp];
			int tmpSNP_chrPos = SNP_chrPosVec[tmp];
			int tmpSNP_chr_length = this->returnChromLength(tmpSNP_chrNameInt);
			if(!((tmpSNP_chrPos - tmpSNPlocInSyntheticSNPseq < 2)||(tmpSNP_chrPos + tmpSNPlocInSyntheticSNPseq > tmpSNP_chr_length - 2)))
			// set as <2 and -2 just to be safe
			{
				SNPmer_ofs << ">" << this->returnChrNameStr(tmpSNP_chrNameInt) << ":" 
					<< tmpSNP_chrPos - tmpSNPlocInSyntheticSNPseq << ":" << tmpSNPmerLength << "M:" 
					<< tmpSNPlocInSyntheticSNPseq << ":" << tmpSNPmerLength << ":" << endl;
				SNPmer_ofs << this->returnChromStrSubstr(tmpSNP_chrNameInt, tmpSNP_chrPos - tmpSNPlocInSyntheticSNPseq, tmpSNPmerLength) << endl;
			}
		}
		SNPmer_ofs.close();
		cout << "end of output SNPmers" << endl;
		log_ofs << "end of output SNPmers" << endl;	
	}
};

int extensionBackwards_errorTolerated(int chrNameInt, int chrMapPos, Index_Info* indexInfo, 
	const string& readSeqInProcess, int extensionStartLocInRead, int extensionEndLocInRead,
	vector<int>& newMismatchPosVec_head, vector<char>& newMismatchCharVec_head, 
	int matchBaseNum_min_perMismatch, int matchBaseNum_min_twoMismatch)  // used to extend to head
{
	//int matchBaseNum_min_perMismatch = 5;
	//int matchBaseNum_min_twoMismatch = 8;
	//cout << "extensionBackwards_errorTolerated starts ..." << endl;

	int maxExtensionLength = extensionStartLocInRead - extensionEndLocInRead + 1; // extensionStartLocInRead > extensionEndLocInRead

	int tmp_loc_extend = 0;
	while(1)
	{
		int tmpLocInChr = chrMapPos - tmp_loc_extend;
		int tmpLocInRead = extensionStartLocInRead - tmp_loc_extend;
		if(tmpLocInRead < extensionEndLocInRead)
			return tmp_loc_extend;
		if(indexInfo->getCharInChromosome(chrNameInt, tmpLocInChr) != readSeqInProcess.at(tmpLocInRead-1)) // mismatch found!
		{
			bool checkNextBasesToSupportMismatchBase_bool = true;
			if(tmp_loc_extend + matchBaseNum_min_perMismatch + 1 > maxExtensionLength)
				return tmp_loc_extend;			
			int newMismatchPosInRead = 0;
			int newMismatchPosInChr = 0;	
			int offset_mismatch = 0;		
			for(int tmp = 1; tmp <= matchBaseNum_min_perMismatch; tmp++)
			{
				if(indexInfo->getCharInChromosome(chrNameInt, tmpLocInChr - tmp) != readSeqInProcess.at(tmpLocInRead - tmp - 1))
				{
					checkNextBasesToSupportMismatchBase_bool = false;
					newMismatchPosInRead = tmpLocInRead - tmp;
					newMismatchPosInChr = tmpLocInChr - tmp;
					offset_mismatch = tmp;
					break;
				}
			}
			if(checkNextBasesToSupportMismatchBase_bool)
			{
				tmp_loc_extend = tmp_loc_extend + 1 + matchBaseNum_min_perMismatch;
				//if(STORE_MISMATCH_POS)
				//{
					newMismatchPosVec_head.push_back(tmpLocInRead);
					//if(STORE_MISMATCH_CHA)
					//{
						newMismatchCharVec_head.push_back(readSeqInProcess.at(tmpLocInRead-1));
					//}
				//}
			}
			else
			{
				if(newMismatchPosInRead - matchBaseNum_min_twoMismatch < extensionEndLocInRead)
				{
					return tmp_loc_extend;
				}
				for(int tmp = 1; tmp <= matchBaseNum_min_twoMismatch; tmp++)
				{
					if(indexInfo->getCharInChromosome(chrNameInt, newMismatchPosInChr - tmp) 
						!= readSeqInProcess.at(newMismatchPosInRead - tmp - 1))
					{
						//checkNextBasesToSupportMismatchBase_bool = false;
						return tmp_loc_extend;
					}						
				}
				//if(STORE_MISMATCH_POS)
				//{		
					newMismatchPosVec_head.push_back(tmpLocInRead);//Xinan: fixed in 11/11/2014
					newMismatchPosVec_head.push_back(newMismatchPosInRead);

					//if(STORE_MISMATCH_CHA)
					//{
						newMismatchCharVec_head.push_back(readSeqInProcess.at(tmpLocInRead-1));//Xinan: fixed in 11/11/2014
						newMismatchCharVec_head.push_back(readSeqInProcess.at(newMismatchPosInRead-1));
					//}
				//}
				
				tmp_loc_extend = tmp_loc_extend + matchBaseNum_min_twoMismatch + offset_mismatch;
			}						
		}
		else
		{
			tmp_loc_extend ++;
		}
	}
	//return currentBestExtensionLength;
}

int extensionBackwards_errorTolerated_finalStepForAligner(int chrNameInt, int chrMapPos, Index_Info* indexInfo, 
	const string& readSeqInProcess, int extensionStartLocInRead, int extensionEndLocInRead,
	vector<int>& newMismatchPosVec_head, vector<char>& newMismatchCharVec_head, 
	int matchBaseNum_min_perMismatch, int matchBaseNum_min_twoMismatch)  // used to extend to head
{
	//int matchBaseNum_min_perMismatch = 5;
	//int matchBaseNum_min_twoMismatch = 8;
	//cout << "extensionBackwards_errorTolerated starts ..." << endl;

	int maxExtensionLength = extensionStartLocInRead - extensionEndLocInRead + 1; // extensionStartLocInRead > extensionEndLocInRead
	//cout << "chrNameInt: " << chrNameInt << endl;
	//cout << "chrMapPos: " << chrMapPos << endl;
	//cout << "maxExtensionLength: " << maxExtensionLength << endl;
	int tmp_loc_extend = 0;
	while(1)
	{
		int tmpLocInChr = chrMapPos - tmp_loc_extend;
		int tmpLocInRead = extensionStartLocInRead - tmp_loc_extend;
		if(tmpLocInRead < extensionEndLocInRead)
			return tmp_loc_extend;
		if(indexInfo->getCharInChromosome(chrNameInt, tmpLocInChr) != readSeqInProcess.at(tmpLocInRead-1)) // mismatch found!
		{
			//cout << "mismatch found in seqBase: " << tmpLocInRead << endl; 
			bool checkNextBasesToSupportMismatchBase_bool = true;
			if(tmp_loc_extend + matchBaseNum_min_perMismatch + 1 > maxExtensionLength)
			{
				// do not require full matchBaseNum_min_perMismatch bases to support one mismatch for the final step of aligner
				//cout << "no enought base to support the mismatch " << endl;
				//cout << "bases to check: " << maxExtensionLength - tmp_loc_extend - 1 << endl;
				for(int tmp = 1; tmp <= maxExtensionLength - tmp_loc_extend - 1; tmp++)
				{
					if(indexInfo->getCharInChromosome(chrNameInt, tmpLocInChr - tmp) != readSeqInProcess.at(tmpLocInRead - tmp - 1))					
					{
						//cout << "found new mismatch at: " << tmpLocInRead - tmp << endl;
						return tmp_loc_extend;
					}
				}
				newMismatchPosVec_head.push_back(tmpLocInRead);
				newMismatchCharVec_head.push_back(readSeqInProcess.at(tmpLocInRead-1));
				return maxExtensionLength;			
			}
			int newMismatchPosInRead = 0;
			int newMismatchPosInChr = 0;	
			int offset_mismatch = 0;		
			for(int tmp = 1; tmp <= matchBaseNum_min_perMismatch; tmp++)
			{
				if(indexInfo->getCharInChromosome(chrNameInt, tmpLocInChr - tmp) != readSeqInProcess.at(tmpLocInRead - tmp - 1))
				{
					checkNextBasesToSupportMismatchBase_bool = false;
					newMismatchPosInRead = tmpLocInRead - tmp;
					newMismatchPosInChr = tmpLocInChr - tmp;
					offset_mismatch = tmp;
					break;
				}
			}
			if(checkNextBasesToSupportMismatchBase_bool)
			{
				tmp_loc_extend = tmp_loc_extend + 1 + matchBaseNum_min_perMismatch;

				newMismatchPosVec_head.push_back(tmpLocInRead);
				newMismatchCharVec_head.push_back(readSeqInProcess.at(tmpLocInRead-1));
			}
			else
			{
				if(newMismatchPosInRead - matchBaseNum_min_twoMismatch < extensionEndLocInRead)
				{
					// do not require full matchBaseNum_min_twoMismatch bases to support two mismatch for the final step of aligner
					for(int tmp = 1; tmp <= newMismatchPosInRead - extensionEndLocInRead; tmp++)
					{
						if(indexInfo->getCharInChromosome(chrNameInt, newMismatchPosInChr - tmp) 
							!= readSeqInProcess.at(newMismatchPosInRead - tmp - 1))					
							return tmp_loc_extend;
					}		
					newMismatchPosVec_head.push_back(tmpLocInRead);
					newMismatchPosVec_head.push_back(newMismatchPosInRead);

					newMismatchCharVec_head.push_back(readSeqInProcess.at(tmpLocInRead-1));
					newMismatchCharVec_head.push_back(readSeqInProcess.at(newMismatchPosInRead-1));
					return maxExtensionLength;
				}
				for(int tmp = 1; tmp <= matchBaseNum_min_twoMismatch; tmp++)
				{
					if(indexInfo->getCharInChromosome(chrNameInt, newMismatchPosInChr - tmp) 
						!= readSeqInProcess.at(newMismatchPosInRead - tmp - 1))
					{
						//checkNextBasesToSupportMismatchBase_bool = false;
						return tmp_loc_extend;
					}						
				}

				newMismatchPosVec_head.push_back(tmpLocInRead);//Xinan: fixed in 11/11/2014
				newMismatchPosVec_head.push_back(newMismatchPosInRead);

				newMismatchCharVec_head.push_back(readSeqInProcess.at(tmpLocInRead-1));//Xinan: fixed in 11/11/2014
				newMismatchCharVec_head.push_back(readSeqInProcess.at(newMismatchPosInRead-1));

				tmp_loc_extend = tmp_loc_extend + matchBaseNum_min_twoMismatch + offset_mismatch;
			}						
		}
		else
		{
			tmp_loc_extend ++;
		}
	}
	//return currentBestExtensionLength;
}

int extensionForwards_errorTolerated(int chrNameInt, int chrMapPos, Index_Info* indexInfo,
	const string& readSeqInProcess, int extensionStartLocInRead, int extensionEndLocInRead,
	vector<int>& newMismatchPosVec_tail, vector<char>& newMismatchCharVec_tail, 
	int matchBaseNum_min_perMismatch, int matchBaseNum_min_twoMismatch)  // used to extend to tail
{
	//cout << "start to fix extensionForwards_errorTolerated ..." << endl;
	//cout << "readSeqInProcess: " << readSeqInProcess << endl;

	int maxExtensionLength = extensionEndLocInRead - extensionStartLocInRead + 1;
	
	//cout << "maxExtensionLength: " << maxExtensionLength << endl;
	int tmp_loc_extend = 0;
	//while(1)
	//{
		while(1) // looking for mismatchPos
		{
			int tmpLocInChr = chrMapPos + tmp_loc_extend;
			int tmpLocInRead = extensionStartLocInRead + tmp_loc_extend;
			if(tmpLocInRead > extensionEndLocInRead)
				return tmp_loc_extend; 
			//cout << "start to check match or mismatch " << endl;
			if(indexInfo->getCharInChromosome(chrNameInt, tmpLocInChr) != readSeqInProcess.at(tmpLocInRead-1)) // mismatch found!
			{
				//cout << "found mismatch ..." << endl;
				bool checkNextBasesToSupportMismatchBase_bool = true;
				//cout << "tmp_loc_extend + matchBaseNum_min_perMismatch: " << tmp_loc_extend + matchBaseNum_min_perMismatch << endl;
				if(tmp_loc_extend + matchBaseNum_min_perMismatch + 1 > maxExtensionLength)
				{
					//cout << "tmp_loc_extend: " << tmp_loc_extend << endl;
					return tmp_loc_extend;
				}
				
				//cout << "start to check match bases to support " << endl;
				int newMismatchPosInRead = 0;
				int newMismatchPosInChr = 0;
				int offset_mismatch = 0;
				for(int tmp = 1; tmp <= matchBaseNum_min_perMismatch; tmp++)
				{
					if(indexInfo->getCharInChromosome(chrNameInt, tmpLocInChr + tmp) != readSeqInProcess.at(tmpLocInRead + tmp-1))
					{
						checkNextBasesToSupportMismatchBase_bool = false;
						newMismatchPosInRead = tmpLocInRead + tmp;
						newMismatchPosInChr = tmpLocInChr + tmp;
						offset_mismatch = tmp;
						break;
					}
				}
				//cout << "checkNextBasesToSupportMismatchBase_bool: " << checkNextBasesToSupportMismatchBase_bool << endl;
				if(checkNextBasesToSupportMismatchBase_bool)
				{
					tmp_loc_extend = tmp_loc_extend + 1 + matchBaseNum_min_perMismatch;
					//if(STORE_MISMATCH_POS)
					//{
						newMismatchPosVec_tail.push_back(tmpLocInRead);
						//if(STORE_MISMATCH_CHA)
						//{
							newMismatchCharVec_tail.push_back(readSeqInProcess.at(tmpLocInRead-1));
						//}
					//}
				}
				else
				{
					if(newMismatchPosInRead + matchBaseNum_min_twoMismatch > extensionEndLocInRead)
					{
						return tmp_loc_extend;
					}
					for(int tmp = 1; tmp <= matchBaseNum_min_twoMismatch; tmp++)
					{
						if(indexInfo->getCharInChromosome(chrNameInt, newMismatchPosInChr + tmp) 
							!= readSeqInProcess.at(newMismatchPosInRead + tmp - 1))
						{
							//checkNextBasesToSupportMismatchBase_bool = false;
							return tmp_loc_extend;
							//break;
						}						
					}
					//if(STORE_MISMATCH_POS)
					//{		
						newMismatchPosVec_tail.push_back(tmpLocInRead);
						newMismatchPosVec_tail.push_back(newMismatchPosInRead);
						//if(STORE_MISMATCH_CHA)
						//{
							newMismatchCharVec_tail.push_back(readSeqInProcess.at(tmpLocInRead-1));
							newMismatchCharVec_tail.push_back(readSeqInProcess.at(newMismatchPosInRead-1));
						//}
					//}					
					tmp_loc_extend = tmp_loc_extend + matchBaseNum_min_twoMismatch + offset_mismatch;
				}
			}
			else
			{
				tmp_loc_extend ++;
			}
		}
	//}	
	//return currentBestExtensionLength;
}

int extensionForwards_errorTolerated_finalStepForAligner(int chrNameInt, int chrMapPos, Index_Info* indexInfo,
	const string& readSeqInProcess, int extensionStartLocInRead, int extensionEndLocInRead,
	vector<int>& newMismatchPosVec_tail, vector<char>& newMismatchCharVec_tail, 
	int matchBaseNum_min_perMismatch, int matchBaseNum_min_twoMismatch)  // used to extend to tail
{
	//cout << "start to fix extensionForwards_errorTolerated ..." << endl;
	//cout << "readSeqInProcess: " << readSeqInProcess << endl;

	int maxExtensionLength = extensionEndLocInRead - extensionStartLocInRead + 1;
	//cout << "chrNameInt: " << chrNameInt << endl;
	//cout << "chrMapPos: " << chrMapPos << endl;
	//cout << "maxExtensionLength: " << maxExtensionLength << endl;	
	//cout << "maxExtensionLength: " << maxExtensionLength << endl;
	int tmp_loc_extend = 0;
	//while(1)
	//{
		while(1) // looking for mismatchPos
		{
			int tmpLocInChr = chrMapPos + tmp_loc_extend;
			int tmpLocInRead = extensionStartLocInRead + tmp_loc_extend;
			if(tmpLocInRead > extensionEndLocInRead)
				return tmp_loc_extend; 
			//cout << "start to check match or mismatch " << endl;
			if(indexInfo->getCharInChromosome(chrNameInt, tmpLocInChr) != readSeqInProcess.at(tmpLocInRead-1)) // mismatch found!
			{
				//cout << "mismatch found in seqBase: " << tmpLocInRead << endl; 
				//cout << "found mismatch ..." << endl;
				bool checkNextBasesToSupportMismatchBase_bool = true;
				//cout << "tmp_loc_extend + matchBaseNum_min_perMismatch: " << tmp_loc_extend + matchBaseNum_min_perMismatch << endl;
				if(tmp_loc_extend + matchBaseNum_min_perMismatch + 1 > maxExtensionLength)
				{
					// do not require full matchBaseNum_min_perMismatch bases to support one mismatch for the final step of aligner
					//cout << "no enough bases to support the mismatch " << endl;
					//cout << "tmp_loc_extend: " << tmp_loc_extend << endl;
					//cout << "bases to check: " << maxExtensionLength - tmp_loc_extend - 1 << endl;
					for(int tmp = 1; tmp <= maxExtensionLength - tmp_loc_extend - 1; tmp++)
					{	
						if(indexInfo->getCharInChromosome(chrNameInt, tmpLocInChr + tmp) 
							!= readSeqInProcess.at(tmpLocInRead + tmp - 1))
						{
							//cout << "found new mismatch in: " << tmpLocInRead + tmp << endl;
							return tmp_loc_extend;
						}
					}
					newMismatchPosVec_tail.push_back(tmpLocInRead);
					newMismatchCharVec_tail.push_back(readSeqInProcess.at(tmpLocInRead-1));
					return maxExtensionLength;
				}
				
				//cout << "start to check match bases to support " << endl;
				int newMismatchPosInRead = 0;
				int newMismatchPosInChr = 0;
				int offset_mismatch = 0;
				for(int tmp = 1; tmp <= matchBaseNum_min_perMismatch; tmp++)
				{
					if(indexInfo->getCharInChromosome(chrNameInt, tmpLocInChr + tmp) != readSeqInProcess.at(tmpLocInRead + tmp-1))
					{
						//cout << "found the 2nd mismatch at base: " << tmpLocInRead + tmp << endl;
						checkNextBasesToSupportMismatchBase_bool = false;
						newMismatchPosInRead = tmpLocInRead + tmp;
						newMismatchPosInChr = tmpLocInChr + tmp;
						offset_mismatch = tmp;
						break;
					}
				}
				//cout << "checkNextBasesToSupportMismatchBase_bool: " << checkNextBasesToSupportMismatchBase_bool << endl;
				if(checkNextBasesToSupportMismatchBase_bool)
				{
					tmp_loc_extend = tmp_loc_extend + 1 + matchBaseNum_min_perMismatch;
					newMismatchPosVec_tail.push_back(tmpLocInRead);
					newMismatchCharVec_tail.push_back(readSeqInProcess.at(tmpLocInRead-1));
				}
				else
				{
					if(newMismatchPosInRead + matchBaseNum_min_twoMismatch > extensionEndLocInRead)
					{
						//cout << "no enough bases to support the two mismatches " << endl;
						//cout << "bases to check: " << extensionEndLocInRead - newMismatchPosInRead << endl;
						// do not require full matchBaseNum_min_twoMismatch bases to support two mismatch for the final step of aligner
						for(int tmp = 1; tmp <= extensionEndLocInRead - newMismatchPosInRead; tmp++)
						{
							if(indexInfo->getCharInChromosome(chrNameInt, newMismatchPosInChr + tmp) 
								!= readSeqInProcess.at(newMismatchPosInRead + tmp - 1))		
								return tmp_loc_extend;
						}
						newMismatchPosVec_tail.push_back(tmpLocInRead);
						newMismatchPosVec_tail.push_back(newMismatchPosInRead);
						newMismatchCharVec_tail.push_back(readSeqInProcess.at(tmpLocInRead-1));
						newMismatchCharVec_tail.push_back(readSeqInProcess.at(newMismatchPosInRead-1));
						return maxExtensionLength;
					}
					for(int tmp = 1; tmp <= matchBaseNum_min_twoMismatch; tmp++)
					{
						if(indexInfo->getCharInChromosome(chrNameInt, newMismatchPosInChr + tmp) 
							!= readSeqInProcess.at(newMismatchPosInRead + tmp - 1))
						{
							//checkNextBasesToSupportMismatchBase_bool = false;
							return tmp_loc_extend;
							//break;
						}						
					}
	
					newMismatchPosVec_tail.push_back(tmpLocInRead);
					newMismatchPosVec_tail.push_back(newMismatchPosInRead);
					newMismatchCharVec_tail.push_back(readSeqInProcess.at(tmpLocInRead-1));
					newMismatchCharVec_tail.push_back(readSeqInProcess.at(newMismatchPosInRead-1));
			
					tmp_loc_extend = tmp_loc_extend + matchBaseNum_min_twoMismatch + offset_mismatch;
				}
			}
			else
			{
				tmp_loc_extend ++;
			}
		}
	//}	
	//return currentBestExtensionLength;
}

int extensionForwards_errorTolerated_finalStepForAligner_wholeGenome_duringSegmengMapping(
	//int chrNameInt, int chrMapPos, 
	unsigned int wholeGenomeMapPos,
	Index_Info* indexInfo,
	//const string& readSeqInProcess, 
	char* read,
	int extensionStartLocInRead, int extensionEndLocInRead,
	vector<int>& newMismatchPosVec_tail, vector<char>& newMismatchCharVec_tail, 
	int matchBaseNum_min_perMismatch, int matchBaseNum_min_twoMismatch)  // used to extend to tail
{
	//cout << "start to fix extensionForwards_errorTolerated ..." << endl;
	//cout << "readSeqInProcess: " << readSeqInProcess << endl;

	int maxExtensionLength = extensionEndLocInRead - extensionStartLocInRead + 1;
	//cout << "chrNameInt: " << chrNameInt << endl;
	//cout << "chrMapPos: " << chrMapPos << endl;
	//cout << "maxExtensionLength: " << maxExtensionLength << endl;	
	//cout << "maxExtensionLength: " << maxExtensionLength << endl;
	int tmp_loc_extend = 0;
	//while(1)
	//{
		while(1) // looking for mismatchPos
		{
			//int tmpLocInChr = chrMapPos + tmp_loc_extend;
			unsigned int tmpLocInWholeGenome = wholeGenomeMapPos + tmp_loc_extend;
			int tmpLocInRead = extensionStartLocInRead + tmp_loc_extend;
			if(tmpLocInRead > extensionEndLocInRead)
				return tmp_loc_extend; 
			//cout << "start to check match or mismatch " << endl;
			if(//indexInfo->getCharInChromosome(chrNameInt, tmpLocInChr) 
				indexInfo->getCharInWholeGenome(tmpLocInWholeGenome)
				!= 
				//readSeqInProcess.at(tmpLocInRead-1)
				*(read + tmpLocInRead-1)
				) // mismatch found!
			{
				//cout << "mismatch found in seqBase: " << tmpLocInRead << endl; 
				//cout << "found mismatch ..." << endl;
				bool checkNextBasesToSupportMismatchBase_bool = true;
				//cout << "tmp_loc_extend + matchBaseNum_min_perMismatch: " << tmp_loc_extend + matchBaseNum_min_perMismatch << endl;
				if(tmp_loc_extend + matchBaseNum_min_perMismatch + 1 > maxExtensionLength)
				{
					// do not require full matchBaseNum_min_perMismatch bases to support one mismatch for the final step of aligner
					//cout << "no enough bases to support the mismatch " << endl;
					//cout << "tmp_loc_extend: " << tmp_loc_extend << endl;
					//cout << "bases to check: " << maxExtensionLength - tmp_loc_extend - 1 << endl;
					for(int tmp = 1; tmp <= maxExtensionLength - tmp_loc_extend - 1; tmp++)
					{	
						if(//indexInfo->getCharInChromosome(chrNameInt, tmpLocInChr + tmp) 
							indexInfo->getCharInWholeGenome(tmpLocInWholeGenome+tmp)
							!= *(read + tmpLocInRead + tmp - 1))
						{
							//cout << "found new mismatch in: " << tmpLocInRead + tmp << endl;
							return tmp_loc_extend;
						}
					}
					newMismatchPosVec_tail.push_back(tmpLocInRead);
					newMismatchCharVec_tail.push_back(*(read + tmpLocInRead-1));
					return maxExtensionLength;
				}
				
				//cout << "start to check match bases to support " << endl;
				int newMismatchPosInRead = 0;
				//int newMismatchPosInChr = 0;
				unsigned int newMismatchPosInWholeGenome = 0;
				int offset_mismatch = 0;
				for(int tmp = 1; tmp <= matchBaseNum_min_perMismatch; tmp++)
				{
					if(//indexInfo->getCharInChromosome(chrNameInt, tmpLocInChr + tmp) 
						indexInfo->getCharInWholeGenome(tmpLocInWholeGenome + tmp)
						!= *(read + tmpLocInRead + tmp-1))
					{
						//cout << "found the 2nd mismatch at base: " << tmpLocInRead + tmp << endl;
						checkNextBasesToSupportMismatchBase_bool = false;
						newMismatchPosInRead = tmpLocInRead + tmp;
						//newMismatchPosInChr = tmpLocInChr + tmp;
						newMismatchPosInWholeGenome = tmpLocInWholeGenome + tmp;
						offset_mismatch = tmp;
						break;
					}
				}
				//cout << "checkNextBasesToSupportMismatchBase_bool: " << checkNextBasesToSupportMismatchBase_bool << endl;
				if(checkNextBasesToSupportMismatchBase_bool)
				{
					tmp_loc_extend = tmp_loc_extend + 1 + matchBaseNum_min_perMismatch;
					newMismatchPosVec_tail.push_back(tmpLocInRead);
					newMismatchCharVec_tail.push_back(*(read + tmpLocInRead-1));
				}
				else
				{
					if(newMismatchPosInRead + matchBaseNum_min_twoMismatch > extensionEndLocInRead)
					{
						//cout << "no enough bases to support the two mismatches " << endl;
						//cout << "bases to check: " << extensionEndLocInRead - newMismatchPosInRead << endl;
						// do not require full matchBaseNum_min_twoMismatch bases to support two mismatch for the final step of aligner
						for(int tmp = 1; tmp <= extensionEndLocInRead - newMismatchPosInRead; tmp++)
						{
							if(//indexInfo->getCharInChromosome(chrNameInt, newMismatchPosInChr + tmp) 
								indexInfo->getCharInWholeGenome(newMismatchPosInWholeGenome + tmp)
								!= *(read + newMismatchPosInRead + tmp - 1))		
								return tmp_loc_extend;
						}
						newMismatchPosVec_tail.push_back(tmpLocInRead);
						newMismatchPosVec_tail.push_back(newMismatchPosInRead);
						newMismatchCharVec_tail.push_back(*(read + tmpLocInRead-1));
						newMismatchCharVec_tail.push_back(*(read + newMismatchPosInRead-1));
						return maxExtensionLength;
					}
					for(int tmp = 1; tmp <= matchBaseNum_min_twoMismatch; tmp++)
					{
						if(//indexInfo->getCharInChromosome(chrNameInt, newMismatchPosInChr + tmp) 
							indexInfo->getCharInWholeGenome(newMismatchPosInWholeGenome + tmp)
							!= *(read + newMismatchPosInRead + tmp - 1))
						{
							//checkNextBasesToSupportMismatchBase_bool = false;
							return tmp_loc_extend;
							//break;
						}						
					}
	
					newMismatchPosVec_tail.push_back(tmpLocInRead);
					newMismatchPosVec_tail.push_back(newMismatchPosInRead);
					newMismatchCharVec_tail.push_back(*(read + tmpLocInRead-1));
					newMismatchCharVec_tail.push_back(*(read + newMismatchPosInRead-1));
			
					tmp_loc_extend = tmp_loc_extend + matchBaseNum_min_twoMismatch + offset_mismatch;
				}
			}
			else
			{
				tmp_loc_extend ++;
			}
		}
	//}	
	//return currentBestExtensionLength;
}
#endif