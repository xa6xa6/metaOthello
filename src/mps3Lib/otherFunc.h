#ifndef OTHERFUNC_H
#define OTHERFUNC_H

//#include "index_info.h"
#include <string>
#include <string.h>

using namespace std;

int selectTheSmallestAmong4values(
	int v1, int v2, int v3, int v4)
{
	int tmpSmallest = v1;
	if(v2 < tmpSmallest)
		tmpSmallest = v2;
	if(v3 < tmpSmallest)
		tmpSmallest = v3;
	if(v4 < tmpSmallest)
		tmpSmallest = v4;
	return tmpSmallest;
}

int selectTheLargestAmong4values(
	int v1, int v2, int v3, int v4)
{
	int tmpLargest = v1;
	if(v2 > tmpLargest)
		tmpLargest = v2;
	if(v3 > tmpLargest)
		tmpLargest = v3;
	if(v4 > tmpLargest)
		tmpLargest = v4;
	return tmpLargest;		
}

char reverseComplement(char ch)
{
	if(ch == 'A')
		return 'T';
	else if(ch == 'C')
		return 'G';
	else if(ch == 'G')
		return 'C';
	else if(ch == 'T')
		return 'A';
	else if(ch == 'N')
	{
		return 'N';
	}
	else
	{
		cout << "incorrect Ori_Char in reverseComplement" << endl;
		exit(1);
		return 'X';
	}
}

string reverseComplementStr(char ch)
{
	if(ch == 'A')
		return "T";
	else if(ch == 'C')
		return "G";
	else if(ch == 'G')
		return "C";
	else if(ch == 'T')
		return "A";
	else if(ch == 'N')
		return "N";
	else
	{
		cout << "incorrect Ori_Char in reverseComplementStr" << endl;
		exit(1);
		return "X";		
	}
}

string convertCharArrayToReverseCompletmentStr(char* readChar, int readLength)
{
	string rcmReadStr = "";
	for(int tmp = 0; tmp < readLength; tmp++)
	{
		rcmReadStr += reverseComplementStr(*(readChar + readLength - 1 - tmp));
	}
	return rcmReadStr;
}

string covertCharToReverseComplement(const string& Ori_Char)
{
	if(Ori_Char == "A")
	{
		return "T";
	}
	else if(Ori_Char == "T")
	{
		return "A";
	}
	else if(Ori_Char == "G")
	{
		return "C";
	}
	else if(Ori_Char == "C")
	{
		return "G";
	}
	else if(Ori_Char == "N")
	{
		return "N";
	}
	else
	{
		cout << "incorrect Ori_Char in covertCharToReverseComplement" << endl;
		exit(1);
		return "X";
	}
}

string covertStringToReverseComplement(const string& originalString)
{
	int stringLength = originalString.size();
	string resultString = covertCharToReverseComplement(originalString.substr(stringLength-1, 1));
	for (int tmp = 1; tmp < stringLength; tmp++)
	{
		resultString = resultString + covertCharToReverseComplement(
			originalString.substr(stringLength-1-tmp, 1));
	}
	return resultString;
}

string convertStringToReverseComplement(const string& originalString)
{
	int stringLength = originalString.size();
	string resultString = covertCharToReverseComplement(originalString.substr(stringLength-1, 1));
	for (int tmp = 1; tmp < stringLength; tmp++)
	{
		resultString = resultString + covertCharToReverseComplement(
			originalString.substr(stringLength-1-tmp, 1));
	}
	return resultString;
}

string convertQualityScoreString2Reverse(const string& originalQualityScoreString)
{
	int stringLength = originalQualityScoreString.size();
	string resultString = originalQualityScoreString.substr(stringLength-1, 1);//covertCharToReverseComplement(originalString.substr(stringLength-1, 1));
	for (int tmp = 1; tmp < stringLength; tmp++)
	{
		resultString = resultString + originalQualityScoreString.substr(stringLength-1-tmp, 1);
			//covertCharToReverseComplement(originalString.substr(stringLength-1-tmp, 1));
	}
	return resultString;
}

int extendBackInChromSeq(int readLoc, const string& readSeq, int chromLoc, const string& chromSeq, int extendBackLengthMax)
{
	int tmp = 1;
	for (tmp = 1; tmp <= extendBackLengthMax; tmp++)
	{
		if(readSeq.at(readLoc - tmp - 1) != chromSeq.at(chromLoc - tmp - 1))
		{
			return tmp - 1;
		}
	}
	return extendBackLengthMax;
}


int skipOneMismatchAndExtendRead_backward(char* read, char* genome, 
	int startLocInRead, unsigned int startPosInGenome, int readLength)
{
	int tmpMismatchNum = 0;
	for(int tmp = startLocInRead; tmp <= readLength; tmp++)
	{	
		char tmpReadChar = *(read + startLocInRead - 1 + tmp);
		char tmpGenomeChar = *(genome + startPosInGenome - 1 + tmp);
		if(tmpReadChar != tmpGenomeChar)
		{
			tmpMismatchNum ++;
			if(tmpMismatchNum > 1)
			{
				return tmp - 1;
			}
		}
	}
	return readLength;
}

bool skipOneMismatchAndExtendBackwards(char *read_start, char* genome_start, 
	int startLoc, int endLoc, int& stopLoc)
{
	int tmpMismatchNum = 0;
	for(int tmp = startLoc; tmp <= endLoc; tmp++)
	{
		char tmpReadChar = *(read_start + tmp);
		char tmpGenomeChar = *(genome_start + tmp);
		if(tmpReadChar != tmpGenomeChar)
		{
			tmpMismatchNum ++;
			if(tmpMismatchNum > 1)
			{
				stopLoc = tmp;
				return false;
			}
		}
	}
	stopLoc = endLoc;
	return true;
}

int checkIH(string& tmpSamStr)
{
	int IH_startLoc = tmpSamStr.find("IH:i:");
	int HI_startLoc = tmpSamStr.find("HI:i:");
	int IH_int_startLoc = IH_startLoc + 5;
	int IH_int_endLoc = HI_startLoc - 2;
	string IH_int_str = tmpSamStr.substr(IH_int_startLoc, IH_int_endLoc - IH_int_startLoc + 1);
	//cout << "IH_int_str: " << IH_int_str << endl;
	int IH_int = atoi(IH_int_str.c_str());
	return IH_int;
}

bool forwardOrReverseComplement(int tmpFlag)
{
	if(tmpFlag & 0x16)
		return false;
	else
		return true;
}


bool mappedOrNot(int tmpFlag)
{
	if(tmpFlag & 0x4)
		return false;
	else
		return true;
}

bool primaryOrNot(int tmpFlag)
{
	if(tmpFlag & 0x100)
		return false;
	else
		return true;
}


inline char getCharRevComp(char ch)
{
	int chInt = ch - 'A';
	static const char alphatChar[26] = {'T', 'N', 'G', 'N', 'N', 'N', 'C',
		'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'N', 'N', 'N', 'N', 'N', 'A',
		'N', 'N', 'N', 'N', 'N', 'N'};
	return alphatChar[chInt];
}

string getRcmSeq(const string& readSeq)
{
	int readSeqLength = readSeq.length();

	char readRcmSeqChar[readSeqLength];

	readRcmSeqChar[0] = getCharRevComp(readSeq.at(readSeqLength-1));

	for(int tmp = 1; tmp < readSeqLength; tmp ++)
	{
		readRcmSeqChar[tmp] = getCharRevComp((readSeq.at(readSeqLength - tmp - 1)));
	}
	string rcmSeq = readRcmSeqChar;
	return rcmSeq.substr(0, readSeqLength);
}

string getRevSeq(const string& qualitySeq)
{
	int qualitySeqLength = qualitySeq.length();
	char qualityRevSeqChar[qualitySeqLength];
	qualityRevSeqChar[0] = qualitySeq.at(qualitySeqLength-1);
	for(int tmp = 1; tmp < qualitySeqLength; tmp++)
	{
		qualityRevSeqChar[tmp] = qualitySeq.at(qualitySeqLength - tmp -1);
	}
	string revSeq = qualityRevSeqChar;
	return revSeq.substr(0, qualitySeqLength);
}


class OtherFunc
{
public:
	vector<char> char2CharRcmVec;
	//vector<string> str2StrRcmVec;


	OtherFunc()
	{
		for(int tmp = 'A'-'A'; tmp <= 'Z' - 'A'; tmp++)
		{
			char2CharRcmVec.push_back('N');
		}
		char2CharRcmVec['A'-'A'] = 'T';
		char2CharRcmVec['C'-'A'] = 'G';
		char2CharRcmVec['G'-'A'] = 'C';
		char2CharRcmVec['T'-'A'] = 'A';

	}

	/*string stringRcm(const string& s)
	{
		string t;
		for(string::reverse_iterator iter = s.rbegin();
			iter != s.rend(); iter++)
		{
			t = t + char2CharRcmVec(*iter);
		}
		return t;		
	}*/

};

#endif