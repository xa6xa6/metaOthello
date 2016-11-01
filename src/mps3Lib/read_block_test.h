#ifndef READS_FILE_H
#define READS_FILE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
//#include <hash_map>
#include <map>
#include <set>
#include <sstream>
#include <omp.h>
#include <time.h>

using namespace std;

inline string int_to_str(int numerical)
{
		char c[100];
		sprintf(c,"%d",numerical);
		string str(c);
		return str;
}

inline string int_to_str_sstream(int numerical)
{
	//char c[12];
	//char* intStr = itoa(numerical);
	//string str = string(intStr);
	//return str;
	stringstream ss;
	ss << numerical;
	string str = ss.str();
	return str;
}

inline char complement(int i) 
{
	static const int b2c_size = 20;
	static const char b2c[] = {'T','N','G','N','N','N','C','N','N','N','N','N','N','N','N','N','N','N','N','A'};
	static const char b2cl[] = {'t','n','g','n','n','n','c','n','n','n','n','n','n','n','n','n','n','n','n','a'};
	if (i - 'A' >= 0 && i - 'A' < b2c_size)
		return b2c[i - 'A'];
	else if (i - 'a' >= 0 && i - 'a' < b2c_size)
		return b2cl[i - 'a'];
	else return 'N';
}
/* //X:does not work
inline string revcomp(const string& s) 
{
	string r;
	transform(s.begin(), s.end(), back_inserter(r), complement);
	reverse(r.begin(), r.end());
	return r;
}
*/
inline string revcomp(string s) //X: rewrite 06/08 
{
	string t;
	for(string::reverse_iterator iter = s.rbegin();
		iter != s.rend(); iter++)
	{
		t = t + complement(*iter);
	}
	return t;
}

class Read_Block
{
public:
	Read_Block(bool _paired)
	{ 
		paired = _paired;
	}
	
	~Read_Block()
	{
		clear();
	}

	string get_read_id()
	{
		return read_id;
	}

	string get_seg_seq(int seg_no)
	{
		return read_seq[seg_no-1];
	}

	string get_revcom_seg_seq(int seg_no)
	{
		return read_seq_revcom[seg_no-1];
	}

	string get_seg_qual(int seg_no)
	{
		return read_quality[seg_no-1];
	}

	string get_revcom_seg_qual(int seg_no)
	{
		return read_quality_revcom[seg_no-1];
	}

	int get_seg_len(int seg_no)
	{
		return (int)(read_seq[seg_no - 1 ].length());
	}
	
	int get_seg_len(int seg_no1, int seg_no2)
	{
		if(seg_no1 < 1 || seg_no2 < 1 || seg_no1 > (int)seg_num || seg_no2 > (int)seg_num)
			return 0;
		if(seg_no1 == (int)seg_num || seg_no2 == (int)seg_num)
			return (abs(seg_no1 - seg_no2)) * get_seg_len(1) + get_seg_len(seg_num);
		else		
		  return (abs(seg_no1 - seg_no2) + 1) * get_seg_len(1);
	}
	
   /////////size_t version
	string get_seg_seq(size_t seg_no)
	{
		return read_seq[seg_no-1];
	}

	string get_revcom_seg_seq(size_t seg_no)
	{
		return read_seq_revcom[seg_no-1];
	}

	string get_seg_qual(size_t seg_no)
	{
		return read_quality[seg_no-1];
	}

	string get_revcom_seg_qual(size_t seg_no)
	{
		return read_quality_revcom[seg_no-1];
	}

	int get_seg_len(size_t seg_no)
	{
		return (int)(read_seq[seg_no - 1 ].length());
	}
	
	int get_seg_len(size_t seg_no1, size_t seg_no2)
	{
		if(seg_no1 < 1 || seg_no2 < 1 || seg_no1 > seg_num || seg_no2 > seg_num)
			return 0;
		if(seg_no1 == seg_num || seg_no2 == seg_num)
			return (abs((int)seg_no1 - (int)seg_no2)) * get_seg_len(1) + get_seg_len(seg_num);
		else		
		  return (abs((int)seg_no1 - (int)seg_no2) + 1) * get_seg_len(1);
	}
	
	void clear()
	{
		read_id.clear();
		read_strand.clear();
		read_seq.clear();
		read_seq_revcom.clear();
		read_quality.clear();
		read_quality_revcom.clear();
		seg_num=0;
	}

	void set_seg_num()
	{
	 	seg_num = read_seq.size();
	}
	
	int get_seg_num()
	{
		return (int)seg_num;
	}
	
	void get_full_read_seq()
	{
		for(size_t i = 0; i < read_seq.size(); i++)
		{
			full_read_seq.append(read_seq[i]);
		}	
	}
	
	string read_id;
	string read_strand;
	string full_read_seq;
	vector<string> read_seq;
	vector<string> read_seq_revcom;
	vector<string> read_quality;
	vector<string> read_quality_revcom;
	size_t seg_num;
	bool paired;
};


#endif