/*!
 * \file io_helper.h
 * Contains IO utilities.
 */
#pragma once
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <errno.h>
#include <queue>
#include <algorithm>
#include <cstring>
/*!
 * \brief interface for converting a key from its raw format to a keyTypeype. Split key into groups.
 */
template <typename keyType, typename valueType>
class IOHelper {
public:
    /*!
       \brief convert a input-style line to key/value pair.
       \param [in] char *s, a line from the input.
       \param [out] keyType* T, converted key.
       \param [out] valueType* V, converted value.
       \retval boolean return true if convert is success.
      */
    virtual bool convert(char *s, keyType *T, valueType *V)  = 0;
    //! \brief skip the value.
    virtual bool convert(char *s, keyType *T)  = 0;

    /*!
     \brief split a keyTypeype value into two parts: groupID/keyInGroup by the highest *splitbit* bits.
    */
    virtual void splitgrp(const keyType &key, uint32_t &grp, keyType &keyInGroup) = 0;
    /*!
     \brief combine groupID/keyInGroup to the origional key
     */
    virtual void combgrp(keyType &key, uint32_t &grp, keyType &keyInGroup) = 0;
};

template<typename keyType, typename valueType>
class FileReader {
public:
    IOHelper<keyType, valueType> *helper;
    virtual bool getFileIsSorted() = 0;
    virtual bool getNext(keyType *T, valueType *V) = 0;
    virtual void finish() =0;
    virtual void reset() = 0;
    virtual ~FileReader() {
    }
};

/*!
 * \brief IOHelper for Constant-Length Kmers.
 * \note Each Kmer is a string of length KMERLENGTH, \n
 * We consider string as a base-4 number. [A=0,C=1,G=2,T=3]. \n
 * Kmers are grouped according to the highest *splitbit* bits.
 */
template <typename keyType, typename valueType>
class ConstantLengthKmerHelper : public IOHelper<keyType,valueType> {
public:
    uint8_t kmerlength; //!< Assume all kmers are of the same length.
    uint8_t splitbit;   //!< group the keys according to the highest bits.
    ConstantLengthKmerHelper(uint8_t _kmerlength, uint8_t _splitbit): kmerlength(_kmerlength),splitbit(_splitbit) {};
    ~ConstantLengthKmerHelper() {};
    inline bool convert(char *s, keyType *k, valueType *v) {
        char *s0;
        s0 = s;
        switch (*s) {
        case 'A':
        case 'T':
        case 'G':
        case  'C':
            keyType ret = 0;
            while (*s == 'A' || *s == 'C' || *s =='T' || *s =='G') {
                ret <<=2;
                switch (*s) {
                case 'T':
                    ret++;
                case 'G':
                    ret++;
                case 'C':
                    ret++;
                }
                s++;
            }
            *k = ret;
            valueType tv;
            sscanf(s,"%d",&tv);
            *v = tv;
            return true;

        }
        return false;

    }
    inline bool convert( char *s, keyType *k) {
        valueType v;
        return convert(s,k,&v);
    }
    void splitgrp(const keyType &key, uint32_t &grp, keyType &keyInGroup) {
        int mvcnt = 2 * kmerlength - splitbit;
        keyType high = (key >> mvcnt);
        grp = high;
        keyType lowmask = 1;
        lowmask <<= mvcnt;
        keyInGroup = (key & (lowmask-1));
    }

    void combgrp(keyType &key, uint32_t &grp, keyType &keyInGroup) {
        key = grp;
        key <<= (2*kmerlength - splitbit);
        key |= (keyInGroup);
    }

};



//! convert a 64-bit Integer to human-readable format in K/M/G. e.g, 102400 is converted to "100K".
std::string human(uint64_t word) {
    std::stringstream ss;
    if (word <= 1024) ss << word;
    else if (word <= 10240) ss << std::setprecision(2) << word*1.0/1024<<"K";
    else if (word <= 1048576) ss << word/1024<<"K";
    else if (word <= 10485760) ss << word*1.0/1048576<<"M";
    else if (word <= (1048576<<10)) ss << word/1048576<<"M";
    else ss << word*1.0/(1<<30) <<"G";
    std::string s;
    ss >>s;
    return s;
}

//! split a c-style string with delimineter chara.
std::vector<std::string> split(const char * str, char deli) {
    std::istringstream ss(str);
    std::string token;
    std::vector<std::string> ret;
    while(std::getline(ss, token, deli)) {
        if (token.size()>=1)
            ret.push_back(token);
    }
    return ret;
}
#include <cstring>

template <typename keyType, typename valueType>
class KmerFileReader : public FileReader<keyType,valueType> {
    FILE *f;
    bool fIsSorted;
public:
    KmerFileReader(const char *fname, IOHelper<keyType,valueType> *_helper, bool b)  {
        fIsSorted = b;
        FileReader<keyType,valueType>::helper = _helper;
        char buf[1024];
        strcpy(buf,fname);
        if (buf[strlen(buf)-1]=='\n')
            buf[strlen(buf)-1] = '\0';
        f=fopen(buf,"r");
        printf("OpenFile %s %x\n",fname,f);
    }
    void finish() {
        fclose(f);
    }
    void reset() {
        rewind(f);
    }
    bool getFileIsSorted() {
        return fIsSorted;
    }
    ~KmerFileReader() {
        finish();
    }
    bool getNext(keyType *T, valueType *V) {
        char buf[1024];
        if (fgets(buf,sizeof(buf),f)==NULL) return false;
        return FileReader<keyType,valueType>::helper->convert(buf,T,V);
    }
};

template <typename keyType, typename valueType>
struct KVpair {
    keyType k;
    valueType v;
    bool friend operator <( const KVpair &a, const KVpair &b) {
        return a.k > b.k;
    }
} __attribute__((packed));



template <typename keyType, typename valueType>
class compressFileReader : public FileReader <keyType, valueType> {
    FILE *f;
    bool fIsSorted;
    static const int buflen = 1024;
    int curr = 0;
    int  max = 0;
    unsigned char buf[1024*64];
    uint32_t kl, vl;
public:
    compressFileReader( const char * fname, IOHelper <keyType, valueType> * _helper, uint32_t klength, uint32_t valuelength, bool _fIsSorted = true) {
        kl = klength;
        vl = valuelength;
        FileReader<keyType,valueType> :: helper = _helper;
        fIsSorted = _fIsSorted;
        char buf[1024];
        strcpy(buf,fname);
        if (buf[strlen(buf)-1]=='\n')
            buf[strlen(buf)-1] = '\0';
        f=fopen(buf,"rb");
        printf("OpenFile to binary read Kmers %s %x\n",fname,f);
    }
    void finish() {
        fclose(f);
    }
    void reset() {
        printf("Do not support reset()\n");
    }
    bool getFileIsSorted() {
        return fIsSorted;
    }
    ~ compressFileReader() {
        finish();
    }
    bool getNext(keyType *k, valueType *v) {
        if (curr == max) {
            max = fread(buf,kl+vl,buflen,f);
            if (max == 0) return false;
            curr = 0;
        }
        *k =0;
        *v =0;
        memcpy( (void *) k, buf + curr*(kl+vl), kl);
        memcpy( (void *) v, buf + curr*(kl+vl)+kl, vl);
        curr++;
        return true;
    }

};


template <typename keyType, typename valueType>
class MultivalueFileReaderWriter : public FileReader <keyType, valueType> {
    FILE *f;
    static const int buflen = 8192;
    int curr = 0;
    int max = 0;
    unsigned char buf[buflen * 2];
    uint32_t kl,vl;
    bool isRead;
    bool isclosed = false;
public:
    static const valueType EMPTYVALUE = ~0;
    bool valid(uint32_t value) {
        if (vl == 1) return value!=0xFF;
        if (vl == 2) return value!=0xFFFF;
        if (vl == 4) return value!=0xFFFFFFFFUL;
    }
    MultivalueFileReaderWriter( const char * fname, uint32_t klength, uint32_t valuelength, bool _isRead) {
        kl = klength;
        vl = valuelength;
        char buf[1024];
        strcpy(buf,fname);
        if (buf[strlen(buf)-1]=='\n')
            buf[strlen(buf)-1] = '\0';
        if (isRead = _isRead)
            f=fopen(buf,"rb");
        else f = fopen(buf,"wb");
        memset(buf,0,sizeof(buf));
        printf("OpenFile to binary read/write Multivalue Kmers file %s %x\n",fname,f);
        curr = 0;
        max = 0;
    }
    void finish() {
        if (!isclosed) {
        if (!isRead) {
            fwrite(buf,sizeof(buf[0]), curr, f);
        }
        fclose(f);
        }
        isclosed = true;
    }
    void reset() {
        rewind(f);
    }
    bool getFileIsSorted() {
        return false;
    }
    ~MultivalueFileReaderWriter() {
        finish();
    }
    void getmore() {
            memmove(buf, buf+curr, max - curr);
            max -= curr;
            curr = 0;
            if (max < buflen) {
                max += fread(buf+max,1,buflen, f);
            }
    }
    bool get(void * mem, uint32_t l) {
        if (curr + l >= max) getmore();
        if (curr + l > max) return false;
        memcpy(mem,buf+curr,l);
        curr+=l;
        return true;
    }
    bool getNext(keyType *k, valueType *v) {
        if (!get(k,kl)) return false;
        while (true) {
            get(v,vl);
            if (!valid(*v)) return true;
            v++;
        }
        /*
        if (curr + kl+vl >= max) {
            getmore();
            if (max < kl) return false;
        }
        *k = 0;
        memcpy( (void *) k, buf + curr, kl);
        curr+=kl;
        memcpy( (void *) v, buf + curr, vl);
        curr += vl;
        while (valid(*v)) {
            v++;
            if (curr+ vl + vl >= max) getmore();
            memcpy( (void *) v, buf + curr, vl);
            curr += vl;
        }
        return true;
        */
    }
    void add(void * mem, uint32_t l) {
        memcpy(buf+curr, mem, l);
        curr += l;
        if (curr >= buflen) {
            fwrite( buf, 1, buflen, f);
            if (curr > buflen) {
                memcpy(buf,buf+buflen,curr-buflen);
            }
            curr -= buflen;
        }
    }
    void write(keyType *k ,valueType *v) {
        add((void *) k, kl);
        while (valid(*v)) {
            add((void *) v, vl);
            v++;
        }
    }
    void write(keyType *k, std::vector<valueType> &vv) {
        add ( (void *) k, kl);
        void *p; uint8_t a8; uint16_t a16; uint32_t a32;
        if (vl == 1) p = &a8;
        if (vl == 2) p = &a16;
        if (vl == 4) p = &a32;
        for (auto v: vv) {
            a8 = a16 = a32 = v;
            add(p,vl);
        }
        a8 =0xFF;  a16 = 0xFFFF; a32 = 0xFFFFFFFF;//~0ULL;
        add (p,vl);
    }
};

using namespace std;

template <typename keyType>
class KmerReader {
public:
    virtual void finish() =0;
    virtual bool getNext(keyType *k) =0;
};

template <typename KVpair>
class BinaryKmerReader: public KmerReader<KVpair> {
    FILE * f;
    static const int buflen = 16;
    KVpair buff[1024];
    int curr = 0;
    int max = 0;
    bool isclosed =false;
public:
    BinaryKmerReader(const char * fname) {
        char buf[1024];
        strcpy(buf,fname);
        if (buf[strlen(buf)-1]=='\n')
            buf[strlen(buf)-1] = '\0';
        f=fopen(buf,"rb");
        printf("OpenFile to read Kmers %s %x\n",fname,f);
        if (f==0) {
            printf("errno %d %s\n",errno,strerror(errno));
        }
        curr = 0;
    }
    void finish() {
        if (!isclosed) {
        fclose(f);
        }
        isclosed = true;
    }
    ~BinaryKmerReader() {        finish();    }
    bool getNext(KVpair *ret) {
        if (curr == max) {
            max = fread(buff,sizeof(buff[0]),buflen,f);
            *ret = (KVpair) ~0ULL;
            if (max == 0) return false;
            curr = 0;
        }
        memcpy(ret, &buff[curr], sizeof(buff[curr]));
        curr++;
        return true;
    }
};

template <typename KVpair>
class BinaryKmerWriter {
    FILE *f;
    int curr = 0;
public:
    KVpair buf[1024];
    BinaryKmerWriter( const char * fname) {
        char buf[1024];
        strcpy(buf,fname);
        if (buf[strlen(buf)-1]=='\n')
            buf[strlen(buf)-1] = '\0';
        f=fopen(buf,"wb");
        printf("OpenFile to write Kmers %s %x\n",fname,f);
        curr = 0;
        memset(buf,0,sizeof(buf));
    }
    static const int buflen = 16;
    void write(KVpair *p) {
        memcpy(&buf[curr],p,sizeof(buf[curr]));
        curr++;
        if (curr == buflen) {
            fwrite(buf,sizeof(buf[0]),buflen,f);
            curr = 0;
        }
    }
    void finish() {
        fwrite(buf,sizeof(buf[0]),curr,f);
        curr = 0;
        fclose(f);
    }
};

//! read kmer from unsorted txt file and sort .
template <typename keyType> 
class SortedKmerTxtReader : public KmerReader<keyType> {
    BinaryKmerReader<keyType> * binaryReader = NULL; 
    uint32_t pointer;
    vector<keyType> * vK;
    public: 
    SortedKmerTxtReader(const char * fname, uint32_t kmerlength, const char *tmpfilename) {
        ConstantLengthKmerHelper<keyType, uint64_t> helper(kmerlength, 0);
        FileReader<keyType,uint64_t> *reader;
        reader = new KmerFileReader<keyType,uint64_t>(fname, &helper, false);
        keyType k; uint64_t v;
        vK = new vector<keyType>();
        while (reader->getNext(&k, &v)) {
            vK->push_back(k);
        }
        sort(vK->begin(),vK->end());
        delete reader;
        if (tmpfilename != NULL) {
            string binaryfilename (tmpfilename);
            BinaryKmerWriter<keyType> writer(binaryfilename.c_str());
            for (uint64_t k:*vK)    
                writer.write(&k);
            writer.finish();
            binaryReader = new BinaryKmerReader<keyType> (binaryfilename.c_str());
        }
        else pointer = 0;
    }
    ~SortedKmerTxtReader() {
        finish();
        if (binaryReader)
            delete binaryReader;
    }
    bool getNext(keyType *k) {
        if (binaryReader!= NULL)
            return binaryReader->getNext(k);
        else {
            if (pointer == vK->size()) {
                delete vK;
                return false;
            }
            *k = (*vK)[pointer++];
            return true;
        }
    }
    void finish() {
        if (binaryReader)
            binaryReader->finish();
    }
};


template <typename keyType, typename valueType>
class taxoTreeBuilder: public FileReader <keyType, valueType> {
    vector< FILE *> fV;
    vector<compressFileReader <keyType, valueType> *> readerV;
    vector< vector< int > > NCBI; //texonomyToNCBIID;
    vector< int > stID;
    struct KIDpair {
        keyType k;
        uint32_t id;
        bool finished;
        bool friend operator <( const KIDpair &a, const KIDpair &b) {
            if (a.finished != b.finished) return (((int) a.finished) > ((int) b.finished));
            return a.k>b.k;
        }
    };
public:
    void finish() {
        for (auto f: fV) fclose(f);
    }
    void reset() {
        printf(" Do not support reset() \n");
    }
    int levelcount;
    vector<vector<int> > NCBI_local;
    vector<int> localshift;
    vector<vector<string> > NCBI_ID;
    vector<KmerReader<uint64_t> *> readers;
    vector<MultivalueFileReaderWriter<uint64_t, uint16_t> *> grpreaders; //must be 64-bit kmers, and 16-bit grpids.
    priority_queue<KIDpair> PQ;
    bool combineMode = false; //used when there are >=800 files;
    uint32_t combineCount; // split the file into combineCount groups,
    bool getFileIsSorted() {
        return true;
    }
    void groupFile(string fname, vector<string> lf, string prefix, string suffix, int32_t idshift, bool useBinaryKmerFile,uint32_t KmerLength, const char * tmpfolder) {
        vector<KmerReader<keyType> *> readers;
        priority_queue<KIDpair> PQN;
        for (string s: lf) {
            string fname = prefix + s + suffix;
            if (useBinaryKmerFile) 
                readers.push_back(new BinaryKmerReader<keyType>(fname.c_str()));
            else {
                string tmpfname(tmpfolder); tmpfname = tmpfname + s + ".bintmp";
                readers.push_back(new SortedKmerTxtReader<keyType>(fname.c_str(),KmerLength,NULL));
            }
            keyType key; 
            readers[readers.size()-1]->getNext(&key);
            KIDpair kid = {key, idshift+readers.size()-1, false};
            PQN.push(kid);
        }
        
        MultivalueFileReaderWriter<keyType,uint16_t> * writer = new MultivalueFileReaderWriter<keyType,uint16_t> (fname.c_str(),8,2,false); 
        // Loop key for these files;
        while (true) {
            keyType key = PQN.top().k;
            uint32_t id = PQN.top().id;
            vector<uint16_t> ret;
            if (PQN.top().finished) {
                for (auto r: readers) { 
                    r->finish();
                    delete r;
                }
                writer->finish();
                delete writer;
                return;
            }
            while (PQN.top().k == key && !PQN.top().finished) {
                int tid = PQN.top().id;
                ret.push_back(tid);
                keyType nextk;
                bool finish = !readers[tid-idshift]->getNext(&nextk);
                PQN.pop();
                KIDpair kid = {nextk, tid, finish};
                PQN.push(kid);
            }
            writer->write(&key, ret); 
        }
    }
    vector< vector<uint16_t> > grpTmpValue;

    taxoTreeBuilder(const char * NCBIfname, const char * fnameprefix, const char * fnamesuffix, const char * tmpFileDirectory, uint32_t KmerLength, uint32_t splitbit, bool useBinaryKmerFile = true ) {
        FileReader<keyType,valueType>::helper = new ConstantLengthKmerHelper<keyType,valueType> (KmerLength,splitbit);
        FILE * fNCBI;
        string prefix ( fnameprefix);
        string suffix (fnamesuffix);
        fNCBI = fopen(NCBIfname, "r");
        //Assuming the file is tab-splited,
        //Species_index	Species_ID	Species_name	Genus_index	Genus_ID	Genus_name	Family_index	Family_ID	Family_name	Order_index	Order_ID	Order_name	Class_index	Class_ID	Class_name	Phylum_index	Phylum_ID	Phylum_name
        char buf[4096];
        fgets(buf, 4096, fNCBI); // skip the first line
        vector<string> vv = split(buf, '\t');
        levelcount = vv.size()/3;
        NCBI_local.clear();
        NCBI_ID.clear();
        NCBI_local.resize(levelcount);
        NCBI_ID.resize(levelcount);
        readers.clear();
        vector<string> fnames;
        while (true) {
            if (fgets(buf, 4096, fNCBI) == NULL) break; // read a Species
            vector<string> vv = split(buf, '\t');
            if (vv.size()<2) break;
            for (int i = 0 ; i*3 < vv.size(); i++) {
                int localID = atoi(vv[i*3].c_str());
                NCBI_local[i].push_back(localID);
                NCBI_ID[i].push_back(vv[i*3+1]);
            }
            fnames.push_back(vv[1]);
        }
        localshift.clear();
        localshift.push_back(1);
        for (int i = 0; i < levelcount; i++)
            localshift.push_back(localshift[i] + *max_element(NCBI_local[i].begin(), NCBI_local[i].end())+1);

        int nn = 50;
        combineMode = (fnames.size()>nn);
        if (combineMode) {
            int curr = 0;
            int combineCount = 0;
            vector<string> * fnamesInThisgrp ;
            vector<string> grpfnames;
            while (curr < fnames.size()) {
                if (curr + nn < fnames.size())
                    fnamesInThisgrp = new vector<string> (fnames.begin()+curr, fnames.begin()+curr+nn);
                else
                    fnamesInThisgrp = new vector<string> (fnames.begin()+curr, fnames.end());
                stringstream ss;
                string tmpFolder(tmpFileDirectory);

                ss<<tmpFolder<<"TMP"<<grpfnames.size();
                string fnamegrp;
                ss>> fnamegrp;
                grpfnames.push_back(fnamegrp);
                printf("merge kmer files %d %d to grp %s\n", curr, curr+fnamesInThisgrp->size()-1, fnamegrp.c_str());
                groupFile(fnamegrp, *fnamesInThisgrp, prefix, suffix, curr, useBinaryKmerFile,KmerLength,tmpFileDirectory);
                curr += fnamesInThisgrp->size();
                delete fnamesInThisgrp;
            }
            combineCount = grpfnames.size();
            for (string v: grpfnames) {
                grpreaders.push_back( new MultivalueFileReaderWriter<uint64_t, uint16_t>(v.c_str(), 8,2, true));
                keyType key;
                uint16_t valuebuf[1024];
                grpreaders[grpreaders.size()-1]->getNext(&key, valuebuf);
                vector<uint16_t> Vvaluebuf;
                for (int i = 0 ; grpreaders[0]->valid(valuebuf[i]); i++)
                   Vvaluebuf.push_back(valuebuf[i]);
                grpTmpValue.push_back(Vvaluebuf);
                KIDpair kid = {key, grpreaders.size()-1, false};
                PQ.push(kid);
            }
        }
        else
            for (int i = 0 ; i < NCBI_ID.size(); i++) {
                string fname = prefix + fnames[i] + suffix;
                if (useBinaryKmerFile)
                    readers.push_back(new BinaryKmerReader<keyType>(fname.c_str()));
                else {
                    string tmpfname(tmpFileDirectory); tmpfname = tmpfname + fnames[i] + ".bintmp";
                    readers.push_back(new SortedKmerTxtReader<keyType>(fname.c_str(),KmerLength,tmpfname.c_str()));
                }
                keyType key;
                readers[readers.size()-1]->getNext(&key);
                KIDpair kid = {key, readers.size()-1, false};
                PQ.push(kid);
            }
        fclose(fNCBI);
        string IDLfname(tmpFileDirectory); IDLfname+= "IDList.txt";
        FILE * IDLf; IDLf = fopen(IDLfname.c_str(),"w");
        for (int t : localshift) {
            fprintf(IDLf,"%d\n",t);
        } 
        fclose(IDLf);
    }
    ~taxoTreeBuilder() {
        if (combineMode)  {
            for (int i = 0 ; i < grpreaders.size(); i++)
                delete grpreaders[i];
        }
        else 
            for (int i = 0 ; i < readers.size(); i++)
                delete readers[i];
        delete FileReader<keyType,valueType>::helper;
    }
    bool getNext( keyType *k, valueType *v) {
        int anslevel = 0;
        keyType key = PQ.top().k;
        vector<int> ret;
        if (PQ.top().finished) {
            finish();
            return false;
        }
       // printf("Find key %llx:", key);
        while (PQ.top().k == key && !PQ.top().finished) {
            int tid;
            tid = PQ.top().id;
            keyType nextk;
            bool finish;
            if (combineMode) {
                ret.insert(ret.end(),grpTmpValue[tid].begin(),grpTmpValue[tid].end());
                int ll = grpTmpValue[tid].size();
         //       printf("   %d keys: (from %d)\t", ll, tid);
           //     for (int i: grpTmpValue[tid])
           //          printf("%x\t",i);
                uint16_t valuebuf[1024];
                finish = !grpreaders[tid]->getNext(&nextk, valuebuf);
                grpTmpValue[tid].clear();
                for (int i = 0; grpreaders[tid]->valid(valuebuf[i]); i++)
                    grpTmpValue[tid].push_back(valuebuf[i]);
            //    printf("Next Has ::%d::", grpTmpValue[tid].size());
            }
            else {
                ret.push_back(tid);
              //  printf(" %x\t",PQ.top().id);
                finish = !readers[tid]->getNext(&nextk);
            }
            PQ.pop();
            KIDpair kid = {nextk, tid, finish};
            PQ.push(kid);
        }
        *k = key;
         
        for (int i = 0; i< levelcount; i++) {
            bool flag = true;
            for (int j = 0; j < ret.size() && flag; j++)
                flag = (NCBI_local[i][ret[j]]==NCBI_local[i][ret[0]]);
            if (flag) {
                *v = localshift[i] + NCBI_local[i][ret[0]];
                return true;
            }
        }
        *v = localshift[levelcount];
        return true;
    }
};
