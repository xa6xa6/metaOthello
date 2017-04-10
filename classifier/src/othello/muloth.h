#pragma once
/*!
 * \file muloth.h
 * Describes the data structure *Grouped l-Othello*.
 * */
#include "othello.h"
#include "io_helper.h"
#include <cstdio>
#include <cstdlib>

#include <cstring>
using namespace std;
const char * VERSION = GITVERSION;
//#pragma message "VERSION: " GITVERSION

// split(k,plow,phigh,spl) (((plow)=(k & ((1ULL << spl)-1))), ((phigh)=(k >> spl)))

/*!
 * \brief Describes the data structure *Grouped l-Othello*. It classifies keys of *keyType* into *2^L* classes. \n
 * The keys are first classifed into groups, each group is then maintained by one *l-Othello*. 
 * \note Query a key of keyType always return uint64_t, however, only the lowest L bits are meaningful. \n
 */
template <typename keyType, typename IOvalueType = uint16_t>
class MulOth {
    typedef uint64_t valueType;
    uint32_t L;
    vector<Othello<keyType> *> vOths; //!< list of *l-Othellos*
    unsigned char split;
    template<typename VT>
    bool addOth(unsigned int groupid, vector<keyType> &keys, vector<VT> &values) {
        Othello<keyType> *poth;
        poth = new Othello<keyType>(L, &keys[0], keys.size(), true, &values[0], sizeof(values[0]));
        if (!poth->build) {
            printf("Build Halt!\n");
            return false;
        }
        vOths[groupid] = poth;
        if ((poth -> removedKeys).size()>0) {
            for (keyType k : poth->removedKeys) {
                keyType fullk; helper->combgrp(fullk,groupid,k);
                printf("Key in Grp %x : %llx is removed \n", groupid, fullk);
                removedKeys.push_back(fullk);
            }
        }
        return true;
    }
public:
    vector<keyType> removedKeys; //!< list of skipped keys for all underlying *Othello*.
    bool buildsucc; 
    IOHelper<keyType,IOvalueType> *helper;

    //!\brief This generates toy data for test purpose.
    /*
    MulOth( uint32_t _L, uint32_t NN) {
        split = 0;
        L = _L;
        vector<keyType> keys;
        vector<uint8_t> values;
        for (int i = 0 ; i < NN; i ++) {
            keys.push_back((((uint64_t) i) << 32) + i + 1);
            values.push_back(i & 0xFF);
        }
        vOths.clear();
        vOths.resize(1);
        addOth(0,keys,values);
    }
    */
    /*!
     \brief Construct a Grouped l-Othello from a file.
     \param uint32_t _L
     \param char * fname  the file contains key/value pairs, each in a line.
     \param unsigned char _split  the keys are first classifed according to their highest *_split* bits. There are in total 2^_split groups. Classification method as described in *splitgrp*.
	 \param IOHelper _helper identifies how to convert a raw data to keytype and group.
     \param bool fileIsSorted When fileIsSorted, assume that the file is sorted so that the keys appear in the ascending order of groups.
     * */
    MulOth(uint32_t _L, unsigned char _split, class FileReader<keyType, uint16_t> * _reader) {
        helper = _reader->helper;
        L = _L;
        vOths.clear();
        vOths.resize(1<<_split);
        char buf[1024];
        split = _split;
        if (split ==0) {
            keyType k; IOvalueType v; 
            vector<keyType> keys;
            vector<IOvalueType> values;
            while (true) {
                if (!_reader->getNext(&k,&v)) break;
                keys.push_back(k);
                values.push_back(v);
            }
            buildsucc = addOth(0,keys,values);
            _reader->finish();
            return;
        }
        
        if (_reader->getFileIsSorted())  {
            uint32_t grpid = 0;
            printf("Reading file for keys in group %2x/%2x\n", grpid,(1<<split)-1);
            vector<keyType> keys;
            vector<IOvalueType> values;
            while (true) {
                keyType k; IOvalueType v;  keyType keyingroup;  uint32_t groupid;
                if (!_reader->getNext(&k,&v)) break;
                _reader->helper->splitgrp(k,groupid,keyingroup);
                if (groupid != grpid) {
                    if (!addOth(grpid,keys,values)) 
                        {_reader->finish(); return;}
                    grpid = groupid;
                    printf("Reading file for keys in group %2x/%0x\n", grpid,(1<<split)-1);
                    keys.clear();
                    values.clear();
                }
                keys.push_back(keyingroup);
                values.push_back(v);
            }
            if (!addOth(grpid,keys,values)) 
                { _reader->finish(); return;}
        }
        else
            for (uint32_t grpid = 0; grpid < (1<<_split); grpid++) {
                printf("Reading file for keys in group %2x/%2x\n", grpid,(1<<split)-1);
                vector<keyType> keys;
                vector<IOvalueType> values;
                _reader->reset();
                while (true) {
                    keyType k; IOvalueType v;  keyType keyingroup;  uint32_t groupid;
                    if (!_reader->getNext(&k,&v)) break;
                    _reader->helper->splitgrp(k,groupid,keyingroup);
                    if (groupid != grpid) continue;
                    keys.push_back(keyingroup);
                    values.push_back(v);
                }
                printf("keycount %d ", keys.size());
                if (keys.size()>0)
                    if (!addOth(grpid,keys,values)) 
                        {_reader->finish(); return;}
            }

        buildsucc = true;
        _reader->finish();
    }
    MulOth(uint32_t _L, const char * fname, unsigned char _split, class IOHelper<keyType,IOvalueType> * _helper, bool fileIsSorted = false) :
        MulOth(_L, _split, new KmerFileReader<keyType,IOvalueType>( fname, _helper, fileIsSorted))
    {

    }

    /*!
       \brief returns a *L*-bit integer query value for a key.
      */
    inline IOvalueType query(const keyType &k) {
        uint32_t grp;
        keyType kgrp;
        if (split ==0) 
            return vOths[0]->queryInt(k);
        else {
            helper->splitgrp(k,grp,kgrp);
            return (vOths[grp]!=NULL)?vOths[grp]->queryInt(kgrp):0;
        }
    }
    void printall () {
        printf("Printall ...\n");
        for (auto V : vOths)
            V.printall();
    }
    /*! \brief write Grouped l-Othello to a file.
     \note file structure. \n
          0x00 32bit split \n
          0x04 version (reserved); \n
          0x20 32Byte info for othello 0 \n
          0x40 32Byte info for othello 1 \n
          ... \n
          ... \n
          Array for othello 0 \n
          .., \n
          32bit nRemovedKeys \n
          (64bit) keyType &removedKeys* 0\n

     * */
    void writeToFile(const char* fname) {
        FILE *pFile;
        pFile = fopen (fname, "wb");
        unsigned char buf0x20[0x20];
        char *p;
        p = (char *) &(buf0x20[0]);
        memset(buf0x20,0, sizeof(buf0x20));
        strcpy(p+0x4,VERSION);
        uint32_t split32 = split;
        memcpy(buf0x20, &split32, sizeof(uint32_t));
        fwrite(buf0x20,sizeof(buf0x20),1,pFile);
        for (int i = 0; i <(1<<split); i++) 
        {
            if (vOths[i]!=NULL) 
                vOths[i]->exportInfo(buf0x20);
            else 
                memset(buf0x20,0, sizeof(buf0x20));
            fwrite(buf0x20,sizeof(buf0x20),1,pFile);
        }
        for (int i = 0; i <(1<<split); i++) {
            if (vOths[i]!=NULL)
            vOths[i]->writeDataToBinaryFile(pFile);
//           fwrite(&(vOths[i]->mem[0]),sizeof(vOths[i]->mem[0]), vOths[i]->mem.size(), pFile);
        }
        
        uint32_t nRemovedKeys = removedKeys.size();
        fwrite(&nRemovedKeys,sizeof(nRemovedKeys),1,pFile);
        for (int i = 0 ; i < nRemovedKeys; i++)
            fwrite(&removedKeys[i],sizeof(keyType),1,pFile);

        fclose(pFile);
    }

    //! \brief construct a Grouped l-Othello from a file.
    MulOth(const char* fname, IOHelper<keyType,IOvalueType> * _helper): helper(_helper) {
        buildsucc = false;
        printf("Read from binary Othello file %s\n", fname);
        FILE *pFile;
        pFile = fopen (fname, "rb");
        uint32_t compversion;
        unsigned char buf0x20[0x20];
        fread(buf0x20,sizeof(buf0x20),1,pFile);
        char *p;
        p = (char *) &(buf0x20[0]);
#ifndef NO_VERSION_CHECK
        if (strcmp(p+4,VERSION)) {
            printf("Wrong version number\n");
            fclose(pFile);
            return;
        }
#endif
        uint32_t split32;
        memcpy(&split32, buf0x20, sizeof(uint32_t));
        split = split32;
        vOths.clear();
        for (int i = 0; i < (1<<split); i++) {
            fread(buf0x20,sizeof(buf0x20),1,pFile);
            Othello<keyType> *poth = NULL;
            uint64_t *ppoth; ppoth = (uint64_t*) &buf0x20;
            if (*ppoth) 
                poth = new Othello<keyType>(buf0x20);
            vOths.push_back(poth);
        }
        for (int i = 0; i < (1<< split); i++) {
            if (vOths[i] != NULL) {
                vOths[i]->loadDataFromBinaryFile(pFile);
                L= vOths[i]->L;
            }
            //fread(&(vOths[i]->mem[0]),sizeof(vOths[i]->mem[0]), vOths[i]->mem.size(), pFile);
        }
        uint32_t nRemovedKeys = removedKeys.size();
        fread(&nRemovedKeys,sizeof(nRemovedKeys),1,pFile);
        removedKeys.resize(nRemovedKeys);
        for (int i = 0 ; i < nRemovedKeys; i++)
            fread(&removedKeys[i],sizeof(keyType),1,pFile);

        fclose(pFile);
        buildsucc = true;
    }
    ~MulOth() {
        for (auto p: vOths)
            delete p;
    }
};


//MulOthello binary file descriptor
//0x00 : uint32_t splitbit
//0x04 : uint32_t signature MulOth version
//0x20 : OthInfo1
//0x40 : OthInfo2
//...
//offset1 : Oth[0].m
//offset2 = offset1 + Oth[0].hashupperbound = Oth[1].m
//...

//OthInfo: 32 Byte
//0x00 : uint64_t hash1
//0x08 : uint64_t hash2
//0x10 : uint32_t mask1
//0x14 : uint32_t mask2
//0x18 : uint64_t m.offset

