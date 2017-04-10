#pragma once
#include "othelloindex.h"
const char * VERSION = GITVERSION;
/*!
 \file mulothindex.h
 Describes the data structure *MulOthIndex*

 */

/*!
 \brief MulOthIndex is a extended version of OthelloIndex, it supports more than 2^29 keys. \n
  The keys are first classifed into groups, each group is then maintained by one *OthelloIndex*.
  \note Query a key of keyType always return uint64_t, however, only the lowest bits are meaningful.
*/

template <class keyType>
class MulOthIndex {
    typedef uint64_t valueType;
    vector<OthelloIndex<keyType> *> vOthIdxs; //!< list of *OthelloIndex*s.
    vector<valueType> shift; 
    unsigned char split;
    bool addOthIndex(unsigned int groupid, vector<keyType> &keys) {
        OthelloIndex<keyType> *poth;
        poth = new OthelloIndex<keyType>(&keys[0], keys.size());
        if (!poth->build) {
            printf("Build Halt!\n");
            return false;
        }
        vOthIdxs[groupid] = poth;
        for (int t = groupid + 1; t< shift.size(); t++)
           shift[t] += keys.size();
		if ((poth-> removedKeys).size()>0){
			for (keyType k : poth->removedKeys) {
                keyType fullk; helper->combgrp(fullk,groupid,k);
                printf("Key in Grp %x : %llx is removed \n", groupid, fullk);
                removedKeys.push_back(fullk);
			}
		}
        return true;
    }
    IOHelper<keyType,valueType> *helper;
public:
    bool buildsucc;
	vector<keyType> removedKeys; //!< list of skipped keys for all underlying *OthelloIndex*.
    /*!
     \brief Construct a Mul-OthelloIndex from a file.
     \param char * fname  the file contains key/value pairs, each in a line. Although only the keys are used.
     \param unsigned char _split  the keys are first classifed according to their highest *_split* bits. There are in total 2^_split groups. Classification method as described in *splitgrp*.
	 \param IOHelper _helper identifies how to convert a raw data to keytype and group.
     \param bool fileIsSorted When fileIsSorted, assume that the file is sorted so that the keys appear in the ascending order of groups.
     * */
    MulOthIndex(const char * fname, unsigned char _split, class IOHelper<keyType,valueType> * _helper, bool fileIsSorted = false) : helper(_helper) {
        split = _split;
        printf("Building MulOthelloIndex from file %s\n", fname);
        FILE *pFile;
        pFile = fopen (fname, "r");
        vOthIdxs. clear();
        vOthIdxs.resize(1<<_split);
        shift. clear();
        shift.resize(1<<_split);
        char buf[1024];
        if (split ==0) {
            keyType k;
            vector<keyType> keys;
            while (true) {
                if (fgets(buf,1024,pFile)==NULL) break;
                if (!helper->convert(buf,&k)) break;
                keys.push_back(k);
            }
            buildsucc = addOthIndex(0,keys);
            fclose(pFile);
            return;
        }
        if (fileIsSorted) {
            uint32_t grpid = 0;
            printf("Reading file for keys in group %2x/%2x\n", grpid,(1<<split)-1);
            vector<keyType> keys;
            while (true) {
                keyType k;
                void * pp;
                if (fgets(buf,1024,pFile)==NULL) break;
                if (!helper->convert(buf, &k)) break;
                keyType keyingroup;
                uint32_t groupid;
                helper->splitgrp(k,groupid,keyingroup);
                if (groupid != grpid) {
                    if (!addOthIndex(grpid,keys))
                    {
                        fclose(pFile);
                        return;
                    }
                    grpid = groupid;
                    printf("Reading file for keys in group %2x/%0x\n", grpid,(1<<split)-1);
                    keys.clear();
                }
                keys.push_back(keyingroup);
            }
            if (!addOthIndex(grpid,keys))
            {
                fclose(pFile);
                return;
            }
        }
        else {
            for (uint32_t grpid = 0; grpid < (1<<_split); grpid++) {
                printf("Reading file for keys in group %2x/%2x\n", grpid,(1<<split)-1);
                vector<keyType> keys;
                rewind(pFile);
                while (true) {
                    keyType k;
                    if (fgets(buf,1024,pFile)==NULL) break;
                    if (!helper->convert(buf, &k)) break;
                    keyType keyingroup;
                    uint32_t groupid;
                    helper->splitgrp(k,groupid,keyingroup);
                    if (groupid != grpid) continue;
                    keys.push_back(keyingroup);
                }
                printf("keycount %d ", keys.size());
                if (keys.size()>0)
                    if (!addOthIndex(grpid,keys)) {
                        fclose(pFile);
                        return;
                    }

            }
        }
        buildsucc = true;
        fclose(pFile);
    }

    inline valueType query(const keyType &k) {
        uint32_t grp;
        keyType kgrp;
        if (split ==0)
            return vOthIdxs[0]->query(k);
        else {
            helper->splitgrp(k,grp,kgrp);
            return (vOthIdxs[grp]!=NULL)?shift[grp]+(vOthIdxs[grp]->query(kgrp)):0;
        }
    }
    /*! \brief write MulOthelloIndex to a file.
      \note  file structure. \n
      0x00 32bit split \n
      0x04 version (reserved); \n
      0x20 32Byte info for *OthelloIndex* 0 \n
      0x40 32Byte info for *OthelloIndex* 1 \n
      ..\n
      ..\n
      Array for OthelloIndex 0 \n
      ..\n
      32bit nRemovedKeys \n
      (64bit) keyType *removedKeys* 0\n

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
            if (vOthIdxs[i]!=NULL) 
                vOthIdxs[i]->exportInfo(buf0x20);
            else 
                memset(buf0x20,0, sizeof(buf0x20));
            fwrite(buf0x20,sizeof(buf0x20),1,pFile);
        }
        fwrite(&shift[0],sizeof(shift[0]),1<<split,pFile);
        for (int i = 0; i <(1<<split); i++) {
            if (vOthIdxs[i]!=NULL)
            vOthIdxs[i]->writeDataToBinaryFile(pFile);
//           fwrite(&(vOths[i]->mem[0]),sizeof(vOths[i]->mem[0]), vOths[i]->mem.size(), pFile);
        }
        
        uint32_t nRemovedKeys = removedKeys.size();
        fwrite(&nRemovedKeys,sizeof(nRemovedKeys),1,pFile);
        for (int i = 0 ; i < nRemovedKeys; i++)
            fwrite(&removedKeys[i],sizeof(keyType),1,pFile);
        fclose(pFile);
    }

    //! \brief construct a Grouped l-Othello Index from a file.
    MulOthIndex(const char* fname, IOHelper<keyType,valueType> * _helper): helper(_helper) {
        buildsucc = false;
        printf("Read from binary file %s\n", fname);
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
        vOthIdxs.clear();
        for (int i = 0; i < (1<<split); i++) {
            fread(buf0x20,sizeof(buf0x20),1,pFile);
            OthelloIndex<keyType> *poth = NULL;
            uint64_t *ppoth; ppoth = (uint64_t*) &buf0x20;
            if (*ppoth) 
                poth = new OthelloIndex<keyType>(buf0x20);
            vOthIdxs.push_back(poth);
        }
		shift.resize(1<<split);
        fread(&shift[0],sizeof(shift[0]),1<<split,pFile);
        for (int i = 0; i < (1<< split); i++) {
            if (vOthIdxs[i] != NULL) {
                vOthIdxs[i]->loadDataFromBinaryFile(pFile);
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
    ~MulOthIndex() {
		for (auto p: vOthIdxs)
			delete p;
	}
};
