#pragma once
#include "othello.h"
/*!
  \file othelloindex.h
  Describes the data structure *OthelloIndex*
*/

/*!
    \brief OthelloIndex is a data structure that stores a minimum perfect hash function. 
    \n For n keys, the query value of keys range in [0..n-1]
    Query always returns uint32_t.
    \note n should not exceed 2^29 for memory consideration. For even larger n, use the grouped version *MulOthIndex*. \n
    For better performance, the we use POPCNT instructions, which is supported by most mordern CPUs.\n
    The index of a key k is either Sum {cntfill[0..Ha(k)-1]} or Sum {cntfill[0..Hb(k)-1]}. This is decided by query result on *oth*.

 */
template <class keyType>
class OthelloIndex : public Othello<keyType> {
    vector<uint32_t> offset;
    uint32_t sum(uint32_t h) {
        uint32_t p = (h & 0x1F);
        return offset[h>>5] + __builtin_popcount(Othello<keyType>::fillcount[h>>5] & ((1<<p)-1));
    }
public:    
    /*! 
     \brief Construct *OthelloIndex*
     \param [in] keyType * keys, pointer to array of keys.
     \param [in] uint32_t keycount, number of keys.
     \note For n keys, the query value of keys range in [0..n-1]
     */
    OthelloIndex(keyType *_keys, uint32_t keycount) :
        Othello<keyType>(1,_keys, keycount)
        {
        offset.resize(Othello<keyType>::fillcount.size());
        offset[0] = 0;
        for (int i = 1 ; i < offset.size(); i++) {
            offset[i] = offset[i-1] + __builtin_popcount(Othello<keyType>::fillcount[i-1]);
        }
          
    }
    OthelloIndex(unsigned char *v): 
        Othello<keyType>(v) {
        Othello<keyType>::fillcount.resize((Othello<keyType>::ma+Othello<keyType>::mb)/32);
        offset.resize(Othello<keyType>::fillcount.size());
        }
    /*!
     \brief return the index of key k, in range [0..n-1]
     \param [in] keyType k
     \retval uint32_t index
     */
    uint32_t query(const keyType &k) {
        uint32_t ha,hb;
        uint64_t v = Othello<keyType>::query(k,ha,hb);
        uint32_t h = v?hb:ha;
        return sum(h);
    }

    /*!
     \brief load the array from file. 
     \note  the arrayA and B are loaded, in addition, read the array *fillcount* and array *offset*. This must be called after using constructor Othello<keyType>::Othello(unsigned char *)
     */
    void loadDataFromBinaryFile(FILE *pF) {
        fread(&(Othello<keyType>::mem[0]),sizeof(Othello<keyType>::mem[0]), Othello<keyType>::mem.size(), pF);
        fread(&(Othello<keyType>::fillcount[0]),sizeof(Othello<keyType>::fillcount[0]), Othello<keyType>::fillcount.size(), pF);
        fread(&(offset[0]),sizeof(offset[0]),offset.size(),pF);
    }
    /*!
     \brief write array to binary file.
     */
    void writeDataToBinaryFile(FILE *pF) {
        fwrite(&(Othello<keyType>::mem[0]),sizeof(Othello<keyType>::mem[0]), Othello<keyType>::mem.size(), pF);
        fwrite(&(Othello<keyType>::fillcount[0]),sizeof(Othello<keyType>::fillcount[0]), Othello<keyType>::fillcount.size(), pF);
        fwrite(&(offset[0]),sizeof(offset[0]),offset.size(),pF);
    }
};



