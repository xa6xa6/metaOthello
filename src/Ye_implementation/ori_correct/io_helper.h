/*!
 * \file io_helper.h
 * Contains IO utilities.
 */
#pragma once
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
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
        case 'A':        case 'T':        case 'G':        case  'C':
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
    inline bool convert( char *s, keyType *k){
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


