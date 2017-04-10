#pragma once
/*!
 \file othello.h
 Describes the data structure *l-Othello*...
*/

#include <vector>
#include <iostream>
#include <array>
#include <queue>
#include <cstring>
#include <list>
#include "hash.h"
#include "io_helper.h"
#include "disjointset.h"
using namespace std;

/*!
 * \brief Describes the data structure *l-Othello*. It classifies keys of *keyType* into *2^L* classes.
 * \note Query a key of keyType always return uint64_t, however, only the lowest L bits are meaningful. \n
 * The array are all stored in an array of uint64_t. There are actually m_a+m_b cells in this array, each of length L.
 */
template< class keyType>
class Othello
{
    typedef uint64_t valueType;
#define MAX_REHASH 20 //!< Maximum number of rehash tries before report an error. If this limit is reached, Othello build fails. 
public:
    vector<valueType> mem; //!< actual memory space for arrayA and arrayB.
    uint32_t L; //!< the length of return value.
//    uint32_t LMASK;  return value must be within [0..2^L-1], i.e., LMASK==((1<<L)-1);
#define LMASK ((1ULL<<L)-1)
    uint32_t ma; //!< length of arrayA.
    uint32_t mb; //!< length of arrayB
    Hasher32<keyType> Ha; //<! hash function Ha
    Hasher32<keyType> Hb; //<! hash function Hb
    bool build = false; //!< true if Othello is successfully built.
    uint32_t trycount = 0; //!< number of rehash before a valid hash pair is found.
    DisjointSet disj; //!< Disjoint Set data structure that helps to test the acyclicity.
    vector<uint32_t> fillcount; //!< Enabled only when the values of the query is not pre-defined. This supports othelloIndex query. fillcount[x]==1 if and only if there is an edge pointed to x,
#define FILLCNTLEN (sizeof(uint32_t)*8)
    uint32_t allowed_conflicts; //!< The number of keys that can be skipped during construction.
    vector<keyType> removedKeys; //!< The list of removed keys.
private:
    bool autoclear = false; //!<  clears the memory allocated during construction automatically.
    keyType *keys;
    /*!
     * \brief Get the consecutive L bits starting from location loc*L bit.
     * \warning May return garbage info on the higher bits, needs (& LMASK) afterwards.
     */
    inline valueType get(uint32_t loc) { //
        uint64_t st = ((uint64_t) loc) * L;
        uint32_t mx = st & 0x3F; //mx = st % 64
        if ( L & (L-1) ) {
            //L &(L-1) !=0 ==> L is NOT some power of 2 ==> value may cross the uint64_t barrier.
            if ( (mx + L) > 64 ) {
                //mx+L > 64 means the value cross the barrier.
                //The highest (64-mx) bits of mem[st>>6],
                //the lowest L-(64-mx) bits of mem[(st>>6)+1], shl (64-mx)
                return ((mem[st>>6]>>mx) | (mem[(st>>6)+1]<<(64-mx)));
            };
        }
        return (mem[st>>6]>>mx);
    }



    valueType inline set(uint32_t loc, valueType &value) {
        value &= LMASK;
        uint64_t st = ((uint64_t) loc) * L;
        uint32_t mx = st & 0x3F; //mx = st % 64
        if ( L & (L-1) ) {
            if ( (mx + L) > 64 ) {
                mem[st>>6] &=  ((UINT64_MAX) >> (64-mx));
                mem[st>>6] |=  (value << mx);
                mem[(st>>6)+1] &= (UINT64_MAX << ((mx+L)-64));
                mem[(st>>6)+1] |= (value >> (64-mx));
                return value;
            };
        }
        mem[st>>6] &= (~(LMASK << mx));
        mem[st>>6] |= (value << mx);
        return value;
    }


    inline valueType getrand(valueType &v) {
        valueType va = rand();
        va <<=28;
        va ^= rand();
        va <<=28;
        v = (va ^ rand());
    }

    void newHash() {
        uint32_t s1 = rand();
        uint32_t s2 = rand();
#ifdef HASHSEED1
        s1 = HASHSEED1;
        s2 = HASHSEED2;
#endif
        Ha.setMaskSeed(ma-1,s1);
        Hb.setMaskSeed(mb-1,s2);
        trycount++;
        if (trycount>1) printf("NewHash for the %d time\n", trycount);
    }

    vector<int32_t> *first, *nxt1, *nxt2;
    bool testHash(uint32_t keycount);
    vector<bool> filled;

    /*!
     \brief Fill *Othello* so that the query returns values as defined.
     \param [in] void * values, pointer to the array of values. Each value is of *valuesize* bytes.
     \param [in] uint32_t keycount.
     \param [in] size_t valuesize.
     \note  When *values* is NULL, mark edges as 0 or 1 according to their direction.

    */
    void fillvalue(void *values, uint32_t keycount, size_t valuesize);
    bool trybuild(void *values, uint32_t keycount, size_t valuesize) {
        bool succ;
        disj.setLength(ma+mb);
        if (succ = testHash(keycount)) {
            fillvalue(values, keycount,valuesize);
        }
        if (autoclear || (!succ))
            finishBuild();
        return succ;
    }
public:
    /*!
     \brief Construct *l-Othello*.
     \param [in] keyType *_keys, pointer to array of keys.
     \param [in] uint32_t keycount, number of keys, array size must match this.
     \param [in] bool _autoclear, clear memory used during construction once completed. Forbid future modification on the structure.
     \param [in] void * _values, Optional, pointer to array of values. When *_values* is empty, fill othello values such that the query result classifies keys to 2 sets. See more in notes.
     \param [in] uint32_t valuesize, must be specifed when *_values* is not NULL. This indicates the length of a valueType;
     \param [in] _allowed_conflicts, default value is 10. during construction, Othello will remove at most this number of keys, instead of rehash.
     \note keycount should not exceed 2^29 for memory consideration.
     \n when *_values* is empty, classify keys into two sets X and Y, defined as follow: for each connected compoenents in G, select a node as the root, mark all edges in this connected compoenent as pointing away from the root. for all edges from U to V, query result is 1 (k in Y), for all edges from V to u, query result is 0 (k in X).

    */
    Othello(uint8_t _L, keyType *_keys,  uint32_t keycount, bool _autoclear = true, void *_values = NULL, size_t _valuesize = 0, int32_t _allowed_conflicts = -1 ) {
        allowed_conflicts = _allowed_conflicts;
        L = _L;
        autoclear = _autoclear;
        keys = _keys;
        int hl1 = 8; //start from ma=64
        int hl2 = 7; //start from mb=64
        while ((1UL<<hl2) <  keycount * 1) hl2++;
        while ((1UL<<hl1) < keycount* 1.333334) hl1++;
        ma = (1UL<<hl1);
        mb = (1UL<<hl2);
        mem.resize(1);
        if (allowed_conflicts<0) {
            allowed_conflicts = (ma<20)?0:ma-20;
        }
        uint64_t bitcnt = (ma+mb);
        bitcnt *= L;
        while (bitcnt % (sizeof(valueType)*8)) bitcnt++;
        bitcnt /= (sizeof(valueType)*8);
        mem.resize(bitcnt,0ULL);
        while ( ((uint64_t)mem.size())*sizeof(mem[0])*8ULL<(ma+mb)*((uint64_t)L) )
            mem.push_back(0);


        cout << "Building" <<endl;
        trycount = 0;
        while ( (!build) && (hl1<=31&&hl2<=31)) {
            while ((!build) && (trycount<MAX_REHASH)) {
                newHash();
                build = trybuild( _values, keycount, _valuesize);
            }
            if (!build) {
                if (hl1 == hl2) hl1++;
                else hl2++;
                ma = (1UL<<hl1);
                mb = (1UL<<hl2);
                mem.resize(1);
                uint64_t bitcnt = (ma+mb);
                bitcnt *= L;
                while (bitcnt % (sizeof(valueType)*8)) bitcnt++;
                bitcnt /= (sizeof(valueType)*8);
                mem.resize(bitcnt,0ULL);
                while ( ((uint64_t)mem.size())*sizeof(mem[0])*8ULL<(ma+mb)*((uint64_t)L) )
                    mem.push_back(0);

                cout << "Extend Othello Length to" << human(keycount) <<" Keys, ma/mb = " << human(ma) <<"/"<<human(mb) <<" keyT"<< sizeof(keyType)*8<<"b  valueT" << sizeof(valueType)*8<<"b"<<" L="<<(int) L<<endl;
                trycount = 0;
            }
        }
        printf("%08x %08x\n", Ha.s, Hb.s);
        if (build)
            cout << "Succ " << human(keycount) <<" Keys, ma/mb = " << human(ma) <<"/"<<human(mb) <<" keyT"<< sizeof(keyType)*8<<"b  valueT" << sizeof(valueType)*8<<"b"<<" L="<<(int) L <<" After "<<trycount << "tries"<< endl;

        else
            cout << "Build Fail!" << endl;
    }
    //!\brief Construct othello with vectors.
    template<typename VT>
    Othello(uint8_t _L, vector<keyType> &keys, vector<VT> &values, bool _autoclear = true, int32_t allowed_conflicts = -1) :
        Othello(_L, & (keys[0]),keys.size(), _autoclear, &(values[0]), sizeof(VT), allowed_conflicts)
    {
    }

    //! \brief release memory space used during construction and forbid future modification of arrayA and arrayB.
    void finishBuild() {
        delete nxt1;
        delete nxt2;
        delete first;
        disj.finish();
        filled.clear();
    }

    /*!
        \brief export the information of the *Othello*, not including the array, to a memory space.
        \note memory space length = 0x20. \n
              Exported infomation contains:  \n
              0x00, 32bit, L; \n
              0x04, 32bit, seedB; \n
              0x08, 32bit, seedA; \n
              0x10, 8bit, hlA;
              0x14, 8bit, hlB;
              seedA, seedB, ma , mb (represented as 1<<hl1 and 1<<hl2).
    */
    void exportInfo(unsigned char * v) {
        memset(v,0,0x20);
        uint32_t s1 = Ha.s;
        uint32_t s2 = Hb.s;
        memcpy(v,&L, sizeof(uint32_t));
        memcpy(v+4,&s1,sizeof(uint32_t));
        memcpy(v+8,&s2,sizeof(uint32_t));
        int hl1 = 0, hl2 = 0;
        if (ma == 0 || mb ==0) {
            hl1 = hl2 = 0;
        }
        else {
            while ((1<<hl1)!= ma) hl1++;
            while ((1<<hl2)!= mb) hl2++;
        }
        memcpy(v+0x10,&hl1, sizeof(uint32_t));
        memcpy(v+0x14,&hl2,sizeof(uint32_t));
    }
    /*!
       \brief load the infomation of the *Othello* from memory.
       \note info is exported using *ExportInro()*
     */
    Othello(unsigned char *v) {
        int32_t hl1,hl2;
        int32_t s1,s2;
        memcpy(&(L),v,sizeof(uint32_t));
        memcpy(&(s1),v+4,sizeof(uint32_t));
        memcpy(&(s2),v+8,sizeof(uint32_t));
        memcpy(&hl1, v+0x10, sizeof(uint32_t));
        memcpy(&hl2, v+0x14, sizeof(uint32_t));
        printf("seeds%08x %08x\n", s1,s2);
        if (hl1 > 0 && hl2 >0) {
            ma = (1<<hl1);
            mb = (1<<hl2);
            mem.resize(1);
            uint64_t bitcnt = (ma+mb);
            bitcnt *= L;
            while (bitcnt % (sizeof(valueType)*8)) bitcnt++;
            bitcnt /= (sizeof(valueType)*8);
            mem.resize(bitcnt,0ULL);
            while ( ((uint64_t)mem.size())*sizeof(mem[0])*8ULL<(ma+mb)*((uint64_t)L) )
                mem.push_back(0);
            Ha.setMaskSeed(ma-1,s1);
            Hb.setMaskSeed(mb-1,s2);
        }
        else
            ma = mb =0;
    }

    /*!
       \brief returns a 64-bit integer query value for a key.
    */
    uint64_t queryInt(const keyType &k) {
        uint32_t ha,hb;
        return query(k,ha,hb);
    }



    void printValueTSize() {
        std::cout << sizeof(valueType) << std::endl;
    }


    vector<uint32_t> getCnt(); //!< \brief returns vector, length = 2L. position x: the number of 1s on the x-th lowest bit, for array A, if x<L; otherwise, for arrayB.
    vector<double> getRatio(); //!<\brief returns vector, length = L. position x: the probability that query return 1 on the x-th lowest bit.

    void randomflip(); //!<\brief adjust the array so that for random alien query, returns 0 or 1 with equal probability on each bit.
    /*!
      \brief adjust the array so that for random alien query, return 1 with probability that is close to the *ideal* value.
      \param [in] double ideal. \n ideal = 1.0 means return 1 with higher probability. \n ideal = 0.0 means return 1 with loest probability.
      \note This function is usually able to tune the probabilty within 0.2 ~ 0.8.
     */
    void setAlienPreference(double ideal = 1.0);
private:
    void inline get_hash_1(const keyType &v, uint32_t &ret1) {
        ret1 = (Ha)(v);
    }
public:
    /*!
     \brief compute the hash value of key and return query value.
     \param [in] keyType &k key
     \param [out] uint32_t &ha  computed hash function ha
     \param [out] uint32_t &hb  computed hash function hb
     \retval valueType
     */
    inline valueType query(const keyType &k, uint32_t &ha, uint32_t &hb) {
        get_hash_1(k,ha);
        valueType aa = get(ha);
        get_hash_2(k,hb);
        valueType bb = get(hb);
        //printf("%llx   [%x] %x ^ [%x] %x = %x\n", k,ha,aa&LMASK,hb,bb&LMASK,(aa^bb)&LMASK);
        return LMASK & (aa^bb);
    }
    void inline get_hash_2(const keyType &v, uint32_t &ret1) {
        ret1 = (Hb)(v);
        ret1 += ma;
    }
    void inline get_hash(const keyType &v, uint32_t &ret1, uint32_t &ret2) {
        get_hash_1(v,ret1);
        get_hash_2(v,ret2);
    }
    /*!
     \brief load the array from file.
     \note only the arrayA and B are loaded. This must be called after using constructor Othello<keyType>::Othello(unsigned char *)
     */
    void loadDataFromBinaryFile(FILE *pF) {
        fread(&(mem[0]),sizeof(mem[0]), mem.size(), pF);
    }
    /*!
     \brief write array to binary file.
     */
    void writeDataToBinaryFile(FILE *pF) {
        fwrite(&(mem[0]),sizeof(mem[0]), mem.size(), pF);
    }
private:
    void padd (vector<int32_t> &A, valueType &t) {
        const valueType one = 1;
        for (int i = 0; i <L; i++)
            if  (t & (one<<i))
                A[i]++;
    }

    void pdiff(vector<int32_t> &A, valueType &t) {
        const valueType one = 1;
        for (int i = 0; i <L; i++)
            if  (t & (one<<i))
                A[i]--;
            else A[i]++;
    }
};

/*
template<size_t L, class valueType>

template<size_t L>
std::array<uint8_t,L> operator ^ (const std::array<uint8_t,L>  &A,const std::array<uint8_t,L>  &B) {
    std::array<uint8_t, L> v = A;
    for (int i = 0; i < L; i++) v[i]^=B[i];
    return v;

}

*/

template< class keyType>
bool Othello<keyType>::testHash(uint32_t keycount) {
    uint32_t ha, hb;
    nxt1  = new vector<int32_t> (keycount);
    nxt2  = new vector<int32_t> (keycount);
    first = new vector<int32_t> (ma+mb, -1);
    removedKeys.clear();
    disj.clear();
    for (int i = 0; i < keycount; i++) {
        if (i>1048576) if ((i & (i-1)) == 0) printf("Tesing keys # %d\n",i);
        get_hash(keys[i], ha, hb);
        
        if (disj.sameset(ha,hb)) {
            removedKeys.push_back(keys[i]);
            if (removedKeys.size()> allowed_conflicts) 
                return false;
            continue;
        }
        (*nxt1)[i] = (*first)[ha];
        (*first)[ha] = i;
        (*nxt2)[i] = (*first)[hb];
        (*first)[hb] = i;
        disj.merge(ha,hb);
    }
    return true;
}

template< class keyType>
void Othello<keyType>::fillvalue(void *values, uint32_t keycount, size_t valuesize) {
    list<uint32_t> Q;
    vector<int32_t> *nxt;
    filled.resize(ma+mb);
    fill(filled.begin(), filled.end(), false);
    if (values == NULL) {
        fillcount.resize((ma+mb)/32);
        fill(fillcount.begin(),fillcount.end(),0);
    }
    for (int i = 0; i< ma+mb; i++)
        if (disj.isroot(i)) {
            while (!Q.empty()) Q.pop_front();
            Q.push_back(i);
            valueType vv;
            getrand(vv);
            set(i,vv);
            filled[i] = true;
            while (!Q.empty()) {
                uint32_t nodeid = (*Q.begin());
                Q.pop_front();
                if (nodeid < ma) nxt = nxt1;
                else nxt = nxt2;
                int32_t kid = first->at(nodeid);
                while (kid >=0) {
                    uint32_t ha,hb;
                    get_hash(keys[kid],ha,hb);
                    if (filled[ha] && filled[hb]) {
                        kid = nxt->at(kid);
                        continue;
                    }
                    int helse = filled[ha] ? hb : ha;
                    int hthis = filled[ha] ? ha : hb;
                    //! m[hthis] is already filled, now fill m[helse].
                    valueType valueKid = 0;
                    if (values != NULL) {
                        uint8_t * loc = (uint8_t *) values;
                        loc += (kid * valuesize);
                        memcpy(&valueKid, loc, valuesize);
                    }
//                        valueKid = values[kid];
                    else {
                        //! when hthis == ha, this is a edge pointing from U to V, i.e., value of this edge shall be set as 1.
                        valueKid = (hthis == ha)?1:0;
                        fillcount[helse/FILLCNTLEN] |= (1<<(helse % FILLCNTLEN));
                    }

                    valueType newvalue = valueKid ^ get(hthis);
//                    printf("%x %x %x===", valueKid, get(hthis), newvalue);
                    set(helse, newvalue);
//                    printf("k%llx ha/hb %lx %lx: %x ^ %x = %x ^ %x = %x (%x), helse%x ,i%x, %d\n",
//                       keys[kid], ha, hb, get(hthis) & LMASK, newvalue & LMASK, get(ha) &LMASK, get(hb) &LMASK, (get(ha)^get(hb)) & LMASK, valueKid &LMASK, helse ,hthis ,(bool)((valueKid &LMASK)==((get(ha)^(get(hb)))&LMASK)));
                    Q.push_back(helse);
                    filled[helse] = true;
                    kid = nxt->at(kid);
                }

            }
        }
}

template< class keyType>
vector<uint32_t> Othello<keyType>::getCnt() {
    vector<uint32_t> cnt;
    cnt.resize(L+L);
    for (int i = 0; i < ma; i++) {
        valueType gv = get(i);
        uint8_t *vv;
        vv = (uint8_t*) ((void *) &gv);
        for (int p = 0; p < L; p++)  {
            uint8_t tv;
            tv =  ((*vv) >> (p & 7));
            if (tv & 1) cnt[p]++;
            if ((p & 7) == 7) vv++;
        }
    }
    for (int i = ma; i < ma+mb; i++) {
        valueType gv = get(i);
        uint8_t *vv;
        vv = (uint8_t*) ((void *) &gv);
        for (int p = 0; p < L; p++)  {
            uint8_t tv;
            tv =  ((*vv) >> (p & 7));
            if (tv & 1) cnt[L+p]++;
            if ((p & 7) == 7) vv++;
        }
    }
    return cnt;
}


template< class keyType>
vector<double> Othello<keyType>::getRatio() {
    vector<uint32_t> cnt = getCnt();
    vector<double> ret;
    ret.resize(L);
    for (int i = 0; i < L; i++) {
        double p1 = 1.0 * cnt[i] / ma;
        double p2 = 1.0 * cnt[i+L] / mb;
        ret[i] = p1*(1-p2)+p2*(1-p1);
    }
    return ret;

}

template< class keyType>
void Othello<keyType>::randomflip() {
    if (filled.size() <=1) return;
    printf("Random flip\n");
    valueType vv;
    vector<list<uint32_t> > VL(ma+mb, list<uint32_t>());
    for (int i = 0; i< ma+mb; i++)
        if (!filled[i]) {
            valueType vv;
            getrand(vv);
            set(i,vv);
        }
        else {
            VL[disj.getfa(i)].push_back(i);
        }
    for (int i = 0; i < ma+mb; i++) {
        getrand(vv);
        for (auto j = VL[i].begin(); j!=VL[i].end(); j++) {
            valueType newvalue =  get(*j)^vv;
            set(*j, newvalue);
        }

    }
}
double getrate(uint32_t ma, uint32_t mb, uint32_t da, uint32_t db) {
    double pa = da*1.0/ma;
    double pb = db*1.0/mb;
    double ret = pa*(1-pb)+pb*(1-pa);
    return ret;
}


template< class keyType>
void Othello<keyType>::setAlienPreference(double ideal) {
    int da[] = {1,1,0,-1,-1,-1,0,1};
    int db[] = {0,1,1,1,0,-1,-1,-1};
    vector< array<int32_t,8> > sa (L, array<int32_t,8>());
    vector< array<int32_t,8> > sb (L, array<int32_t,8>());
    for (int i = 0; i < 8; i++)
        for (int j =0; j <L; j++)
            sa[j][i] = sb[j][i] =0;
    vector<int32_t> na(L,0),nb(L,0);
    //for (int j =0; j <L; j++) na[j]=nb[j] =0;
    int emptyA = 0;
    int emptyB = 0;
    if (filled.size() <=1) return;
    printf("Adjust bitmap goal: return 1 with rate %.3lf\n",ideal);
    valueType vv;
    vector<list<uint32_t> > VL(ma+mb, list<uint32_t>());
    for (int i = 0; i< ma+mb; i++) {
        if (!filled[i]) {
            valueType vv;
            if (i<ma) emptyA++;
            else emptyB++;
        }
        else {
            VL[disj.getfa(i)].push_back(i);
            valueType cur = get(i);
            if (i<ma) padd(na, cur);
            else padd(nb,cur);
        }
    }
    for (int i = 0; i < ma+mb; i++) {
        vector<int32_t> diffa(L,0),diffb(L,0);
//        for (int j = 0; j< L; j++) diffa[j]=diffb[j] =0;
        for (auto j = VL[i].begin(); j!=VL[i].end(); j++) {
            valueType cur = get(*j);
            if (*j < ma) pdiff(diffa,cur);
            else pdiff(diffb, cur);
        }
        for (int bitID = 0; bitID<L; bitID++)
            for (int j = 0; j <8; j++)
                if (da[j]*diffa[bitID] + db[j] *diffb[bitID] > 0)
                {
                    sa[bitID][j]+=diffa[bitID];
                    sb[bitID][j]+=diffb[bitID];
                }

    }
    vector<int32_t> direction(L,0);
    for (int bitID = 0; bitID <L; bitID++) {
        double ratemin = 1.0;
        for (int dir = 0; dir <8*4; dir++) {
            int setA = ((dir & 8) !=0);
            int setB = ((dir & 16)!=0);
            double rate = getrate(ma,mb,na[bitID]+setA*emptyA+sa[bitID][dir & 0x7], nb[bitID]+setB*emptyB+sb[bitID][dir & 0x7]);
            if ((rate-ideal)*(rate-ideal) < ratemin) {
                ratemin = (rate-ideal)*(rate-ideal);
                direction[bitID] = dir;
            }
        }
    }
    valueType veA=0, veB=0;

    for (int bitID = 0; bitID<L; bitID++) {
        veA |= ((direction[bitID] & 8) ? (1<<bitID) : 0);
        veB |= ((direction[bitID] & 16) ? (1<<bitID) : 0);
    }
    for (int i = 0; i < ma; i++) if (!filled[i])
            set(i,veA);
    for (int i = ma; i <ma+mb; i++) if (!filled[i])
            set(i,veB);

    for (int i = 0; i < ma+mb; i++) {
        vector<int32_t> diffa(L,0),diffb(L,0);
        for (auto j = VL[i].begin(); j!=VL[i].end(); j++) {
            valueType cur = get(*j);
            if (*j < ma) pdiff(diffa,cur);
            else pdiff(diffb, cur);
        }
        valueType vv =0;
        for (int bitID = 0; bitID<L; bitID++)
            if (da[direction[bitID] &7 ] * diffa[bitID] + db[direction[bitID] &7] * diffb[bitID] > 0)
                vv |= (1<<bitID);

        for (auto j = VL[i].begin(); j!=VL[i].end(); j++) {
            valueType newvalue =  get(*j)^vv;
            set(*j, newvalue);
        }

    }




}
