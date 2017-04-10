#pragma once
#include <vector>
using namespace std;
/*! \file disjointset.h
 *  Disjoint Set data structure.  
 */
/*!
 * \brief Disjoint Set data structure. Helps to test the acyclicity of the graph during construction. 
 * */ 
class DisjointSet {
    int siz;
    vector<int32_t> *fa;
public:
    uint32_t getfa(int i) {
        if ((*fa)[i] < 0) (*fa)[i] = i;
        else if ((*fa)[i]!=i)
            (*fa)[i] = getfa((*fa)[i]);
        return (*fa)[i];
    }
    //! Release the memory to save some space.
    void finish() {
        fa->clear();
        delete fa;
    }
    void setLength(int n) {
        fa = new vector<int32_t> (n,-1);
    }
    //! re-initilize the disjoint sets.
    void clear() {
        for (auto a : *fa)
            a = -1;
    }
    void merge(int a, int b) {
        if (a==0) swap(a,b);  //a!=0
        (*fa)[getfa(b)] = getfa(a);
    }
    bool sameset(int a, int b) {
        return getfa(a)==getfa(b);
    }
    bool isroot(int a) {
        return ((*fa)[a]==a);
    }
};
