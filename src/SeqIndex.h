/* 
 * File:   SeqIndex.h
 * Author: vujevic
 *
 * Created on April 25, 2014, 11:42 PM
 */

#include<string>
#include<vector>


using namespace std;

#ifndef SEQINDEX_H
#define	SEQINDEX_H

/**
 * Class which represents five pentapeptides for every sequence.
 */
class SeqIndex {
public:
    vector<string> getPentapeptides();
    vector<string> getLowestPentapeptides();
    void setPentapeptides(vector<string> pent, int listSize);
    void setLowestPentapeptides(vector<string>pent,int listSize);    
    /**
     * Method which check similarty of two @pentaptides. Vectors are similar if they contains three or more same 
     * pentapeptides.
     * @param other 
     * @return @true if pentapeptides are similar, else @false
     */
    bool checkSimilarity(SeqIndex other);
    
    /**
     * Methods which returns numbers of hits with pentapeptide, return frequency of @pentapeptid in @pentapeptides.
     * @return number of hits.
     */
    int numberOfHits(string pentapeptid);

		int regSize;
private:
    /**
     * Indexes
     */
    vector<string> pentapeptides;
    vector<string> lowestPentapeptides;
};

#endif	/* SEQINDEX_H */
