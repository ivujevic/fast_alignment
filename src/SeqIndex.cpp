
/* 
 * File:   SeqIndex.cpp
 * Author: vujevic
 * 
 * Created on April 25, 2014, 11:42 PM
 */

#include "SeqIndex.h"
#include <vector>
#include <algorithm>

        
using namespace std;

vector<string> SeqIndex::getPentapeptides() {
    return pentapeptides;
}
vector<string> SeqIndex::getLowestPentapeptides() {
	return lowestPentapeptides;
}

void SeqIndex::setPentapeptides(vector<string> pent,int listSize) {
  if(pent.size() > 0 ) {

    //get pentapeptide with the largest frequency
    string mostFrequent = pent.at(0);

    int diffToListSize = listSize - pent.size();

    // if number of pentapeptides are less than five.
    while(diffToListSize > 0) {
      pent.insert(pent.begin(),mostFrequent);
      diffToListSize--;
    }

    pentapeptides = pent;
  }//if(pent.size() == 0) pentapeptides.assign(listSize,"Q");
}
void SeqIndex::setLowestPentapeptides(vector<string> pent,int listSize) {
  if(pent.size() > 0 ) {

    //get pentapeptide with the largest frequency
    string mostFrequent = pent.at(0);

    int diffToListSize =listSize - pent.size();

    // if number of pentapeptides are less than five.
    while(diffToListSize > 0) {
      pent.insert(pent.begin(),mostFrequent);
      diffToListSize--;
    }

    lowestPentapeptides = pent;
  }//if(pent.size() == 0) pentapeptides.assign(listSize,"Q");
}
bool SeqIndex::checkSimilarity(SeqIndex other) {
    vector<string> otherPentapeptides = other.getPentapeptides();
    
    int numberOfSame =0;
    
    for (vector <string>::iterator it = otherPentapeptides.begin(); it != otherPentapeptides.end(); it++) {

        if(std::find(pentapeptides.begin(),pentapeptides.end(),*it) != pentapeptides.end()) {
            numberOfSame ++;
        }
    }
    
    return numberOfSame > 2;
}

int SeqIndex::numberOfHits(string pentapeptid) {
    int counter = 0;
    
    for(vector<string>::iterator it = pentapeptides.begin(); it!= pentapeptides.end(); it++) {
        if(*it == pentapeptid) counter ++;
    }
    
    return counter;
}
