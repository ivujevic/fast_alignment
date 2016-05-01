/* 
 * File:   Seg.h
 * Author: vujevic
 *
 * Created on October 15, 2014, 8:24 PM
 */

#ifndef SEG_H
#define	SEG_H

#include <string>
#include <math.h>
#include <vector>
#include "Constants.h"
using namespace std;

enum class Type{PROTEINS,NUCLEOTIDES}; // Protein sequence or nucleotide

class Seg{
private:
	Type type;
	int L = windowSize;
	double k1 = lowerEntropy;
	double k2 = upperEntropy;
	char maskChar = 'X';
	int maxtrim = 100;
	int N;
	int downset;
	int upset;

	void getEntropies(string& seq,int len,vector<double>& H);
	double logFactorial(int n);
	int findLo(int i, int limit,vector<double>& H);
	int findHi(int i,int limit,vector<double>& H);
	double getProb(string& seq,int L);

	void trim(string& seq1,int left,int right,int& leftEnd,int& rightEnd);
	void segSeq(string& seq,vector< pair<int,int>>& segs, int offset);
public:
	Seg();
	Seg(Type type);
	void mask(string& seq);
};

#endif	/* SEG_H */
