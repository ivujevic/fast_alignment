#include "Constants.h"
#include <boost/unordered_map.hpp>
			// A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X   Y   Z
char reductionMap [26] = {'K','L','B','A','A','G','C','E','E','L','A','E','F','A','L','J','A','A','K','K','L','E','I','L','H','L'};

/**
 * 3D array that represents codon table. First dimension represents first nucleotide 
 * in codons, second dimension second nucleotide and so on
 */

char codonTable1[4][4][4] = {
		//						T								   C									A									 G
	/*T*/ { {'F','F','L','L'}, {'S','S','S','S'}, {'Y','Y','*','*'}, {'C','C','*','W'} },
	/*C*/ { {'L','L','L','L'}, {'P','P','P','P'}, {'H','H','Q','Q'}, {'R','R','R','R'} },
	/*A*/ { {'I','I','I','M'}, {'T','T','T','T'}, {'N','N','K','K'}, {'S','S','R','R'} },
	/*G*/ { {'V','V','V','V'}, {'A','A','A','A'}, {'D','D','E','E'}, {'G','G','G','G'}}
};
	
int BLOSUM_45_TABLE[26 * 26] = {
    5, -1, -1, -2, -1, -2, 0, -2, -1, -1, -1, -1, -1, -1, -5, -1, -1, -2, 1, 0, -5, 0, -2, -1, -2, -1,
    -1, 5, -2, 6, 1, -3, -1, 0, -3, -3, 0, -3, -2, 5, -5, -2, 0, -1, 0, 0, -5, -3, -4, -1, -2, 1,
    -1, -2, 12, -3, -3, -2, -3, -3, -3, -2, -3, -2, -2, -2, -5, -4, -3, -3, -1, -1, -5, -1, -5, -1, -3, -3,
    -2, 6, -3, 7, 2, -4, -1, 0, -4, -3, 0, -3, -3, 2, -5, -1, 0, -1, 0, -1, -5, -3, -4, -1, -2, 1,
    -1, 1, -3, 2, 6, -3, -2, 0, -3, -3, 1, -2, -2, 0, -5, 0, 2, 0, 0, -1, -5, -3, -3, -1, -2, 5,
    -2, -3, -2, -4, -3, 8, -3, -2, 0, 1, -3, 1, 0, -2, -5, -3, -4, -2, -2, -1, -5, 0, 1, -1, 3, -3,
    0, -1, -3, -1, -2, -3, 7, -2, -4, -4, -2, -3, -2, 0, -5, -2, -2, -2, 0, -2, -5, -3, -2, -1, -3, -2,
    -2, 0, -3, 0, 0, -2, -2, 10, -3, -2, -1, -2, 0, 1, -5, -2, 1, 0, -1, -2, -5, -3, -3, -1, 2, 0,
    -1, -3, -3, -4, -3, 0, -4, -3, 5, 4, -3, 2, 2, -2, -5, -2, -2, -3, -2, -1, -5, 3, -2, -1, 0, -3,
    -1, -3, -2, -3, -3, 1, -4, -2, 4, 4, -3, 4, 2, -3, -5, -3, -2, -3, -2, -1, -5, 2, -2, -1, 0, -2,
    -1, 0, -3, 0, 1, -3, -2, -1, -3, -3, 5, -3, -1, 0, -5, -1, 1, 3, -1, -1, -5, -2, -2, -1, -1, 1,
    -1, -3, -2, -3, -2, 1, -3, -2, 2, 4, -3, 5, 2, -3, -5, -3, -2, -2, -3, -1, -5, 1, -2, -1, 0, -2,
    -1, -2, -2, -3, -2, 0, -2, 0, 2, 2, -1, 2, 6, -2, -5, -2, 0, -1, -2, -1, -5, 1, -2, -1, 0, -1,
    -1, 5, -2, 2, 0, -2, 0, 1, -2, -3, 0, -3, -2, 6, -5, -2, 0, 0, 1, 0, -5, -3, -4, -1, -2, 0,
    -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 1, -5, -5, -5, -5, -5, 1, -5, -5, -5, -5, -5,
    -1, -2, -4, -1, 0, -3, -2, -2, -2, -3, -1, -3, -2, -2, -5, 9, -1, -2, -1, -1, -5, -3, -3, -1, -3, -1,
    -1, 0, -3, 0, 2, -4, -2, 1, -2, -2, 1, -2, 0, 0, -5, -1, 6, 1, 0, -1, -5, -3, -2, -1, -1, 4,
    -2, -1, -3, -1, 0, -2, -2, 0, -3, -3, 3, -2, -1, 0, -5, -2, 1, 7, -1, -1, -5, -2, -2, -1, -1, 1,
    1, 0, -1, 0, 0, -2, 0, -1, -2, -2, -1, -3, -2, 1, -5, -1, 0, -1, 4, 2, -5, -1, -4, -1, -2, 0,
    0, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, 0, -5, -1, -1, -1, 2, 5, -5, 0, -3, -1, -1, -1,
    -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 1, -5, -5, -5, -5, -5, 1, -5, -5, -5, -5, -5,
    0, -3, -1, -3, -3, 0, -3, -3, 3, 2, -2, 1, 1, -3, -5, -3, -3, -2, -1, 0, -5, 5, -3, -1, -1, -3,
    -2, -4, -5, -4, -3, 1, -2, -3, -2, -2, -2, -2, -2, -4, -5, -3, -2, -2, -4, -3, -5, -3, 15, -1, 3, -2,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5, -1, -1, -1, -1, -1, -5, -1, -1, -1, -1, -1,
    -2, -2, -3, -2, -2, 3, -3, 2, 0, 0, -1, 0, 0, -2, -5, -3, -1, -1, -2, -1, -5, -1, 3, -1, 8, -2,
    -1, 1, -3, 1, 5, -3, -2, 0, -3, -2, 1, -2, -1, 0, -5, -1, 4, 1, 0, -1, -5, -3, -2, -1, -2, 5
};

int BLOSUM_50_TABLE[26 * 26] = {
    5, -2, -1, -2, -1, -3, 0, -2, -1, -8, -1, -2, -1, -1, -8, -1, -1, -2, 1, 0, -8, 0, -3, -1, -2, -1,
    -2, 5, -3, 5, 1, -4, -1, 0, -4, -8, 0, -4, -3, 4, -8, -2, 0, -1, 0, 0, -8, -4, -5, -1, -3, 2,
    -1, -3, 13, -4, -3, -2, -3, -3, -2, -8, -3, -2, -2, -2, -8, -4, -3, -4, -1, -1, -8, -1, -5, -2, -3, -3,
    -2, 5, -4, 8, 2, -5, -1, -1, -4, -8, -1, -4, -4, 2, -8, -1, 0, -2, 0, -1, -8, -4, -5, -1, -3, 1,
    -1, 1, -3, 2, 6, -3, -3, 0, -4, -8, 1, -3, -2, 0, -8, -1, 2, 0, -1, -1, -8, -3, -3, -1, -2, 5,
    -3, -4, -2, -5, -3, 8, -4, -1, 0, -8, -4, 1, 0, -4, -8, -4, -4, -3, -3, -2, -8, -1, 1, -2, 4, -4,
    0, -1, -3, -1, -3, -4, 8, -2, -4, -8, -2, -4, -3, 0, -8, -2, -2, -3, 0, -2, -8, -4, -3, -2, -3, -2,
    -2, 0, -3, -1, 0, -1, -2, 10, -4, -8, 0, -3, -1, 1, -8, -2, 1, 0, -1, -2, -8, -4, -3, -1, 2, 0,
    -1, -4, -2, -4, -4, 0, -4, -4, 5, -8, -3, 2, 2, -3, -8, -3, -3, -4, -3, -1, -8, 4, -3, -1, -1, -3,
    -8, -8, -8, -8, -8, -8, -8, -8, -8, 1, -8, -8, -8, -8, 1, -8, -8, -8, -8, -8, 1, -8, -8, -8, -8, -8,
    -1, 0, -3, -1, 1, -4, -2, 0, -3, -8, 6, -3, -2, 0, -8, -1, 2, 3, 0, -1, -8, -3, -3, -1, -2, 1,
    -2, -4, -2, -4, -3, 1, -4, -3, 2, -8, -3, 5, 3, -4, -8, -4, -2, -3, -3, -1, -8, 1, -2, -1, -1, -3,
    -1, -3, -2, -4, -2, 0, -3, -1, 2, -8, -2, 3, 7, -2, -8, -3, 0, -2, -2, -1, -8, 1, -1, -1, 0, -1,
    -1, 4, -2, 2, 0, -4, 0, 1, -3, -8, 0, -4, -2, 7, -8, -2, 0, -1, 1, 0, -8, -3, -4, -1, -2, 0,
    -8, -8, -8, -8, -8, -8, -8, -8, -8, 1, -8, -8, -8, -8, 1, -8, -8, -8, -8, -8, 1, -8, -8, -8, -8, -8,
    -1, -2, -4, -1, -1, -4, -2, -2, -3, -8, -1, -4, -3, -2, -8, 10, -1, -3, -1, -1, -8, -3, -4, -2, -3, -1,
    -1, 0, -3, 0, 2, -4, -2, 1, -3, -8, 2, -2, 0, 0, -8, -1, 7, 1, 0, -1, -8, -3, -1, -1, -1, 4,
    -2, -1, -4, -2, 0, -3, -3, 0, -4, -8, 3, -3, -2, -1, -8, -3, 1, 7, -1, -1, -8, -3, -3, -1, -1, 0,
    1, 0, -1, 0, -1, -3, 0, -1, -3, -8, 0, -3, -2, 1, -8, -1, 0, -1, 5, 2, -8, -2, -4, -1, -2, 0,
    0, 0, -1, -1, -1, -2, -2, -2, -1, -8, -1, -1, -1, 0, -8, -1, -1, -1, 2, 5, -8, 0, -3, 0, -2, -1,
    -8, -8, -8, -8, -8, -8, -8, -8, -8, 1, -8, -8, -8, -8, 1, -8, -8, -8, -8, -8, 1, -8, -8, -8, -8, -8,
    0, -4, -1, -4, -3, -1, -4, -4, 4, -8, -3, 1, 1, -3, -8, -3, -3, -3, -2, 0, -8, 5, -3, -1, -1, -3,
    -3, -5, -5, -5, -3, 1, -3, -3, -3, -8, -3, -2, -1, -4, -8, -4, -1, -3, -4, -3, -8, -3, 15, -3, 2, -2,
    -1, -1, -2, -1, -1, -2, -2, -1, -1, -8, -1, -1, -1, -1, -8, -2, -1, -1, -1, 0, -8, -1, -3, -1, -1, -1,
    -2, -3, -3, -3, -2, 4, -3, 2, -1, -8, -2, -1, 0, -2, -8, -3, -1, -1, -2, -2, -8, -1, 2, -1, 8, -2,
    -1, 2, -3, 1, 5, -4, -2, 0, -3, -8, 1, -3, -1, 0, -8, -1, 4, 0, 0, -1, -8, -3, -2, -1, -2, 5
};

int BLOSUM_62_TABLE[26 * 26] = {
    4, -2, 0, -2, -1, -2, 0, -2, -1, -4, -1, -1, -1, -2, -4, -1, -1, -1, 1, 0, -4, 0, -3, 0, -2, -1,
    -2, 4, -3, 4, 1, -3, -1, 0, -3, -4, 0, -4, -3, 3, -4, -2, 0, -1, 0, -1, -4, -3, -4, -1, -3, 1,
    0, -3, 9, -3, -4, -2, -3, -3, -1, -4, -3, -1, -1, -3, -4, -3, -3, -3, -1, -1, -4, -1, -2, -2, -2, -3,
    -2, 4, -3, 6, 2, -3, -1, -1, -3, -4, -1, -4, -3, 1, -4, -1, 0, -2, 0, -1, -4, -3, -4, -1, -3, 1,
    -1, 1, -4, 2, 5, -3, -2, 0, -3, -4, 1, -3, -2, 0, -4, -1, 2, 0, 0, -1, -4, -2, -3, -1, -2, 4,
    -2, -3, -2, -3, -3, 6, -3, -1, 0, -4, -3, 0, 0, -3, -4, -4, -3, -3, -2, -2, -4, -1, 1, -1, 3, -3,
    0, -1, -3, -1, -2, -3, 6, -2, -4, -4, -2, -4, -3, 0, -4, -2, -2, -2, 0, -2, -4, -3, -2, -1, -3, -2,
    -2, 0, -3, -1, 0, -1, -2, 8, -3, -4, -1, -3, -2, 1, -4, -2, 0, 0, -1, -2, -4, -3, -2, -1, 2, 0,
    -1, -3, -1, -3, -3, 0, -4, -3, 4, -4, -3, 2, 1, -3, -4, -3, -3, -3, -2, -1, -4, 3, -3, -1, -1, -3,
    -4, -4, -4, -4, -4, -4, -4, -4, -4, 1, -4, -4, -4, -4, 1, -4, -4, -4, -4, -4, 1, -4, -4, -4, -4, -4,
    -1, 0, -3, -1, 1, -3, -2, -1, -3, -4, 5, -2, -1, 0, -4, -1, 1, 2, 0, -1, -4, -2, -3, -1, -2, 1,
    -1, -4, -1, -4, -3, 0, -4, -3, 2, -4, -2, 4, 2, -3, -4, -3, -2, -2, -2, -1, -4, 1, -2, -1, -1, -3,
    -1, -3, -1, -3, -2, 0, -3, -2, 1, -4, -1, 2, 5, -2, -4, -2, 0, -1, -1, -1, -4, 1, -1, -1, -1, -1,
    -2, 3, -3, 1, 0, -3, 0, 1, -3, -4, 0, -3, -2, 6, -4, -2, 0, 0, 1, 0, -4, -3, -4, -1, -2, 0,
    -4, -4, -4, -4, -4, -4, -4, -4, -4, 1, -4, -4, -4, -4, 1, -4, -4, -4, -4, -4, 1, -4, -4, -4, -4, -4,
    -1, -2, -3, -1, -1, -4, -2, -2, -3, -4, -1, -3, -2, -2, -4, 7, -1, -2, -1, -1, -4, -2, -4, -2, -3, -1,
    -1, 0, -3, 0, 2, -3, -2, 0, -3, -4, 1, -2, 0, 0, -4, -1, 5, 1, 0, -1, -4, -2, -2, -1, -1, 3,
    -1, -1, -3, -2, 0, -3, -2, 0, -3, -4, 2, -2, -1, 0, -4, -2, 1, 5, -1, -1, -4, -3, -3, -1, -2, 0,
    1, 0, -1, 0, 0, -2, 0, -1, -2, -4, 0, -2, -1, 1, -4, -1, 0, -1, 4, 1, -4, -2, -3, 0, -2, 0,
    0, -1, -1, -1, -1, -2, -2, -2, -1, -4, -1, -1, -1, 0, -4, -1, -1, -1, 1, 5, -4, 0, -2, 0, -2, -1,
    -4, -4, -4, -4, -4, -4, -4, -4, -4, 1, -4, -4, -4, -4, 1, -4, -4, -4, -4, -4, 1, -4, -4, -4, -4, -4,
    0, -3, -1, -3, -2, -1, -3, -3, 3, -4, -2, 1, 1, -3, -4, -2, -2, -3, -2, 0, -4, 4, -3, -1, -1, -2,
    -3, -4, -2, -4, -3, 1, -2, -2, -3, -4, -3, -2, -1, -4, -4, -4, -2, -3, -3, -2, -4, -3, 11, -2, 2, -3,
    0, -1, -2, -1, -1, -1, -1, -1, -1, -4, -1, -1, -1, -1, -4, -2, -1, -1, 0, 0, -4, -1, -2, -1, -1, -1,
    -2, -3, -2, -3, -2, 3, -3, 2, -1, -4, -2, -1, -1, -2, -4, -3, -1, -2, -2, -2, -4, -1, 2, -1, 7, -2,
    -1, 1, -3, 1, 4, -3, -2, 0, -3, -4, 1, -3, -1, 0, -4, -1, 3, 0, 0, -1, -4, -2, -3, -1, -2, 4
};

int BLOSUM_80_TABLE[26 * 26] = {
    5, -2, -1, -2, -1, -3, 0, -2, -2, -2, -1, -2, -1, -2, -6, -1, -1, -2, 1, 0, -6, 0, -3, -1, -2, -1,
    -2, 5, -4, 5, 1, -4, -1, -1, -4, -4, -1, -4, -3, 5, -6, -2, 0, -1, 0, -1, -6, -4, -5, -1, -3, 0,
    -1, -4, 9, -4, -5, -3, -4, -4, -2, -2, -4, -2, -2, -3, -6, -4, -4, -4, -2, -1, -6, -1, -3, -1, -3, -4,
    -2, 5, -4, 6, 1, -4, -2, -2, -4, -5, -1, -5, -4, 1, -6, -2, -1, -2, -1, -1, -6, -4, -6, -1, -4, 1,
    -1, 1, -5, 1, 6, -4, -3, 0, -4, -4, 1, -4, -2, -1, -6, -2, 2, -1, 0, -1, -6, -3, -4, -1, -3, 5,
    -3, -4, -3, -4, -4, 6, -4, -2, -1, 0, -4, 0, 0, -4, -6, -4, -4, -4, -3, -2, -6, -1, 0, -1, 3, -4,
    0, -1, -4, -2, -3, -4, 6, -3, -5, -5, -2, -4, -4, -1, -6, -3, -2, -3, -1, -2, -6, -4, -4, -1, -4, -3,
    -2, -1, -4, -2, 0, -2, -3, 8, -4, -4, -1, -3, -2, 0, -6, -3, 1, 0, -1, -2, -6, -4, -3, -1, 2, 0,
    -2, -4, -2, -4, -4, -1, -5, -4, 5, 3, -3, 1, 1, -4, -6, -4, -3, -3, -3, -1, -6, 3, -3, -1, -2, -4,
    -2, -4, -2, -5, -4, 0, -5, -4, 3, 3, -3, 3, 2, -4, -6, -4, -3, -3, -3, -1, -6, 2, -3, -1, -2, -3,
    -1, -1, -4, -1, 1, -4, -2, -1, -3, -3, 5, -3, -2, 0, -6, -1, 1, 2, -1, -1, -6, -3, -4, -1, -3, 1,
    -2, -4, -2, -5, -4, 0, -4, -3, 1, 3, -3, 4, 2, -4, -6, -3, -3, -3, -3, -2, -6, 1, -2, -1, -2, -3,
    -1, -3, -2, -4, -2, 0, -4, -2, 1, 2, -2, 2, 6, -3, -6, -3, 0, -2, -2, -1, -6, 1, -2, -1, -2, -1,
    -2, 5, -3, 1, -1, -4, -1, 0, -4, -4, 0, -4, -3, 6, -6, -3, 0, -1, 0, 0, -6, -4, -4, -1, -3, 0,
    -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, 1, -6, -6, -6, -6, -6, 1, -6, -6, -6, -6, -6,
    -1, -2, -4, -2, -2, -4, -3, -3, -4, -4, -1, -3, -3, -3, -6, 8, -2, -2, -1, -2, -6, -3, -5, -1, -4, -2,
    -1, 0, -4, -1, 2, -4, -2, 1, -3, -3, 1, -3, 0, 0, -6, -2, 6, 1, 0, -1, -6, -3, -3, -1, -2, 4,
    -2, -1, -4, -2, -1, -4, -3, 0, -3, -3, 2, -3, -2, -1, -6, -2, 1, 6, -1, -1, -6, -3, -4, -1, -3, 0,
    1, 0, -2, -1, 0, -3, -1, -1, -3, -3, -1, -3, -2, 0, -6, -1, 0, -1, 5, 1, -6, -2, -4, -1, -2, 0,
    0, -1, -1, -1, -1, -2, -2, -2, -1, -1, -1, -2, -1, 0, -6, -2, -1, -1, 1, 5, -6, 0, -4, -1, -2, -1,
    -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, 1, -6, -6, -6, -6, -6, 1, -6, -6, -6, -6, -6,
    0, -4, -1, -4, -3, -1, -4, -4, 3, 2, -3, 1, 1, -4, -6, -3, -3, -3, -2, 0, -6, 4, -3, -1, -2, -3,
    -3, -5, -3, -6, -4, 0, -4, -3, -3, -3, -4, -2, -2, -4, -6, -5, -3, -4, -4, -4, -6, -3, 11, -1, 2, -3,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -6, -1, -1, -1, -1, -1, -6, -1, -1, -1, -1, -1,
    -2, -3, -3, -4, -3, 3, -4, 2, -2, -2, -3, -2, -2, -3, -6, -4, -2, -3, -2, -2, -6, -2, 2, -1, 7, -3,
    -1, 0, -4, 1, 5, -4, -3, 0, -4, -3, 1, -3, -1, 0, -6, -2, 4, 0, 0, -1, -6, -3, -3, -1, -3, 5
};

int BLOSUM_90_TABLE[26 * 26] = {
    5, -2, -1, -3, -1, -3, 0, -2, -2, -6, -1, -2, -2, -2, -6, -1, -1, -2, 1, 0, -6, -1, -4, -1, -3, -1,
    -2, 4, -4, 4, 0, -4, -2, -1, -5, -6, -1, -5, -4, 4, -6, -3, -1, -2, 0, -1, -6, -4, -6, -2, -4, 0,
    -1, -4, 9, -5, -6, -3, -4, -5, -2, -6, -4, -2, -2, -4, -6, -4, -4, -5, -2, -2, -6, -2, -4, -3, -4, -5,
    -3, 4, -5, 7, 1, -5, -2, -2, -5, -6, -1, -5, -4, 1, -6, -3, -1, -3, -1, -2, -6, -5, -6, -2, -4, 0,
    -1, 0, -6, 1, 6, -5, -3, -1, -4, -6, 0, -4, -3, -1, -6, -2, 2, -1, -1, -1, -6, -3, -5, -2, -4, 4,
    -3, -4, -3, -5, -5, 7, -5, -2, -1, -6, -4, 0, -1, -4, -6, -4, -4, -4, -3, -3, -6, -2, 0, -2, 3, -4,
    0, -2, -4, -2, -3, -5, 6, -3, -5, -6, -2, -5, -4, -1, -6, -3, -3, -3, -1, -3, -6, -5, -4, -2, -5, -3,
    -2, -1, -5, -2, -1, -2, -3, 8, -4, -6, -1, -4, -3, 0, -6, -3, 1, 0, -2, -2, -6, -4, -3, -2, 1, 0,
    -2, -5, -2, -5, -4, -1, -5, -4, 5, -6, -4, 1, 1, -4, -6, -4, -4, -4, -3, -1, -6, 3, -4, -2, -2, -4,
    -6, -6, -6, -6, -6, -6, -6, -6, -6, 1, -6, -6, -6, -6, 1, -6, -6, -6, -6, -6, 1, -6, -6, -6, -6, -6,
    -1, -1, -4, -1, 0, -4, -2, -1, -4, -6, 6, -3, -2, 0, -6, -2, 1, 2, -1, -1, -6, -3, -5, -1, -3, 1,
    -2, -5, -2, -5, -4, 0, -5, -4, 1, -6, -3, 5, 2, -4, -6, -4, -3, -3, -3, -2, -6, 0, -3, -2, -2, -4,
    -2, -4, -2, -4, -3, -1, -4, -3, 1, -6, -2, 2, 7, -3, -6, -3, 0, -2, -2, -1, -6, 0, -2, -1, -2, -2,
    -2, 4, -4, 1, -1, -4, -1, 0, -4, -6, 0, -4, -3, 7, -6, -3, 0, -1, 0, 0, -6, -4, -5, -2, -3, -1,
    -6, -6, -6, -6, -6, -6, -6, -6, -6, 1, -6, -6, -6, -6, 1, -6, -6, -6, -6, -6, 1, -6, -6, -6, -6, -6,
    -1, -3, -4, -3, -2, -4, -3, -3, -4, -6, -2, -4, -3, -3, -6, 8, -2, -3, -2, -2, -6, -3, -5, -2, -4, -2,
    -1, -1, -4, -1, 2, -4, -3, 1, -4, -6, 1, -3, 0, 0, -6, -2, 7, 1, -1, -1, -6, -3, -3, -1, -3, 4,
    -2, -2, -5, -3, -1, -4, -3, 0, -4, -6, 2, -3, -2, -1, -6, -3, 1, 6, -1, -2, -6, -3, -4, -2, -3, 0,
    1, 0, -2, -1, -1, -3, -1, -2, -3, -6, -1, -3, -2, 0, -6, -2, -1, -1, 5, 1, -6, -2, -4, -1, -3, -1,
    0, -1, -2, -2, -1, -3, -3, -2, -1, -6, -1, -2, -1, 0, -6, -2, -1, -2, 1, 6, -6, -1, -4, -1, -2, -1,
    -6, -6, -6, -6, -6, -6, -6, -6, -6, 1, -6, -6, -6, -6, 1, -6, -6, -6, -6, -6, 1, -6, -6, -6, -6, -6,
    -1, -4, -2, -5, -3, -2, -5, -4, 3, -6, -3, 0, 0, -4, -6, -3, -3, -3, -2, -1, -6, 5, -3, -2, -3, -3,
    -4, -6, -4, -6, -5, 0, -4, -3, -4, -6, -5, -3, -2, -5, -6, -5, -3, -4, -4, -4, -6, -3, 11, -3, 2, -4,
    -1, -2, -3, -2, -2, -2, -2, -2, -2, -6, -1, -2, -1, -2, -6, -2, -1, -2, -1, -1, -6, -2, -3, -2, -2, -1,
    -3, -4, -4, -4, -4, 3, -5, 1, -2, -6, -3, -2, -2, -3, -6, -4, -3, -3, -3, -2, -6, -3, 2, -2, 8, -3,
    -1, 0, -5, 0, 4, -4, -3, 0, -4, -6, 1, -4, -2, -1, -6, -2, 4, 0, -1, -1, -6, -3, -4, -1, -3, 4
};

int PAM_30_TABLE[26 * 26] = {
    6, -3, -6, -3, -2, -8, -2, -7, -5, -6, -7, -6, -5, -4, -17, -2, -4, -7, 0, -1, -17, -2, -13, -1, -8, -3,
    -3, 6, -12, 6, 1, -10, -3, -1, -6, -8, -2, -9, -10, 6, -17, -7, -3, -7, -1, -3, -17, -8, -10, -1, -6, 0,
    -6, -12, 10, -14, -14, -13, -9, -7, -6, -9, -14, -15, -13, -11, -17, -8, -14, -8, -3, -8, -17, -6, -15, -1, -4, -14,
    -3, 6, -14, 8, 2, -15, -3, -4, -7, -10, -4, -12, -11, 2, -17, -8, -2, -10, -4, -5, -17, -8, -15, -1, -11, 1,
    -2, 1, -14, 2, 8, -14, -4, -5, -5, -7, -4, -9, -7, -2, -17, -5, 1, -9, -4, -6, -17, -6, -17, -1, -8, 6,
    -8, -10, -13, -15, -14, 9, -9, -6, -2, -2, -14, -3, -4, -9, -17, -10, -13, -9, -6, -9, -17, -8, -4, -1, 2, -13,
    -2, -3, -9, -3, -4, -9, 6, -9, -11, -10, -7, -10, -8, -3, -17, -6, -7, -9, -2, -6, -17, -5, -15, -1, -14, -5,
    -7, -1, -7, -4, -5, -6, -9, 9, -9, -7, -6, -6, -10, 0, -17, -4, 1, -2, -6, -7, -17, -6, -7, -1, -3, -1,
    -5, -6, -6, -7, -5, -2, -11, -9, 8, 5, -6, -1, -1, -5, -17, -8, -8, -5, -7, -2, -17, 2, -14, -1, -6, -6,
    -6, -8, -9, -10, -7, -2, -10, -7, 5, 6, -7, 6, 0, -6, -17, -7, -5, -7, -8, -5, -17, 0, -7, -1, -7, -6,
    -7, -2, -14, -4, -4, -14, -7, -6, -6, -7, 7, -8, -2, -1, -17, -6, -3, 0, -4, -3, -17, -9, -12, -1, -9, -4,
    -6, -9, -15, -12, -9, -3, -10, -6, -1, 6, -8, 7, 1, -7, -17, -7, -5, -8, -8, -7, -17, -2, -6, -1, -7, -7,
    -5, -10, -13, -11, -7, -4, -8, -10, -1, 0, -2, 1, 11, -9, -17, -8, -4, -4, -5, -4, -17, -1, -13, -1, -11, -5,
    -4, 6, -11, 2, -2, -9, -3, 0, -5, -6, -1, -7, -9, 8, -17, -6, -3, -6, 0, -2, -17, -8, -8, -1, -4, -3,
    -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, 1, -17, -17, -17, -17, -17, 1, -17, -17, -17, -17, -17,
    -2, -7, -8, -8, -5, -10, -6, -4, -8, -7, -6, -7, -8, -6, -17, 8, -3, -4, -2, -4, -17, -6, -14, -1, -13, -4,
    -4, -3, -14, -2, 1, -13, -7, 1, -8, -5, -3, -5, -4, -3, -17, -3, 8, -2, -5, -5, -17, -7, -13, -1, -12, 6,
    -7, -7, -8, -10, -9, -9, -9, -2, -5, -7, 0, -8, -4, -6, -17, -4, -2, 8, -3, -6, -17, -8, -2, -1, -10, -4,
    0, -1, -3, -4, -4, -6, -2, -6, -7, -8, -4, -8, -5, 0, -17, -2, -5, -3, 6, 0, -17, -6, -5, -1, -7, -5,
    -1, -3, -8, -5, -6, -9, -6, -7, -2, -5, -3, -7, -4, -2, -17, -4, -5, -6, 0, 7, -17, -3, -13, -1, -6, -6,
    -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, 1, -17, -17, -17, -17, -17, 1, -17, -17, -17, -17, -17,
    -2, -8, -6, -8, -6, -8, -5, -6, 2, 0, -9, -2, -1, -8, -17, -6, -7, -8, -6, -3, -17, 7, -15, -1, -7, -6,
    -13, -10, -15, -15, -17, -4, -15, -7, -14, -7, -12, -6, -13, -8, -17, -14, -13, -2, -5, -13, -17, -15, 13, -1, -5, -14,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -17, -1, -1, -1, -1, -1, -17, -1, -1, -1, -1, -1,
    -8, -6, -4, -11, -8, 2, -14, -3, -6, -7, -9, -7, -11, -4, -17, -13, -12, -10, -7, -6, -17, -7, -5, -1, 10, -9,
    -3, 0, -14, 1, 6, -13, -5, -1, -6, -6, -4, -7, -5, -3, -17, -4, 6, -4, -5, -6, -17, -6, -14, -1, -9, 6
};

int PAM_70_TABLE[26 * 26] = {
    5, -1, -4, -1, -1, -6, 0, -4, -2, -3, -4, -4, -3, -2, -11, 0, -2, -4, 1, 1, -11, -1, -9, -1, -5, -1,
    -1, 5, -8, 5, 2, -7, -1, 0, -4, -5, -1, -6, -6, 5, -11, -4, -1, -4, 0, -1, -11, -5, -7, -1, -4, 1,
    -4, -8, 9, -9, -9, -8, -6, -5, -4, -7, -9, -10, -9, -7, -11, -5, -9, -5, -1, -5, -11, -4, -11, -1, -2, -9,
    -1, 5, -9, 6, 3, -10, -1, -1, -5, -7, -2, -8, -7, 3, -11, -4, 0, -6, -1, -2, -11, -5, -10, -1, -7, 2,
    -1, 2, -9, 3, 6, -9, -2, -2, -4, -5, -2, -6, -4, 0, -11, -3, 2, -5, -2, -3, -11, -4, -11, -1, -6, 5,
    -6, -7, -8, -10, -9, 8, -7, -4, 0, -1, -9, -1, -2, -6, -11, -7, -9, -7, -4, -6, -11, -5, -2, -1, 4, -9,
    0, -1, -6, -1, -2, -7, 6, -6, -6, -7, -5, -7, -6, -1, -11, -3, -4, -6, 0, -3, -11, -3, -10, -1, -9, -3,
    -4, 0, -5, -1, -2, -4, -6, 8, -6, -4, -3, -4, -6, 1, -11, -2, 2, 0, -3, -4, -11, -4, -5, -1, -1, 1,
    -2, -4, -4, -5, -4, 0, -6, -6, 7, 4, -4, 1, 1, -3, -11, -5, -5, -3, -4, -1, -11, 3, -9, -1, -4, -4,
    -3, -5, -7, -7, -5, -1, -7, -4, 4, 5, -5, 5, 2, -4, -11, -5, -3, -5, -5, -3, -11, 1, -5, -1, -4, -4,
    -4, -1, -9, -2, -2, -9, -5, -3, -4, -5, 6, -5, 0, 0, -11, -4, -1, 2, -2, -1, -11, -6, -7, -1, -7, -2,
    -4, -6, -10, -8, -6, -1, -7, -4, 1, 5, -5, 6, 2, -5, -11, -5, -3, -6, -6, -4, -11, 0, -4, -1, -4, -4,
    -3, -6, -9, -7, -4, -2, -6, -6, 1, 2, 0, 2, 10, -5, -11, -5, -2, -2, -3, -2, -11, 0, -8, -1, -7, -3,
    -2, 5, -7, 3, 0, -6, -1, 1, -3, -4, 0, -5, -5, 6, -11, -3, -1, -3, 1, 0, -11, -5, -6, -1, -3, -1,
    -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, 1, -11, -11, -11, -11, -11, 1, -11, -11, -11, -11, -11,
    0, -4, -5, -4, -3, -7, -3, -2, -5, -5, -4, -5, -5, -3, -11, 7, -1, -2, 0, -2, -11, -3, -9, -1, -9, -2,
    -2, -1, -9, 0, 2, -9, -4, 2, -5, -3, -1, -3, -2, -1, -11, -1, 7, 0, -3, -3, -11, -4, -8, -1, -8, 5,
    -4, -4, -5, -6, -5, -7, -6, 0, -3, -5, 2, -6, -2, -3, -11, -2, 0, 8, -1, -4, -11, -5, 0, -1, -7, -2,
    1, 0, -1, -1, -2, -4, 0, -3, -4, -5, -2, -6, -3, 1, -11, 0, -3, -1, 5, 2, -11, -3, -3, -1, -5, -2,
    1, -1, -5, -2, -3, -6, -3, -4, -1, -3, -1, -4, -2, 0, -11, -2, -3, -4, 2, 6, -11, -1, -8, -1, -4, -3,
    -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, 1, -11, -11, -11, -11, -11, 1, -11, -11, -11, -11, -11,
    -1, -5, -4, -5, -4, -5, -3, -4, 3, 1, -6, 0, 0, -5, -11, -3, -4, -5, -3, -1, -11, 6, -10, -1, -5, -4,
    -9, -7, -11, -10, -11, -2, -10, -5, -9, -5, -7, -4, -8, -6, -11, -9, -8, 0, -3, -8, -11, -10, 13, -1, -3, -10,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -11, -1, -1, -1, -1, -1, -11, -1, -1, -1, -1, -1,
    -5, -4, -2, -7, -6, 4, -9, -1, -4, -4, -7, -4, -7, -3, -11, -9, -8, -7, -5, -4, -11, -5, -3, -1, 9, -7,
    -1, 1, -9, 2, 5, -9, -3, 1, -4, -4, -2, -4, -3, -1, -11, -2, 5, -2, -2, -3, -11, -4, -10, -1, -7, 5
};

int PAM_250_TABLE[26 * 26] = {
    2, 0, -2, 0, 0, -3, 1, -1, -1, -8, -1, -2, -1, 0, -8, 1, 0, -2, 1, 1, -8, 0, -6, 0, -3, 0,
    0, 3, -4, 3, 3, -4, 0, 1, -2, -8, 1, -3, -2, 2, -8, -1, 1, -1, 0, 0, -8, -2, -5, -1, -3, 2,
    -2, -4, 12, -5, -5, -4, -3, -3, -2, -8, -5, -6, -5, -4, -8, -3, -5, -4, 0, -2, -8, -2, -8, -3, 0, -5,
    0, 3, -5, 4, 3, -6, 1, 1, -2, -8, 0, -4, -3, 2, -8, -1, 2, -1, 0, 0, -8, -2, -7, -1, -4, 3,
    0, 3, -5, 3, 4, -5, 0, 1, -2, -8, 0, -3, -2, 1, -8, -1, 2, -1, 0, 0, -8, -2, -7, -1, -4, 3,
    -3, -4, -4, -6, -5, 9, -5, -2, 1, -8, -5, 2, 0, -3, -8, -5, -5, -4, -3, -3, -8, -1, 0, -2, 7, -5,
    1, 0, -3, 1, 0, -5, 5, -2, -3, -8, -2, -4, -3, 0, -8, 0, -1, -3, 1, 0, -8, -1, -7, -1, -5, 0,
    -1, 1, -3, 1, 1, -2, -2, 6, -2, -8, 0, -2, -2, 2, -8, 0, 3, 2, -1, -1, -8, -2, -3, -1, 0, 2,
    -1, -2, -2, -2, -2, 1, -3, -2, 5, -8, -2, 2, 2, -2, -8, -2, -2, -2, -1, 0, -8, 4, -5, -1, -1, -2,
    -8, -8, -8, -8, -8, -8, -8, -8, -8, 1, -8, -8, -8, -8, 1, -8, -8, -8, -8, -8, 1, -8, -8, -8, -8, -8,
    -1, 1, -5, 0, 0, -5, -2, 0, -2, -8, 5, -3, 0, 1, -8, -1, 1, 3, 0, 0, -8, -2, -3, -1, -4, 0,
    -2, -3, -6, -4, -3, 2, -4, -2, 2, -8, -3, 6, 4, -3, -8, -3, -2, -3, -3, -2, -8, 2, -2, -1, -1, -3,
    -1, -2, -5, -3, -2, 0, -3, -2, 2, -8, 0, 4, 6, -2, -8, -2, -1, 0, -2, -1, -8, 2, -4, -1, -2, -2,
    0, 2, -4, 2, 1, -3, 0, 2, -2, -8, 1, -3, -2, 2, -8, 0, 1, 0, 1, 0, -8, -2, -4, 0, -2, 1,
    -8, -8, -8, -8, -8, -8, -8, -8, -8, 1, -8, -8, -8, -8, 1, -8, -8, -8, -8, -8, 1, -8, -8, -8, -8, -8,
    1, -1, -3, -1, -1, -5, 0, 0, -2, -8, -1, -3, -2, 0, -8, 6, 0, 0, 1, 0, -8, -1, -6, -1, -5, 0,
    0, 1, -5, 2, 2, -5, -1, 3, -2, -8, 1, -2, -1, 1, -8, 0, 4, 1, -1, -1, -8, -2, -5, -1, -4, 3,
    -2, -1, -4, -1, -1, -4, -3, 2, -2, -8, 3, -3, 0, 0, -8, 0, 1, 6, 0, -1, -8, -2, 2, -1, -4, 0,
    1, 0, 0, 0, 0, -3, 1, -1, -1, -8, 0, -3, -2, 1, -8, 1, -1, 0, 2, 1, -8, -1, -2, 0, -3, 0,
    1, 0, -2, 0, 0, -3, 0, -1, 0, -8, 0, -2, -1, 0, -8, 0, -1, -1, 1, 3, -8, 0, -5, 0, -3, -1,
    -8, -8, -8, -8, -8, -8, -8, -8, -8, 1, -8, -8, -8, -8, 1, -8, -8, -8, -8, -8, 1, -8, -8, -8, -8, -8,
    0, -2, -2, -2, -2, -1, -1, -2, 4, -8, -2, 2, 2, -2, -8, -1, -2, -2, -1, 0, -8, 4, -6, -1, -2, -2,
    -6, -5, -8, -7, -7, 0, -7, -3, -5, -8, -3, -2, -4, -4, -8, -6, -5, 2, -2, -5, -8, -6, 17, -4, 0, -6,
    0, -1, -3, -1, -1, -2, -1, -1, -1, -8, -1, -1, -1, 0, -8, -1, -1, -1, 0, 0, -8, -1, -4, -1, -2, -1,
    -3, -3, 0, -4, -4, 7, -5, 0, -1, -8, -4, -1, -2, -2, -8, -5, -4, -4, -3, -3, -8, -2, 0, -2, 10, -4,
    0, 2, -5, 3, 3, -5, 0, 2, -2, -8, 0, -3, -2, 1, -8, 0, 3, 0, 0, -1, -8, -2, -6, -1, -4, 3
};

int EDNA_FULL_TABLE[26 * 26] = {
    5, -10, -4, -10, -10, -10, -4, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -4, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -4, -10, 5, -10, -10, -10, -4, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -4, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -4, -10, -4, -10, -10, -10, 5, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -4, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -4, -10, -4, -10, -10, -10, -4, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, 5, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10,
    -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10
};
