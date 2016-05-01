#ifndef CONSTANTS_H
#define	CONSTANTS_H

// Constants for TCP connection
#define ip_adress "127.0.0.1" //localhost
#define server_port 9341 // Tachyon listens on this port

// Constants for peptides size and minimum hits
#define proteinsWordLength 5
#define proteinsListSize 7
#define proteinsMinHits 2

// Constants for nucleotides
#define nucleotidesWordLength 15
#define nucleotidesListSize 5
#define nucleotidesMinHits 3

// Constants for SEG
#define windowSize 12
#define lowerEntropy 2.2d
#define upperEntropy 2.5d

extern int BLOSUM_45_TABLE[26 * 26];
extern int BLOSUM_50_TABLE[26 * 26];
extern int BLOSUM_62_TABLE[26 * 26];
extern int BLOSUM_80_TABLE[26 * 26];
extern int BLOSUM_90_TABLE[26 * 26];

extern int PAM_30_TABLE[26 * 26];
extern int PAM_70_TABLE[26 * 26];
extern int PAM_250_TABLE[26 * 26];

extern int EDNA_FULL_TABLE[26 * 26];
//extern char codonTable[4][4][4];
#endif

