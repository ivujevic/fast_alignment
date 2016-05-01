#include <stdlib.h>
#include <string.h>

#include "Scorer.h"
#include "Constants.h"

#define SCORERS_LEN (sizeof(scorers) / sizeof(ScorerEntry))

struct Scorer {

    char* name;
    int nameLen;
    
    int gapOpen;
    int gapExtend;
    
    int* table;


};


typedef struct ScorerEntry {
    const char* name;
    int (*table)[26 * 26];
} ScorerEntry;

// to register a scorer just add his name and corresponding table to this array
static ScorerEntry scorers[] = {
    { "BLOSUM_62", &BLOSUM_62_TABLE }, // default one
    { "BLOSUM_45", &BLOSUM_45_TABLE },
    { "BLOSUM_50", &BLOSUM_50_TABLE },
    { "BLOSUM_80", &BLOSUM_80_TABLE },
    { "BLOSUM_90", &BLOSUM_90_TABLE },
    { "PAM_30", &PAM_30_TABLE },
    { "PAM_70", &PAM_70_TABLE },
    { "PAM_250", &PAM_250_TABLE },
    { "EDNA_FULL", &EDNA_FULL_TABLE }
};


static const char CODER[] = {
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  26,  27, 
     28,  29,  30,  31,  32,  33,  34,  35,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,   0,   1,   2,   3,   4, 
      5,   6,   7,   8,   9,  10,  11,  12,  13,  14, 
     15,  16,  17,  18,  19,  20,  21,  22,  23,  24, 
     25,  -1,  -1,  -1,  -1,  -1,  -1,   0,   1,   2, 
      3,   4,   5,   6,   7,   8,   9,  10,  11,  12, 
     13,  14,  15,  16,  17,  18,  19,  20,  21,  22, 
     23,  24,  25,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1
};

static const char DECODER[] = {
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 
    'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 
    'U', 'V', 'W', 'X', 'Y', 'Z', '0', '1', '2', '3', 
    '4', '5', '6', '7', '8', '9',  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
     -1,  -1,  -1,  -1,  -1
};

Scorer* scorerCreate(const char* name, int* scores, char maxCode, int gapOpen, int gapExtend) {
    for (int i = 0; i < maxCode; ++i) {
        for (int j = i + 1; j < maxCode; ++j) {
            int a = scores[i * maxCode + j];
            int b = scores[j * maxCode + i];
        }
    }
    
    Scorer* scorer =(Scorer*) malloc(sizeof(struct Scorer));
    scorer->nameLen = strlen(name) + 1;
    scorer->name = (char*) malloc(scorer->nameLen * sizeof(char));
    scorer->name[scorer->nameLen - 1] = '\0';
    memcpy(scorer->name, name, (scorer->nameLen - 1) * sizeof(char));
    
     scorer->gapOpen = gapOpen;
    scorer->gapExtend = gapExtend;
    
    size_t tableSize = maxCode * maxCode * sizeof(int);    
    scorer->table = (int*) malloc(tableSize);
    memcpy(scorer->table, scores, tableSize);
    
    return scorer;
}


void scorerCreateMatrix(Scorer** scorer, char* name, int gapOpen, int gapExtend) {
    
    int index = -1;
    for(int i=0; i <SCORERS_LEN ; i++) {
        if(strcmp(name,scorers[i].name) == 0) {
            
            index = i;
            break;
        }
    }
    ScorerEntry* entry = &(scorers[index]);
    *scorer = scorerCreate(entry->name,*(entry->table),26,gapOpen,gapExtend);
}

int scorerGetGapExtend(Scorer* scorer) {
    return scorer->gapExtend;
}

 int scorerGetGapOpen(Scorer* scorer) {
    return scorer->gapOpen;
}
 
 const char* scorerGetName(Scorer* scorer) {
    return scorer->name;
}
