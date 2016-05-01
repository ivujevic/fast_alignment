
#ifndef SCORER_H
#define	SCORER_H

typedef struct Scorer Scorer;

Scorer* scorerCreate(const char* name, int* scores, char maxCode, int gapOpen, int gapExtend);
void scorerCreateMatrix(Scorer** scorer, char* name, int gapOpen, int gapExtend);

const char* scorerGetName(Scorer* scorer);
int scorerGetGapOpen(Scorer* scorer);
int scorerGetGapExtend(Scorer* scorer);
#endif	/* SCORER_H */

