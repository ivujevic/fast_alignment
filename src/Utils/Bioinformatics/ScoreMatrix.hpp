#pragma once


#include <vector>
#include <string>


class ScoreMatrix;

enum class ScoreMatrixType {
    kBlosum45 = 0,
    kBlosum50,
    kBlosum62,
    kBlosum80,
    kBlosum90,
    kPam30,
    kPam70,
    kPam250
};

using namespace std;

class ScoreMatrix {
 private:
    vector<unsigned char> alphabet; //!< letters in same order as columns/rows in matrix
    vector<int> matrix; //!< All rows of matrix concatenated. Has length: alphabetLength * alphabetLength

 public:
    ScoreMatrix(ScoreMatrixType type, int gap_open, int gap_extend);
    ScoreMatrix();
    ScoreMatrix(vector<unsigned char> alphabet, vector<int> matrix);
    /**
     * Format of matrix file:
     *  - first line contains letters from alphabet naming columns separated with spaces.
     *  - each next line is one row of matrix (integers separated with spaces).
     */
    ScoreMatrix(const char* filepath);

    int getAlphabetLength();
    unsigned char* getAlphabet();
    int* getMatrix();

    static ScoreMatrix getBlosum50();
    static ScoreMatrix getBlosum60();
    static ScoreMatrix getBlosum62();
    int gap_open(){return gap_open_;};
    int gap_extend(){return gap_extend_;};

    int score(int row, int column);
    std::string scorerName() const;

 private:
    ScoreMatrixType type_;
    int gap_open_;
    int gap_extend_;
    int num_columns_ = 24;
    static vector<unsigned char> getBlosumAlphabet();
};

