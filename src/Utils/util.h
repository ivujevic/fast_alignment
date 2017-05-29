#pragma once

#include<vector>
#include<string>
#include <set>
#include <bits/unique_ptr.h>
#include "../Tachyon/DatabaseElement.h"
#include "../Tachyon/Base.h"

using ChainSet = std::vector<DatabaseElement>;

extern void convertDNAToAA(const std::string& dnaSeq, std::vector<std::string>& aaSeqs);
/**
 * This function reads fasta file and put all elements in vector of pairs(id,sequence).
 * Return number of readed queries
 */
extern int readFastaFile(std::string file, std::vector<std::pair<std::string,std::string>>& elements);
extern uint64_t readFastaFile(const char* file, ChainSet & elements, uint64_t database_size);

extern std::vector<std::string> split(std::string& text, char delim);

extern void get_codes(const std::string& sequence, std::vector<long>& results,int kmer_len);
extern void get_codes(const std::string& sequence, std::vector<std::pair<long,long>>& results,int kmer_len);
typedef enum {makedb, blastp, blastx, server} Command;

extern Command strToCommand(const std::string &str);
