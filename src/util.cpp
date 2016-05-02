#include<fstream>
#include <sstream>


#include "util.h"
#include "Constants.h"

using namespace std;
constexpr uint32_t kBufferSize = 1024 * 1024;
constexpr uint32_t kArraySize = 65000;

std::unordered_map<std::string,char>codonTable = {
		{"TTT",'F'},{"TTC",'F'},
		{"TTA",'L'},{"TTG",'L'},{"CTT",'L'},{"CTC",'L'},{"CTA",'L'},{"CTG",'L'},
		{"ATT",'I'},{"ATC",'I'},{"ATA",'I'},
		{"ATG",'M'},
		{"GTT",'V'},{"GTC",'V'},{"GTA",'V'},{"GTG",'V'},
		{"TCT",'S'},{"TCC",'S'},{"TCA",'S'},{"TCG",'S'},
		{"CCT",'P'},{"CCC",'P'},{"CCA",'P'},{"CCG",'P'},
		{"ACT",'T'},{"ACC",'T'},{"ACA",'T'},{"ACG",'T'},
		{"GCG",'A'},{"GCA",'A'},{"GCC",'A'},{"GCT",'A'},
		{"TAT",'Y'},{"TAC",'Y'},
		{"CAT",'H'},{"CAC",'H'},
		{"CAA",'Q'},{"CAG",'Q'},
		{"AAT",'N'},{"AAC",'N'},
		{"AAA",'K'},{"AAG",'K'},
		{"GAT",'D'},{"GAC",'D'},
		{"GAA",'E'},{"GAG",'E'},
		{"TGT",'C'},{"TGC",'C'},
		{"TGG",'W'},
		{"CGT",'R'},{"CGC",'R'},{"CGA",'R'},{"CGG",'R'},{"AGA",'R'},{"AGG",'R'},
		{"AGT",'S'},{"AGC",'S'},
		{"GGT",'G'},{"GGC",'G'},{"GGA",'G'},{"GGG",'G'},
		{"TAA",'*'},{"TAG",'*'},{"TGA",'*'}
};


void getComplementChain(const std::string& chain,int len, std::string& complementChain);

int readFastaFile(std::string file, std::vector< std::pair< std::string,std::string > >& elems) {
	std::ifstream queryFile(file.c_str());
	std::string id = "";
	std::string sequence = "";

	std::string line = "";
	int numberOfQueries = 0;
	while(getline(queryFile,line)){
		if(line == "") continue;
		if(line[0] == '>'){
			if(sequence != ""){
				numberOfQueries++;
				std::pair<std::string,std::string> p = make_pair(id,sequence);
				elems.push_back(p);
				sequence = "";
			}
			id = line.substr(1);
		}else{
			sequence += line;
		}
	}

	//add last element from file
	if(sequence !="" && id !=""){
		numberOfQueries++;
		std::pair<std::string,std::string>p = make_pair(id,sequence);
		elems.push_back(p);
	}
	return numberOfQueries;
}



uint64_t readFastaFile(const char* file, ChainSet & elements,uint64_t databaseSize) {

	/* Code taken from SW# (author: Matija Korpar) */

	size_t bytes_read = 0;
	size_t bytes_over = 0;

	auto input_file = fopen(file,"r");
	if(!input_file) {
		return 0;
	}
	bool is_name = true;
	bool is_end = feof(input_file);

	char name[kArraySize];
	uint32_t name_length = 0;

	char data[kArraySize];
	uint32_t data_length = 0;

	std::vector<char> buffer_(kBufferSize,'0');
	long num_chains_read_ = 0;
	while (!is_end) {

		uint32_t read = fread(buffer_.data(), sizeof(char), kBufferSize, input_file);
		is_end = feof(input_file);

		bytes_read += read;

		for (uint32_t i = 0; i < read; ++i) {

			auto c = buffer_[i];

			if (!is_name && (c == '>' || (is_end && i == read - 1))) {

				bytes_over = 0;
				is_name = true;

				elements.push_back(DatabaseElement(num_chains_read_++, name, name_length,data, data_length));

				databaseSize +=elements.back().sequence_len();
				name_length = data_length = 0;
			}

			if (is_name) {
				if (c == '\n') {
					is_name = false;
				} else if (name_length == kArraySize)  {
					continue;
				} else if (!(name_length == 0 && (c == '>' || isspace(c)))) {
					if (c != '\r') {
						name[name_length++] = c;
					}
				}
			} else {
				data[data_length++] = c;
			}

			bytes_over++;
		}
	}

	return databaseSize;
}

char getComplementNucleotide(char nucleotide) {
	switch(nucleotide){
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'T':
			return 'A';
		case 'A':
			return 'T';
		default:
			return 'N';
	}
}
void convertDNAToAA(const std::string& dnaSeq, std::vector<std::string>& aaSeqs) {
	//5' to 3'
	int len = dnaSeq.length();
	std::string complementChain = "";
	getComplementChain(dnaSeq,len,complementChain);


	for(int i = 0; i < 3; i++){
		std::string aaSeq ="";
		std:: string aaSeqCompl = "";
		for(int j = i; j <len- 2; j+=3){
			std::string codon = dnaSeq.substr(j,3);
			std::string codonCompl = complementChain.substr(j,3);
			aaSeq += codonTable[codon];
			aaSeqCompl += codonTable[codonCompl];
		}
		aaSeqs.push_back(aaSeq);
		aaSeqs.push_back(aaSeqCompl);
	}
}

void getComplementChain(const std::string& chain,int len, std::string& complementChain) {
	for(int i = len - 1; i >=0; i--) {
		char c = chain[i];
		complementChain += getComplementNucleotide(c);
	}
}


std::vector<std::string> split(std::string& text, char delim)
{
	std::vector<std::string> tokens;
	std::stringstream ss(text);
	std::string token;
	while(getline(ss,token,delim))
	{
		tokens.push_back(token);
	}
	return tokens;
}

void get_codes(const std::string& sequence, std::vector<long>& results) {
	int len = sequence.length();
	for(int i = 0; i < len - 4; i++) {
		long code = 0;
		bool isValid = true;
		string s = sequence.substr(i,5);

		for(int j = 0; j < 5; j++) {
			char c = sequence[i+j];
			code = code << 5;
			if(c == '*' || c == 'X') {
				isValid = false;
				break;
			}
			code = code | ( toupper(c) - 'A');
		}
		if(isValid) results.push_back(code);
	}
}

void get_codes(const std::string& sequence, std::set<long>& results) {
	int len = sequence.length();
	for(int i = 0; i < len - 4; i++) {
		long code = 0;
		bool isValid = true;
		for(int j = 0; j < 5; j++) {
			char c = sequence[i+j];
			code = code << 5;
			if(c == '*' || c == 'X') {
				isValid = false;
				break;
			}
			code = code | ( toupper(c) - 'A');
		}
		if(isValid) results.insert(code);
	}
}

void get_codes(const string& sequence, vector<pair<long,long>>& results) {
	int len = sequence.length();
	set<long> tempSet;
	for(int i = 0; i < len - 4; i++) {
		long code = 0;
		bool isValid = true;
		for(int j = 0; j < 5; j++) {
			char c = sequence[i+j];
			code = code << 5;
			if(c == '*' || c == 'X') {
				isValid = false;
				break;
			}
			code = code | ( toupper(c) - 'A');
		}
		if(tempSet.find(code) == tempSet.end() && isValid) {
			auto pt = make_pair(code,i);
			results.push_back(pt);
			tempSet.insert(code);
		}
	}
}