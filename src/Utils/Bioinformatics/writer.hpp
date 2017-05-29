
#pragma once

#include <cstdio>
#include <vector>
#include <string>
#include <fstream>
#include <memory>
#include "../../Tachyon/Base.h"
#include "../Ref.hpp"

class DatabaseElement;
class ScoreMatrix;
class Alignment;
class Writer;

using ChainSet = std::vector<DatabaseElement>;
using AlignmentSet = std::vector<Alignment>;

enum class OutputType {
	kBm0, // BLAST m0
	kBm8, // BLAST m8
	kBm9, // BLAST m9
};

class Writer {
public:
	~Writer();
	Writer(string & output_file, OutputType format, ScoreMatrix& scorer);
	Writer(FILE* output_file, OutputType format, ScoreMatrix& scorer);
	void write_alignments(const AlignmentSet& alignments, const ChainSet& queries,
	                      const Ref<Base>& database);


private:

	Writer(const Writer&) = delete;
	const Writer& operator=(const Writer&) = delete;

	void write_bm0(const AlignmentSet& alignments, const ChainSet& queries,
	               const Ref<Base>& database);

	void write_bm8(const AlignmentSet& alignments, const ChainSet& queries,
				   const Ref<Base>& database);

	void write_bm9(const AlignmentSet& alignments, const ChainSet& queries,
				   const Ref<Base>& database);

	FILE* output_file_;
	OutputType format_;
	ScoreMatrix scorer_;
};
