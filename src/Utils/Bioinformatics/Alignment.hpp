#pragma once

#include <stdint.h>
#include <string.h>
#include <vector>
#include <memory>
#include "../../Tachyon/DatabaseElement.h"

enum class OutputType;
class ScoreMatrix;
class EValue;
class Alignment;

enum class AlignmentType {
	kNW = 0, // global alignment (Needleman-Wunsch)
	kHW, // semi-global alignemnt (query - _target_)
	kOV, // semi-global alignment (_query_ - _target_)
	kSW // local alignment (Smith-Waterman)
};


using AlignmentSet = std::vector<Alignment>;

bool compareAlignment(const Alignment& left, const Alignment& right);


class Alignment {
public:


	Alignment();
	Alignment(int32_t score, double evalue, uint32_t queryId, const DatabaseElement* target);

	~Alignment() = default;

	int32_t score() const {
		return score_;
	}

	double evalue() const {
		return evalue_;
	}

	uint32_t query_id() const {
		return query_id_;
	}

	uint32_t query_begin() const {
		return query_begin_;
	}

	uint32_t query_end() const {
		return query_end_;
	}

	const DatabaseElement* target() const {
		return target_;
	}

	uint32_t target_begin() const {
		return target_begin_;
	}

	uint32_t target_end() const {
		return target_end_;
	}

	uint32_t getAlignmentLen() const {
		return alignmentLen_;
	}

	const std::string& alignment() const {
		return alignment_;
	}

	void update(uint32_t query_begin, uint32_t query_end, uint32_t target_begin,
	            uint32_t target_end, const unsigned char* alignment, uint32_t length);


protected:


	int32_t score_;
	double evalue_;

	uint32_t query_id_;
    const DatabaseElement* target_;

	uint32_t query_begin_;
	uint32_t query_end_;
	uint32_t target_begin_;
	uint32_t target_end_;

	uint32_t alignmentLen_;

	std::string alignment_;
};
