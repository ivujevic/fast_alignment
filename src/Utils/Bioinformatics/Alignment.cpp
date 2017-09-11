//
// Created by vujevic on 27.04.16..
//

#include <iostream>
#include "Alignment.hpp"

bool compareAlignment(const Alignment& left,
                      const Alignment& right) {

	if (left.evalue() < right.evalue()) return true;
	if (left.evalue() == right.evalue() && left.score() > right.score()) return true;
	return false;
}

Alignment::Alignment() {

}

Alignment::Alignment(int32_t score, double evalue, uint32_t queryId, const DatabaseElement* target)
	: score_(score), evalue_(evalue), query_id_(queryId), target_(target) {

}

void Alignment::update(uint32_t query_begin, uint32_t query_end, uint32_t target_begin,
                       uint32_t target_end, const unsigned char* alignment, uint32_t length) {

	query_begin_ = query_begin;
	query_end_ = query_end;
	target_begin_ = target_begin;
	target_end_ = target_end;

	alignmentLen_ = length;

	alignment_.clear();
    alignment_.resize(length);

	for (int  i = 0; i < length; i++) {
		alignment_[i] = alignment[i];
	}
}

