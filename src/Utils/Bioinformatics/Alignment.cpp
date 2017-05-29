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


void Alignment::update(uint32_t query_begin, uint32_t query_end, uint32_t target_begin,
                       uint32_t target_end, const unsigned char* alignment, uint32_t length) {

	query_begin_ = query_begin;
	query_end_ = query_end;
	target_begin_ = target_begin;
	target_end_ = target_end;

	alignmentLen_ = length;

//    alignment_ = new char[length + 1];
//
//    for (int i = 0; i < length; i++) {
//        alignment_[i] = alignment[i];
//    }
//	alignment_[length] = '\0';

	alignment_.clear();
	for (int  i = 0; i < length; i++) {
		alignment_ += alignment[i];
	}
}

Alignment::Alignment() {

}

