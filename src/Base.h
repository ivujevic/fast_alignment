/* 
 * File:   Base.h
 * Author: vujevic
 *
 * Created on April 25, 2014, 11:41 PM
 */

#pragma once
#include<string>
#include "DatabaseElement.h"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include "Seg.h"
#include <tuple>

using namespace boost::serialization;


class Base {
public:

	Base(const char *path_to_database, const char* path_to_hf) : pathToDatabase(path_to_database),
	                                                                                     reduced_database_(path_to_hf) { };


	const uint64_t database_size() {return databaseSize_;}

	const DatabaseElement& operator [] (int index) {
		return sets[index];
	}

	/**
	 * Map where we store indexes for every sequence.
	 */
	void make_indexes();
	bool dump_in_memory();
	bool read();
	void readBinary();
	std::unordered_map<long, vector<pair<long, long>>> high_freq_map;
	std::unordered_map<long, vector<pair<long, long>>> low_freq_map;
private:
	vector<DatabaseElement> sets;
	Base(const Base&) = delete;
	const Base& operator=(const Base&) = delete;

	const char *pathToDatabase;
	const char *reduced_database_;



	uint64_t databaseSize_;


	void count(DatabaseElement &elem, std::unordered_map<int, int> &counters);

	void find_indexes(DatabaseElement &elem, std::unordered_map<int, int> &counters,std::unordered_map<long,
			vector<pair<long,long>> >&high_map,std::unordered_map<long, vector<pair<long,long>> >&low_map);

	void writeMaps();

	void read_indexes();
};
