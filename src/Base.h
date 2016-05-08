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

	Base(const char *path_to_database, const char *reduced_database, int kmer_len, int high_numb, int high_olen,
	     int low_numb,
	     int low_olen, int seg_window, double seg_low, double seg_high) : pathToDatabase(path_to_database),
	                                                                      reduced_database_(reduced_database),
	                                                                      kmer_len_(kmer_len),
	                                                                      high_numb_(high_numb), high_olen_(high_olen),
	                                                                      low_numb_(low_numb), low_olen_(low_olen),
	                                                                      seg_window_(seg_window), seg_high_(seg_high),
	                                                                      seg_low_(seg_low) {
		printf("Creating base with this values:\n"
				       "\tkmer len = %d\n"
				       "\tnumber of high frequence kmers = %d\n"
				       "\tregion size for high frequence kmers = %d\n"
				       "\tnumber of low frerquence kmers = %d\n"
				       "\t region size for low frequence kmers = %d\n", kmer_len_, high_numb_, high_olen_, low_numb_,
		       low_olen);
	};


	Base(const char *path_to_database, const char *reduced_database) : pathToDatabase(path_to_database),
	                                                                   reduced_database_(reduced_database) { };


	const uint64_t database_size() { return databaseSize_; }

	const DatabaseElement &operator[](int index) {
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

	int kmer_len() { return kmer_len_; }

private:
	int kmer_len_;
	int high_numb_;
	int high_olen_;
	int low_numb_;
	int low_olen_;

	int seg_window_;
	int seg_low_;
	int seg_high_;

	vector<DatabaseElement> sets;

	Base(const Base &) = delete;

	const Base &operator=(const Base &) = delete;

	const char *pathToDatabase;
	const char *reduced_database_;


	uint64_t databaseSize_;


	void count(DatabaseElement &elem, std::unordered_map<int, int> &counters);

	void find_indexes(DatabaseElement &elem, std::unordered_map<int, int> &counters, std::unordered_map<long,
			vector<pair<long, long>>> &high_map, std::unordered_map<long, vector<pair<long, long>>> &low_map);

	void findInRegion(string &sequence, std::vector<pair<long, long>> &indexes, std::unordered_map<int, int> &counters,
	                  int kmer_numb,int region_start,
	                  bool (*f)(const tuple<long, long, long> &a, const tuple<long, long, long> &b));

	void writeMaps();

	void read_indexes();
};
