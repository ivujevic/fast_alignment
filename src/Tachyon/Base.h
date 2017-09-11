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

#include "../Utils/Bioinformatics/Seg.h"
#include <tuple>

using namespace boost::serialization;


class Base {
public:

	Base(const char *path_to_database, const char *reduced_database, int kmer_len, int high_numb, int high_olen,
	     int low_numb,
	     int low_olen, int seg_window, double seg_low, double seg_high) : pathToDatabase(path_to_database),
	                                                                      reducedDatabase_(reduced_database),
	                                                                      kmerLen_(kmer_len),
	                                                                      highNumb_(high_numb), highOverLen_(high_olen),
	                                                                      lowNumb_(low_numb), lowOverLen_(low_olen),
	                                                                      segWindow_(seg_window), segHigh_(seg_high),
	                                                                      segLow_(seg_low) {
		printf("Creating base with this values:\n"
				       "\tkmer len = %d\n"
				       "\tnumber of high frequence kmers = %d\n"
				       "\tregion size for high frequence kmers = %d\n"
				       "\tnumber of low frerquence kmers = %d\n"
				       "\t region size for low frequence kmers = %d\n", kmerLen_, highNumb_, highOverLen_, lowNumb_,
		       low_olen);
	};


	Base(const char *path_to_database, const char *reduced_database) : pathToDatabase(path_to_database),
	                                                                   reducedDatabase_(reduced_database) { };

	const uint64_t numberOfElem() const { return numberOfElem_; }

	const DatabaseElement& operator[](int index) const {
		return sets[index];
	}

	/**
	 * Map where we store indexes for every sequence.
	 */
	void makeIndexes();

	bool dumpInMemory();

	bool read();

	void readBinary();

	std::unordered_map<long, vector<pair<long, long>>> highFreqMap;
	std::unordered_map<long, vector<pair<long, long>>> lowFreqMap;

	int kmerLen() const { return kmerLen_; }

private:
	int kmerLen_;
	int highNumb_;
	int highOverLen_;
	int lowNumb_;
	int lowOverLen_;

	int segWindow_;
	int segLow_;
	int segHigh_;

	vector<DatabaseElement> sets;

	const char *pathToDatabase;
	const char *reducedDatabase_;


	uint64_t numberOfElem_;


	void count(const DatabaseElement& elem, std::unordered_map<int, int>& counters);

	void findIndexes(const DatabaseElement& elem, const std::unordered_map<int, int>& counters,
                     std::unordered_map<long, vector<pair<long, long>>>& highMap,
                     std::unordered_map<long, vector<pair<long, long>>>& lowMap);

	void findInRegion(const std::string& sequence, std::vector<pair<long, long>>& indexes,
					  const std::unordered_map<int, int> &counters, const int kmerNumb, const int regionStart,
	                  bool (*f)(const tuple<long, long, long>& a, const tuple<long, long, long>& b));

	void writeMaps();

	void readIndexes();
};
