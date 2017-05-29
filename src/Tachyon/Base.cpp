#include "Base.h"
#include<iostream>
#include<fstream>
#include<sstream>

#include<boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


#include "../Utils/util.h"

#include "../Utils/threadpool11/include/threadpool11/threadpool11.hpp"

using namespace std;
using namespace boost::serialization;

void Base::count(DatabaseElement &elem, std::unordered_map<int, int> &counters) {
	vector<long> results;

	Seg seg = Seg(seg_window_, seg_low_, seg_high_);
	std::string seq = elem.getSequence();
	seg.mask(seq);
	get_codes(seq, results, kmer_len_);

	for (auto &&code : results) {
		#pragma omp critical
		{
			if (counters.find(code) == counters.end()) {
				counters[code] = 1;
			} else {
				int cn = counters[code];
				counters[code] = cn + 1;
			}

		}
	}
}

void Base::findInRegion(std::string &region, std::vector<pair<long, long>> &indexes, std::unordered_map<int, int> &counters,
                        int kmer_numb, int region_start,
                        bool (*sortFunction)(const tuple<long, long, long> &, const tuple<long, long, long> &)) {

	vector<pair<long, long> > results;
	get_codes(region, results, kmer_len_);
	vector<tuple<long, long, long> > sorted_pentapeptides;

	for (auto &&p : results) {
		if (counters.find(p.first) == counters.end()) continue;
		int cn = counters[p.first];
		if (cn == 0) continue;
		sorted_pentapeptides.push_back(std::make_tuple(p.first, region_start+p.second, cn));
	}

	stable_sort(sorted_pentapeptides.begin(), sorted_pentapeptides.end(), sortFunction);

	vector<pair<long, long>> added;
	int added_size = 0;
	for (auto &&p : sorted_pentapeptides) {
		long code, position;
		std::tie(code, position, std::ignore) = p;
		if (added_size == high_numb_) break;
		bool isValid = true;
		for (auto &&elem : added)
			if (abs(position - elem.second) < kmer_len_) {
				isValid = false;
				break;
			}
		if (isValid) {
			added.push_back(make_pair(code, position));
			added_size++;
		}
	}

	if (added_size > 0) {
		int diff = kmer_numb - added_size;
		while (diff--) {
			added.push_back(added[0]);
		}
	}

	indexes.insert(indexes.end(), added.begin(), added.end());
}

void Base::find_indexes(DatabaseElement &elem, std::unordered_map<int, int> &counters,
                        std::unordered_map<long, vector<pair<long, long>>> &high_map,
                        std::unordered_map<long, vector<pair<long, long>>> &low_map) {


	vector<pair<long, long>> high_indexes;
	vector<pair<long, long>> low_indexes;

	int regionSize = high_olen_ == 0 ? elem.getSequenceLen() : high_olen_;
	for (int k = 0; k < elem.getSequenceLen(); k += regionSize) {
		string region = elem.getSequence().substr(k, regionSize);
		findInRegion(region, high_indexes, counters, high_numb_,k,
		             [](const tuple<long, long, long> &a, const tuple<long, long, long> &b) -> bool {
			             long freq_a, freq_b;
			             std::tie(std::ignore, std::ignore, freq_a) = a;
			             std::tie(std::ignore, std::ignore, freq_b) = b;

			             return freq_a > freq_b;
		             });
	}

	regionSize = low_olen_ == 0 ? elem.getSequenceLen() : low_olen_;
	for (int k = 0; k < elem.getSequenceLen(); k += regionSize) {
		string region = elem.getSequence().substr(k, regionSize);
		findInRegion(region, low_indexes, counters, low_numb_,k,
		             [](const tuple<long, long, long> &a, const tuple<long, long, long> &b) -> bool {
			             long freq_a, freq_b;
			             std::tie(std::ignore, std::ignore, freq_a) = a;
			             std::tie(std::ignore, std::ignore, freq_b) = b;

			             return freq_a < freq_b;
		             });
	}


	int id = elem.id();

	#pragma omp critical
	{
		for (auto &&p : high_indexes) {
			high_map[p.first].push_back(make_pair(id, p.second));
		}

		for (auto &&p : low_indexes) {
			low_map[p.first].push_back(make_pair(id, p.second));
		}
	}
}

void Base::make_indexes() {
	std::unordered_map<long, std::vector<std::pair<long, long>>> high_map;
	std::unordered_map<long, vector<pair<long, long>>> low_map;
	std::unordered_map<int, int> counters;

	cout << "Readed " << sets.size() << " sequences" << endl;
	cout << "Number of counted sequences:" << endl;
	#pragma omp parallel
	{
		#pragma omp single
		{
			for (int i = 0; i < sets.size(); i++) {
				#pragma omp task shared(counters)
				{
					count(sets[i], counters);
					if (i % 10000 == 0) {
						#pragma omp critical
						cout << i << endl;
					}
				}
			}
		}

		#pragma omp barrier
	}

	cout << "Finished with counting" << endl;
	cout << "Number of finished indexes:" << endl;
	#pragma omp parallel
	{
		#pragma omp single
		{
			for (int i = 0; i < sets.size(); i++) {
				#pragma omp task shared(high_map,low_map)
				{
					find_indexes(sets[i], counters, high_map, low_map);
					if (i % 10000 == 0) {
						#pragma omp critical
						cout << i << endl;
					}
				}
			}
		}
		#pragma omp barrier
	}
	high_freq_map = high_map;
	low_freq_map = low_map;
	writeMaps();
}

void Base::writeMaps() {
	string a = reduced_database_;
	std::ofstream out(a);
	boost::archive::binary_oarchive oa(out);

	oa << kmer_len_;
	oa << high_freq_map;
	oa << low_freq_map;
	out.close();

}

bool Base::read() {
	databaseSize_ = readFastaFile(pathToDatabase, sets, 0);
	if (!databaseSize_) {
		#pragma omp critical
		cerr << "Error: Original database doesn't exist" << endl;
		exit(-1);
	}
	return databaseSize_;
}


void Base::read_indexes() {
	std::ifstream out(reduced_database_);
	if (!out.good()) {
		#pragma omp critical
		cerr << "Error: Reduced database file doesn't exist" << endl;
		exit(-1);
	}
	boost::archive::binary_iarchive ia(out);
	ia >> kmer_len_;
	ia >> high_freq_map;
	ia >> low_freq_map;
	out.close();

	printf("This database was reduced with kmer length %d\n", kmer_len_);
}

bool Base::dump_in_memory() {

	threadpool11::Pool pool(2);

	pool.postWork<void>([&](){
		read();
	});
	pool.postWork<void>([&](){
		read_indexes();
	});

	pool.waitAll();

}
