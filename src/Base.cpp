#pragma once
#include "Base.h"
#include<iostream>
#include<fstream>
#include<sstream>

#include<boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


#include "util.h"


using namespace std;
using namespace boost::serialization;

void Base::count(DatabaseElement& elem, std::unordered_map<int, int>& counters) {
	int len = elem.sequence().length();
	vector<long> results;

	Seg seg = Seg(Type::PROTEINS);
	string seq = elem.sequence();
	seg.mask(seq);
	get_codes(seq,results);

	for(auto&& code : results) {
		#pragma omp critical
		{
			if(counters.find(code) == counters.end()) {
				counters[code] = 1;
			}else {
				int cn = counters[code];
				counters[code] = cn +1;
			}

		}
	}
}

void Base::find_indexes(DatabaseElement& elem, std::unordered_map<int, int>& counters,std::unordered_map<long, vector<pair<long,long>> >& high_map,
                        std::unordered_map<long, vector<pair<long,long>> >& low_map) {
	int len = elem.sequence().length();
	vector< pair<long,long> > results;
	get_codes(elem.sequence(),results);

	vector<tuple<long,long,long> > sorted_pentapeptides;

	for(auto&& p : results) {
		int cn = counters[p.first];
		if(cn == 0) continue;
		sorted_pentapeptides.push_back(std::make_tuple(p.first,p.second,counters[p.first]));
	}

	stable_sort(sorted_pentapeptides.begin(),sorted_pentapeptides.end(), [](const tuple<long,long,long>& a , const tuple<long,long,long> & b)->bool {
		long freq_a,freq_b;
		std::tie(std::ignore,std::ignore,freq_a) =a;
		std::tie(std::ignore,std::ignore,freq_b) =b;

		return freq_a > freq_b;
	});

	vector<pair<long,long>> high_freqs;
	vector<pair<long,long>> low_freqs;

	int ucitano = 0;
	for( auto&& p : sorted_pentapeptides){
		long code,position;
		std::tie(code,position,std::ignore) = p;
		if(ucitano == 5) break;
		bool isValid = true;
		for(auto&& added : high_freqs) if(abs(position - added.second) < 5) {isValid = false; break;}
		if( isValid) {
			high_freqs.push_back(make_pair(code,position));
			ucitano++;
		}
	}

	//add last five pentapeptides
	for( int i = sorted_pentapeptides.size() -1, j = 0; i >= 0 && j < 7; i--){
		auto p = sorted_pentapeptides[i];
		long code,position;
		std::tie(code,position,std::ignore) = p;
		bool isValid = true;
		for(auto&& added : low_freqs) if(abs(position - added.second) < 5) {isValid = false;break;}
		if( isValid) {
			low_freqs.push_back(make_pair(code,position));
			j++;
		}
	}

	int h_size = high_freqs.size();

	if(h_size > 0) {
		int diff = 5 - h_size;
		while(diff--) {
			high_freqs.push_back(high_freqs[0]);
		}
	}

	int l_size = low_freqs.size();
	if(l_size > 0) {
		int diff = 7 - l_size;
		while(diff--) {
			low_freqs.push_back(low_freqs[0]);
		}
	}

	int id = elem.id();

	#pragma omp critical
	{
		for(auto&& p : high_freqs) {
			high_map[p.first].push_back(make_pair(id,p.second));
		}

		for(auto&& p : low_freqs) {
			low_map[p.first].push_back(make_pair(id,p.second));
		}
	}
}
void Base::make_indexes() {
	std::unordered_map<long, std::vector<std::pair<long,long>> > high_map;
	std::unordered_map<long, vector<pair<long,long>> > low_map;
	std::unordered_map<int,int> counters;

	cout<<"Readed " <<sets.size()<<" sequences" << endl;
	#pragma omp parallel
	{
		#pragma omp single
		{
			for(int i =0; i < sets.size();i++){
				#pragma omp task shared(counters)
				{
					count(sets[i],counters);
				}
			}
		}

		#pragma omp barrier
	}

	cout<<"Finished with counting"<<endl;

	#pragma omp parallel
	{
		#pragma omp single
		{
			for(int i =0; i < sets.size();i++){
				#pragma omp task shared(high_map,low_map)
				{
					find_indexes(sets[i],counters,high_map,low_map);
					if(i%10000 == 0) {
						#pragma omp critical
						cout<<i<<endl;
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

	oa << high_freq_map;
	oa << low_freq_map;
	out.close();

}

bool Base::read() {
	databaseSize_ = readFastaFile(pathToDatabase,sets,0);
	if(!databaseSize_) {
		#pragma critical
		cerr<<"Error: Original database doesn't exist" << endl;
		exit(-1);
	}
	return databaseSize_;
}


void Base::read_indexes(std::unordered_map<long, vector<pair<long, long>>> &map) {
	std::ifstream out(reduced_database_);
	if(!out.good()) {
		#pragma critical
		cerr<<"Error: Reduce database file doesn't exist" << endl;
		exit(-1);
	}
	boost::archive::binary_iarchive ia(out);
	ia>>high_freq_map;
	ia>>low_freq_map;
	out.close();
}

bool Base::dump_in_memory() {
	#pragma omp parallel
	{
		#pragma omp single
		{
			#pragma omp task
			read();
			#pragma omp task
			read_indexes(high_freq_map);
		}
		#pragma omp barrier
	}
}