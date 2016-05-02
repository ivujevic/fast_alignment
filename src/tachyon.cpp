#include <string>
#include <unordered_set>
#include "tachyon.h"


void Tachyon::search(ChainSet &queries, Type type, OutSet &results, AlignmentType align_type, double max_evalue,
                     int max_out, EValue &evalue_params, ScoreMatrix &scorer) {

	#pragma omp parallel
	{
		#pragma omp single
		{
			for (int i = 0; i < queries.size(); i++) {
				#pragma omp task shared(results)
				{
					findInDatabase(queries[i], type, results[i], align_type, max_evalue, evalue_params, scorer);
					int resLen = results[i].size();
					results[i].resize(min(resLen, max_out));
				}
			}
		}
		#pragma omp barrier
	}
}

int Tachyon::findInDatabase(DatabaseElement &query, Type type, AlignmentSet &results,
                            AlignmentType align_type, double max_evalue,
                            EValue &evalue_params, ScoreMatrix &scorer) {
	if (type == Type::NUCLEOTIDES) {

				vector<string> queries_AA;
				convertDNAToAA(query.name(),queries_AA);

				int q_len = queries_AA.size();
				OutSet  in_results(q_len);

				for(int i =0 ; i < q_len; i++) {
					DatabaseElement db(query.name(),query.name_len(),queries_AA[i],queries_AA[i].length(),query.id());
					#pragma omp task shared(in_results,evalue_params,scorer)
					{
						findInDatabase(db,Type::PROTEINS,in_results[i],align_type,max_evalue,evalue_params,scorer);
					}
				}

		unordered_map<int,int> max_scores;
		for(int i =0; i < q_len;i++) {
			for(auto& alignment : in_results[i]) {
				results.push_back(alignment);
				auto id = alignment.target_id();
				max_scores[id] = max(max_scores[id],alignment.score());
			}
		}

		//TODO: Group by name and sort
		sort(results.begin(),results.end(),compareAlignment);
	}

	vector<pair<long, long>> coded_pentapeptides;

	get_codes(query.sequence(), coded_pentapeptides);

	std::unordered_set<long> hits_id;

	hit_database(coded_pentapeptides, hits_id);


	ssearch(query, hits_id, results, max_evalue, align_type, evalue_params, scorer);
}

void Tachyon::hit_database(vector<pair<long, long>> &coded_pentapeptides, std::unordered_set<long> &results) {
	std::unordered_map<long, vector<pair<long, long>>> high_positions;
	std::unordered_map<long, vector<pair<long, long>>> low_positions;

	unordered_map<long, int> potentialSequence;
	unordered_map<long, int> lowSequences;

	for (auto &code_pair : coded_pentapeptides) {
		for (auto &&p : base_.high_freq_map[code_pair.first]) {

			high_positions[p.first].push_back(make_pair(code_pair.second, p.second));
			potentialSequence[p.first] += 1;
		}

		for (auto &&p : base_.low_freq_map[code_pair.first]) {
			low_positions[p.first].push_back(make_pair(code_pair.second, p.second));
			lowSequences[p.first] += 1;
		}
	}

	filter_hits(potentialSequence, lowSequences, results, high_positions, low_positions);
}


void Tachyon::filter_hits(unordered_map<long, int> &high_hits, unordered_map<long, int> &low_hits,
                          std::unordered_set<long> &results,
                          unordered_map<long, vector<pair<long, long>>> &positions,
                          unordered_map<long, vector<pair<long, long>>> &lowPositions) {
	for (auto &&p : high_hits) {
		if (p.second >= 3) {
			auto &v = positions[p.first];
			sort(v.begin(), v.end(), [](pair<long, long> a, pair<long, long> b) -> bool { return a.first < b.first; });
			bool indexOk = true;
			for (int i = 0; i < p.second - 1; i++)
				if (v[i + 1].second < v[i].second) {
					indexOk = false;
					break;
				}
			if (indexOk) {
				vector<long> ds;
				for (auto &ps : v) ds.push_back(ps.first - ps.second);
				sort(ds.begin(), ds.end());
				int mx = 0;
				int elem = ds[0];
				int cnt = 0;
				for (auto &di : ds) {
					if (elem == di) cnt++;
					else {
						elem = di;
						mx = max(mx, cnt);
						cnt = 1;
					}
				}
				mx = max(mx, cnt);
				if (mx > 1)
					results.insert(p.first);
			}
		}
	}
	for (auto &&p : low_hits) {
		if (p.second >= 2) {
			auto &v = lowPositions[p.first];
			sort(v.begin(), v.end(), [](pair<long, long> a, pair<long, long> b) -> bool { return a.first < b.first; });
			bool indexOk = true;
			for (int i = 0; i < p.second - 1; i++)
				if (v[i + 1].second < v[i].second) {
					indexOk = false;
					break;
				}
			if (indexOk) {
				vector<long> ds;
				for (auto &ps : v) ds.push_back(ps.first - ps.second);
				sort(ds.begin(), ds.end());
				int mx = 0;
				int elem = ds[0];
				int cnt = 0;
				for (auto &di : ds) {
					if (elem == di) cnt++;
					else {
						elem = di;
						mx = max(mx, cnt);
						cnt = 1;
					}
				}
				mx = max(mx, cnt);
				if (mx > 1) {
					results.insert(p.first);
				}
			}
		}
	}

}


unsigned char *convert_str_seq_to_uchar_pt(const std::string &sequence, int alphabetLength, unsigned char *alphabet) {

	unsigned char letterIdx[128];
	for (int i = 0; i < alphabetLength; i++) {
		if (alphabet[i] == '*') {
			for (int j = 0; j < 128; j++)
				letterIdx[j] = i;
			break;
		}
	}

	for (int i = 0; i < alphabetLength; i++)
		letterIdx[alphabet[i]] = i;


	int size = sequence.size();
	unsigned char *ret = new unsigned char[size];

	int i =0;
	for (char c : sequence)
		if (c != '*')
			ret[i++] = letterIdx[toupper(c)];
		else
			ret[i++] = letterIdx[c];

	return ret;
}

uint32_t convertAlignTypeToOpal(AlignmentType type) {
	switch (type) {
		case AlignmentType::kNW:
			return OPAL_MODE_NW;
		case AlignmentType::kHW:
			return OPAL_MODE_HW;
		case AlignmentType::kOV:
			return OPAL_MODE_OV;
		case AlignmentType::kSW:
			return OPAL_MODE_SW;
	}
}

void Tachyon::ssearch(DatabaseElement &query, std::unordered_set<long> hits_id,
                      AlignmentSet &alignments, double max_evalue, AlignmentType align_type,
                      EValue &evalue_params, ScoreMatrix &scorer) {

	unsigned char *query_ = convert_str_seq_to_uchar_pt(query.sequence(),scorer.getAlphabetLength(),scorer.getAlphabet());
	int query_length = query.sequence_len();


	int db_len = hits_id.size();
	unsigned char *database_[db_len];
	int db_seq_lens[db_len];

	int k = 0;

	for (auto index : hits_id) {
		const auto &target = base_[index];
		database_[k] = convert_str_seq_to_uchar_pt(target.sequence(),scorer.getAlphabetLength(),scorer.getAlphabet());
		db_seq_lens[k] = target.sequence_len();
		k++;
	}


	OpalSearchResult *results[db_len];
	for (uint32_t i = 0; i < db_len; ++i) {
		results[i] = new OpalSearchResult();
		opalInitSearchResult(results[i]);
	}


	auto error = opalSearchDatabase(query_, query_length, database_, db_len, db_seq_lens, scorer.gap_open(),
	                                scorer.gap_extend(),
	                                scorer.getMatrix(), scorer.getAlphabetLength(), results, OPAL_SEARCH_ALIGNMENT,
	                                convertAlignTypeToOpal(align_type), OPAL_OVERFLOW_BUCKETS);

	k = 0;
	for (auto &index : hits_id) {
		if (results[k]->scoreSet == 1) {
			auto evalue = evalue_params.calculate(results[k]->score, query_length, db_seq_lens[k]);
			if (evalue <= max_evalue) {
				alignments.push_back(Alignment(results[k]->score, evalue, query.id(), index));
				alignments.back().update(results[k]->startLocationQuery, results[k]->endLocationQuery,
				                         results[k]->startLocationTarget,
				                         results[k]->endLocationTarget, results[k]->alignment,
				                         results[k]->alignmentLength);
			}
		}
		k++;
	}

	if (error) {
		cout << "Greska " << error << endl;
	}
	std::sort(alignments.begin(), alignments.end(), compareAlignment);
	for (const auto &it: database_) {
		delete[] it;
	}
	delete[] query_;
	for (auto &it : results) {
		delete it;
	}
}