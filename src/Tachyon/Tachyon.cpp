#include <string>
#include <unordered_set>

#include "../Utils/threadpool11/include/threadpool11/threadpool11.hpp"

#include "Tachyon.h"

Tachyon::Tachyon(const char* origPath, const char* reducedPath): base_(Ref<Base>(new Base(origPath, reducedPath))){
    base_->dump_in_memory();
}


void f(const Tachyon& alg, const DatabaseElement& query_, const Type type_, const RunParams& params_, EValue& evalueParams_,
                  ScoreMatrix& scorer_, vector<Alignment>& results_) {
        alg.findInDatabase(query_, type_, params_, evalueParams_, scorer_, results_);
        int resLen = results_.size();
        results_.resize(min(resLen, (int) params_.maxOutPerQuery));
}

void Tachyon::search(const std::vector<DatabaseElement>& queries, RunParams& params,
					 std::vector<std::vector<Alignment>>& results) const {

	params.scoreMatrixType = ScoreMatrixType::kBlosum62;
	params.alignType = AlignmentType::kSW;
    ScoreMatrix scorer(params.scoreMatrixType, params.gapOpen, params.gapExtend);
    EValue evalueParams(base_->database_size(), scorer);

	threadpool11::Pool pool(8);

	std::vector<std::future<void>> futures;

	results.clear();
	results.resize(queries.size());
	for (int i = 0; i < queries.size(); i++) {
		futures.emplace_back(pool.postWork<void>([&, i](){
			f(*this, queries[i], params.inputType, params, evalueParams, scorer, results[i]);
		}));
	}

	pool.joinAll();

//
//    #pragma omp parallel
//    {
//        #pragma omp single
//        {
//            for (int i = 0; i < queries.size(); i++) {
//                #pragma omp task shared(results)
//                {
//                    findInDatabase(queries[i], params.inputType, params, evalueParams, scorer, results[i]);
//                    int resLen = results[i].size();
//                    results[i].resize(min(resLen, (int) params.maxOutPerQuery));
//                }
//            }
//        }
//        #pragma omp barrier
//    }
}

int Tachyon::findInDatabase(const DatabaseElement& query, const Type type, const RunParams& params,
                            EValue& evalueParams, ScoreMatrix& scorer, std::vector<Alignment>& results) const {

//    if (type== Type::NUCLEOTIDES) {
//
//        std::vector<std::string> queriesAA;
//        convertDNAToAA(query.getSequence(), queriesAA);
//
//        int qLen = queriesAA.size();
//        std::vector<std::vector<Alignment>> inResults;
//        inResults.clear();
//        inResults.resize(qLen);
//
//        for (int i = 0; i < qLen; i++) {
//            DatabaseElement db(query.getName(), query.getNameLen(), queriesAA[i], queriesAA[i].length(), query.id());
//            #pragma omp task shared(inResults, evalueParams, scorer)
//            {
//                findInDatabase(db, Type::PROTEINS, params, evalueParams, scorer, inResults[i]);
//            }
//        }
//        #pragma omp taskwait
//        unordered_map<int, int> maxScores;
//        for (int i = 0; i < qLen; i++) {
//            for (const auto& it : inResults[i]) {
//                results.push_back(it);
//                maxScores[it.target_id()] = max(maxScores[it.target_id()], it.score());
//            }
//        }
//
//        //TODO: Group by getName and sort
//        sort(results.begin(), results.end(), compareAlignment);
//        return results.size();
//    }

    vector<pair<long, long>> codedPentapeptides;

    get_codes(query.getSequence(), codedPentapeptides, base_->kmer_len());
    std::unordered_set<long> hitsId;
    hit_database(codedPentapeptides, params.numOfHighMatch, params.numOfLowMatch, hitsId);
    ssearch(query, hitsId, results, params.maxEvalue, params.alignType, evalueParams, scorer);
}


void Tachyon::hit_database(vector<pair<long, long>>& coded_pentapeptides, int highMatch, int lowMatch, std::unordered_set<long> &results) const{

	std::unordered_map<long, vector<pair<long, long>>> high_positions;
	std::unordered_map<long, vector<pair<long, long>>> low_positions;

	unordered_map<long, int> potentialSequence;
	unordered_map<long, int> lowSequences;

	for (auto &code_pair : coded_pentapeptides) {
		if(base_->high_freq_map.find(code_pair.first) == base_->high_freq_map.end()) continue;
		for (auto &&p : base_->high_freq_map[code_pair.first]) {
			high_positions[p.first].push_back(make_pair(code_pair.second, p.second));
			potentialSequence[p.first] += 1;
		}

		if(base_->low_freq_map.find(code_pair.first) == base_->low_freq_map.end()) continue;
		for (auto &&p : base_->low_freq_map[code_pair.first]) {
			low_positions[p.first].push_back(make_pair(code_pair.second, p.second));
			lowSequences[p.first] += 1;
		}
	}

	filter_hits(potentialSequence, lowSequences, results, high_positions, low_positions, highMatch, lowMatch);
}


void Tachyon::filter_hits(unordered_map<long, int> &high_hits, unordered_map<long, int> &low_hits,
                          std::unordered_set<long> &results,
                          unordered_map<long, vector<pair<long, long>>> &positions,
                          unordered_map<long, vector<pair<long, long>>> &lowPositions,
                          int highMatch, int lowMatch) const{
	for (auto &&p : high_hits) {
		if (p.second >= highMatch) {
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
		if (p.second >= lowMatch) {
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


unsigned char *convert_str_seq_to_uchar_pt(const std::string &sequence, const int alphabetLength, const unsigned char *alphabet) {

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

void Tachyon::ssearch(const DatabaseElement& query, const std::unordered_set<long>& hits_id,
                      AlignmentSet& alignments, const double max_evalue, const AlignmentType align_type,
                      EValue& evalue_params, ScoreMatrix& scorer) const{

	unsigned char *query_ = convert_str_seq_to_uchar_pt(query.getSequence(), scorer.getAlphabetLength(), scorer.getAlphabet());
	int query_length = query.getSequenceLen();


	int db_len = hits_id.size();
	unsigned char *database_[db_len];
	int db_seq_lens[db_len];

	int k = 0;

	for (auto index : hits_id) {
		const auto &target = base_->getAt(index);
		database_[k] = convert_str_seq_to_uchar_pt(target.getSequence(),scorer.getAlphabetLength(),scorer.getAlphabet());
		db_seq_lens[k] = target.getSequenceLen();
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
		cout << "Error in opal! " << error << endl;
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

RunParams::RunParams(Type inputType, int numOfHighMatch, int numOfLowMatch, int gapOpen, int gapExtend,
                     double maxEvalue, double maxOutPerQuery, ScoreMatrixType scoreMatrixType, AlignmentType alignType)
        : inputType(inputType), numOfHighMatch(numOfHighMatch), numOfLowMatch(numOfLowMatch), gapOpen(gapOpen),
          gapExtend(gapExtend), maxEvalue(maxEvalue), maxOutPerQuery(maxOutPerQuery), scoreMatrixType(scoreMatrixType),
          alignType(alignType) {}
