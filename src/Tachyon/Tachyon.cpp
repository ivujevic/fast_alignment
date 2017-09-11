#include <string>
#include <unordered_set>

#include "Tachyon.h"

Tachyon::Tachyon(const char* origPath, const char* reducedPath) : base_(Base(origPath, reducedPath)) {
    base_.dumpInMemory();
}

void Tachyon::search(const std::vector<DatabaseElement>& queries, RunParams& params,
                     std::vector<std::vector<Alignment>>& results) const {

    params.scoreMatrixType = ScoreMatrixType::kBlosum62;
    params.alignType = AlignmentType::kSW;
    ScoreMatrix scorer(params.scoreMatrixType, params.gapOpen, params.gapExtend);
    EValue evalueParams(base_.numberOfElem(), scorer);

    #pragma omp parallel
    {
        #pragma omp single
        {
            for (int i = 0; i < queries.size(); i++) {
                #pragma omp task shared(results)
                {
                    findInDatabase(queries[i], params.inputType, params, evalueParams, scorer, results[i]);
                    int resLen = results[i].size();
                    results[i].resize(min(resLen, (int) params.maxOutPerQuery));
                }
            }
        }
        #pragma omp barrier
    }
}

int Tachyon::findInDatabase(const DatabaseElement& query, const Type type, const RunParams& params,
                            EValue& evalueParams, ScoreMatrix& scorer, std::vector<Alignment>& results) const {

    if (type == Type::NUCLEOTIDES) {

        std::vector<std::string> queriesAA;
        convertDNAToAA(query.getSequence(), queriesAA);

        int qLen = queriesAA.size();
        std::vector<std::vector<Alignment>> inResults;
        inResults.clear();
        inResults.resize(qLen);

        for (int i = 0; i < qLen; i++) {
            DatabaseElement db(query.getName(), query.getNameLen(), queriesAA[i], queriesAA[i].length(), query.id());
            #pragma omp task shared(inResults, evalueParams, scorer)
            {
                findInDatabase(db, Type::PROTEINS, params, evalueParams, scorer, inResults[i]);
            }
        }
        #pragma omp taskwait
        for (int i = 0; i < qLen; i++) {
            for (const auto& it : inResults[i]) {
                results.emplace_back(std::move(it));
            }
        }

        //TODO: Group by getName and sort
        sort(results.begin(), results.end(), compareAlignment);
        return results.size();
    }

    vector<pair<long, long>> codedPentapeptides;

    getCodes(query.getSequence(), codedPentapeptides, base_.kmerLen());

    std::vector<const DatabaseElement*> hits;
    hitDatabase(codedPentapeptides, params.numOfHighMatch, params.numOfLowMatch, hits);
    ssearch(query, hits, results, params.maxEvalue, params.alignType, evalueParams, scorer);
}


void Tachyon::hitDatabase(const vector<pair<long, long>>& coded_pentapeptides, const int highMatch, const int lowMatch,
                          std::vector<const DatabaseElement*>& results) const {

    std::unordered_map<long, vector<pair<long, long>>> high_positions;
    std::unordered_map<long, vector<pair<long, long>>> low_positions;

    unordered_map<long, int> potentialSequence;
    unordered_map<long, int> lowSequences;

    for (auto& code_pair : coded_pentapeptides) {
        if (auto it = base_.highFreqMap.find(code_pair.first);
        it != base_.highFreqMap.end()) {
            for (auto& p : it->second) {
                high_positions[p.first].emplace_back(code_pair.second, p.second);
                potentialSequence[p.first]++;
            }
        }
        if (auto it = base_.lowFreqMap.find(code_pair.first);
        it != base_.lowFreqMap.end()) {
            for (auto& p : it->second) {
                low_positions[p.first].emplace_back(code_pair.second, p.second);
                lowSequences[p.first]++;
            }
        }
    }

    std::unordered_set<long> resultsId;
    filterHits(potentialSequence, lowSequences, high_positions, low_positions, highMatch, lowMatch, resultsId);
    results.resize(resultsId.size());
    int i = -1;
    for (const long it : resultsId) results[++i] = &base_[it];
}


void Tachyon::filterHits(const unordered_map<long, int>& high_hits, const unordered_map<long, int>& low_hits,
                         unordered_map<long, vector<pair<long, long>>>& positions,
                         unordered_map<long, vector<pair<long, long>>>& lowPositions,
                         const int highMatch, const int lowMatch, std::unordered_set<long>& results) const {

    for (auto&& p : high_hits) {
        if (p.second >= highMatch) {
            auto& v = positions[p.first];
            sort(v.begin(), v.end(), [](pair<long, long> a, pair<long, long> b) -> bool { return a.first < b.first; });
            bool indexOk = true;
            for (int i = 0; i < p.second - 1; i++)
                if (v[i + 1].second < v[i].second) {
                    indexOk = false;
                    break;
                }
            if (indexOk) {
                vector<long> ds;
                for (auto& ps : v) ds.push_back(ps.first - ps.second);
                sort(ds.begin(), ds.end());
                int mx = 0;
                int elem = ds[0];
                int cnt = 0;
                for (auto& di : ds) {
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
    for (auto&& p : low_hits) {
        if (p.second >= lowMatch) {
            auto& v = lowPositions[p.first];
            sort(v.begin(), v.end(), [](pair<long, long> a, pair<long, long> b) -> bool { return a.first < b.first; });
            bool indexOk = true;
            for (int i = 0; i < p.second - 1; i++)
                if (v[i + 1].second < v[i].second) {
                    indexOk = false;
                    break;
                }
            if (indexOk) {
                vector<long> ds;
                for (auto& ps : v) ds.push_back(ps.first - ps.second);
                sort(ds.begin(), ds.end());
                int mx = 0;
                int elem = ds[0];
                int cnt = 0;
                for (auto& di : ds) {
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


unsigned char*
convert_str_seq_to_uchar_pt(const std::string& sequence, const int alphabetLength, const unsigned char* alphabet) {

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
    unsigned char* ret = new unsigned char[size];

    int i = 0;
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

void Tachyon::ssearch(const DatabaseElement& query, const std::vector<const DatabaseElement*>& hits,
                      AlignmentSet& alignments, const double max_evalue, const AlignmentType align_type,
                      EValue& evalue_params, ScoreMatrix& scorer) const {

    unsigned char* query_ = convert_str_seq_to_uchar_pt(query.getSequence(), scorer.getAlphabetLength(),
                                                        scorer.getAlphabet());
    int query_length = query.getSequenceLen();


    int db_len = hits.size();
    unsigned char* database_[db_len];
    int db_seq_lens[db_len];

    int k = 0;

    for (const auto& it : hits) {
        database_[k] = convert_str_seq_to_uchar_pt(it->getSequence(), scorer.getAlphabetLength(), scorer.getAlphabet());
        db_seq_lens[k] = it->getSequenceLen();
        k++;
    }

    OpalSearchResult* results[db_len];
    for (uint32_t i = 0; i < db_len; ++i) {
        results[i] = new OpalSearchResult();
        opalInitSearchResult(results[i]);
    }


    auto error = opalSearchDatabase(query_, query_length, database_, db_len, db_seq_lens, scorer.gap_open(),
                                    scorer.gap_extend(),
                                    scorer.getMatrix(), scorer.getAlphabetLength(), results, OPAL_SEARCH_ALIGNMENT,
                                    convertAlignTypeToOpal(align_type), OPAL_OVERFLOW_BUCKETS);

    k = 0;
    for (const auto& it : hits) {
        if (results[k]->scoreSet == 1) {
            auto eValue = evalue_params.calculate(results[k]->score, query_length, db_seq_lens[k]);
            if (eValue <= max_evalue) {
                alignments.emplace_back(results[k]->score, eValue, query.id(), it);
                alignments.back().update(results[k]->startLocationQuery, results[k]->endLocationQuery,
                                         results[k]->startLocationTarget,
                                         results[k]->endLocationTarget, results[k]->alignment,
                                         results[k]->alignmentLength);
            }
        }
        ++k;
    }

    if (error) {
        cout << "Error in opal! " << error << endl;
    }

    std::sort(alignments.begin(), alignments.end(), compareAlignment);
    for (const auto& it: database_) {
        delete[] it;
    }
    delete[] query_;
    for (auto& it : results) {
        delete it;
    }
}

RunParams::RunParams(Type inputType, int numOfHighMatch, int numOfLowMatch, int gapOpen, int gapExtend,
                     double maxEvalue, double maxOutPerQuery, ScoreMatrixType scoreMatrixType, AlignmentType alignType)
        : inputType(inputType), numOfHighMatch(numOfHighMatch), numOfLowMatch(numOfLowMatch), gapOpen(gapOpen),
          gapExtend(gapExtend), maxEvalue(maxEvalue), maxOutPerQuery(maxOutPerQuery), scoreMatrixType(scoreMatrixType),
          alignType(alignType) {}
