#pragma once
#include<string>
#include<memory>

#include "Base.h"
#include "util.h"
#include "DatabaseElement.h"
#include "opal.h"
#include "ScoreMatrix.hpp"
#include "Evalue.h"
#include <unordered_set>

#include "alignment.hpp"

using OutSet = std::vector<AlignmentSet>;

using ChainSetP = std::vector<DatabaseElement*>;
class Tachyon{
public:
	Tachyon(Base& base,int high_match,int low_match, int kmer_len):base_(base),high_match_(high_match),
	                                                                low_match_(low_match),kmer_len_(kmer_len){}

	void search(ChainSet& queries,Type type, OutSet & results, AlignmentType align_type,
	            double max_evalue, int max_uout,
	            EValue& evalue_params,ScoreMatrix& scorer);
private:

	int high_match_;
	int low_match_;
	int kmer_len_;

	Base& base_;

	int findInDatabase(DatabaseElement& query, Type type, AlignmentSet & results,
	                   AlignmentType align_type,double max_evalue,
	                   EValue& evalue_params,ScoreMatrix& scorer);

	void hit_database(vector<pair<long,long>>& coded_pentapeptides, std::unordered_set<long>& results);
	void filter_hits(unordered_map<long,int>& high_hits, unordered_map<long,int>& low_hits, std::unordered_set<long>& results,
	                 unordered_map<long, vector<pair<long,long>>> &positions,
	                 unordered_map<long, vector<pair<long,long>>> &lowPositions);

	void ssearch(DatabaseElement& query, std::unordered_set<long> hits_id, AlignmentSet & alignments,
	             double max_evalue,
	             AlignmentType align_type,EValue& evalue_params,ScoreMatrix& scorer);
};
