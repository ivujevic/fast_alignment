#pragma once
#include<string>
#include<memory>
#include <unordered_set>

#include "Base.h"
#include "DatabaseElement.h"
#include "../Utils/util.h"
#include "../Utils/Bioinformatics/Alignment/Opal.h"
#include "../Utils/Bioinformatics/ScoreMatrix.hpp"
#include "../Utils/Bioinformatics/Evalue.h"

#include "../Utils/Bioinformatics/Alignment.hpp"
#include "../Utils/Ref.hpp"

using OutSet = std::vector<AlignmentSet>;

using ChainSetP = std::vector<DatabaseElement*>;

class TachyonWorker;

class RunParams{
public:
    Type inputType;
	int numOfHighMatch;
	int numOfLowMatch;
	int gapOpen;
	int gapExtend;
	double maxEvalue;
    double maxOutPerQuery;
    ScoreMatrixType scoreMatrixType;
    AlignmentType alignType;

    RunParams(Type inputType, int numOfHighMatch, int numOfLowMatch, int gapOpen, int gapExtend, double maxEvalue,
              double maxOutPerQuery, ScoreMatrixType scoreMatrixType, AlignmentType alignType);
};

class Tachyon{
public:

    Tachyon(const char* origPath, const char* reducedPath);
	void search(const std::vector<DatabaseElement>& queries, RunParams& params, std::vector<std::vector<Alignment>>& results) const;

    const Ref<Base>& getBase() const { return base_; }
	int findInDatabase(const DatabaseElement& query, const Type type, const RunParams& params, EValue& evalueParams,
					   ScoreMatrix& scorer, std::vector<Alignment>& results) const;
private:

    Ref<Base> base_;


	void hit_database(vector<pair<long,long>>& coded_pentapeptides, int highMatch, int lowMatch, std::unordered_set<long>& results) const;

	void filter_hits(unordered_map<long,int>& high_hits, unordered_map<long,int>& low_hits, std::unordered_set<long>& results,
	                 unordered_map<long, vector<pair<long,long>>> &positions,
	                 unordered_map<long, vector<pair<long,long>>> &lowPositions,
                     int highMatch, int lowMatch) const;

	void ssearch(const DatabaseElement& query, const std::unordered_set<long>& hits_id, AlignmentSet & alignments,
	             const double max_evalue, const AlignmentType align_type, EValue& evalue_params, ScoreMatrix& scorer) const;

    friend class TachyonWorker;
};
