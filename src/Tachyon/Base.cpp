#include "Base.h"
#include<fstream>
#include<sstream>

#include "../Utils/util.h"

using namespace std;
using namespace boost::serialization;

void Base::count(const DatabaseElement& elem, std::unordered_map<int, int>& counters) {

    vector<long> results;
    Seg seg = Seg(segWindow_, segLow_, segHigh_);
    std::string seq = elem.getSequence();
    seg.mask(seq);
    getCodes(seq, results, kmerLen_);

    for (auto&& code : results) {
        #pragma omp critical
        {
            counters[code]++;
        }
    }
}

void Base::findInRegion(const std::string& region, std::vector<pair<long, long>>& indexes,
                        const std::unordered_map<int, int>& counters, const int kmerNumb, const int regionStart,
                        bool (* sortFunction)(const tuple<long, long, long>&, const tuple<long, long, long>&)) {

    vector<pair<long, long> > results;
    getCodes(region, results, kmerLen_);
    vector<tuple<long, long, long> > sorted_pentapeptides;

    for (auto&& p : results) {
        if (auto it = counters.find(p.first); it != counters.end() && it->second != 0) {
            sorted_pentapeptides.emplace_back(p.first, regionStart + p.second, it->second);
        }
    }

    stable_sort(sorted_pentapeptides.begin(), sorted_pentapeptides.end(), sortFunction);

    vector<pair<long, long>> added;
    int added_size = 0;
    for (auto&& p : sorted_pentapeptides) {
        long code, position;
        std::tie(code, position, std::ignore) = p;
        if (added_size == highNumb_) break;
        bool isValid = true;
        for (auto&& elem : added)
            if (abs(position - elem.second) < kmerLen_) {
                isValid = false;
                break;
            }
        if (isValid) {
            added.emplace_back(code, position);
            added_size++;
        }
    }

    if (added_size > 0) {
        int diff = kmerNumb - added_size;
        while (diff--) {
            added.push_back(added[0]);
        }
    }

    indexes.insert(indexes.end(), added.begin(), added.end());
}

void Base::findIndexes(const DatabaseElement& elem, const std::unordered_map<int, int>& counters,
                       std::unordered_map<long, vector<pair<long, long>>>& highMap,
                       std::unordered_map<long, vector<pair<long, long>>>& lowMap) {


    vector<pair<long, long>> highI;
    vector<pair<long, long>> lowIndexes;

    int regionSize = highOverLen_ == 0 ? elem.getSequenceLen() : highOverLen_;
    for (int k = 0; k < elem.getSequenceLen(); k += regionSize) {
        const string& region = elem.getSequence().substr(k, regionSize);
        findInRegion(region, highI, counters, highNumb_, k,
                     [](const tuple<long, long, long>& a, const tuple<long, long, long>& b) -> bool {
                         long freq_a, freq_b;
                         std::tie(std::ignore, std::ignore, freq_a) = a;
                         std::tie(std::ignore, std::ignore, freq_b) = b;

                         return freq_a > freq_b;
                     });
    }

    regionSize = lowOverLen_ == 0 ? elem.getSequenceLen() : lowOverLen_;
    for (int k = 0; k < elem.getSequenceLen(); k += regionSize) {
        const string& region = elem.getSequence().substr(k, regionSize);
        findInRegion(region, lowIndexes, counters, lowNumb_, k,
                     [](const tuple<long, long, long>& a, const tuple<long, long, long>& b) -> bool {
                         long freq_a, freq_b;
                         std::tie(std::ignore, std::ignore, freq_a) = a;
                         std::tie(std::ignore, std::ignore, freq_b) = b;

                         return freq_a < freq_b;
                     });
    }


    int id = elem.id();

    #pragma omp critical
    {
        for (auto&& p : highI) {
            highMap[p.first].emplace_back(id, p.second);
        }

        for (auto&& p : lowIndexes) {
            lowMap[p.first].emplace_back(id, p.second);
        }
    }
}

void Base::makeIndexes() {
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
                #pragma omp task shared(high_map, low_map)
                {
                    findIndexes(sets[i], counters, high_map, low_map);
                    if (i % 10000 == 0) {
                        #pragma omp critical
                        cout << i << endl;
                    }
                }
            }
        }
        #pragma omp barrier
    }

    highFreqMap = std::move(high_map);
    lowFreqMap = std::move(low_map);
    writeMaps();
}

void Base::writeMaps() {
    string a = reducedDatabase_;
    std::ofstream out(a);
    boost::archive::binary_oarchive oa(out);

    oa << kmerLen_;
    oa << highFreqMap;
    oa << lowFreqMap;
    out.close();

}

bool Base::read() {
    numberOfElem_ = readFastaFile(pathToDatabase, sets, 0);
    if (!numberOfElem_) {
        #pragma omp critical
        cerr << "Error: Original database doesn't exist" << endl;
        exit(-1);
    }
    return numberOfElem_;
}


void Base::readIndexes() {
    std::ifstream out(reducedDatabase_);
    if (!out.good()) {
        #pragma omp critical
        cerr << "Error: Reduced database file doesn't exist" << endl;
        exit(-1);
    }
    boost::archive::binary_iarchive ia(out);
    ia >> kmerLen_;
    ia >> highFreqMap;
    ia >> lowFreqMap;
    out.close();

    printf("This database was reduced with kmer length %d\n", kmerLen_);
}

bool Base::dumpInMemory() {
    #pragma omp parallel
    {
        #pragma omp single
        {
            #pragma omp task
            read();
            #pragma omp task
            readIndexes();
        }
        #pragma omp barrier
    }
}
