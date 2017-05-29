#include "../Tachyon/Tachyon.h"
#include <iostream>

#ifdef _WIN32
#   define DLL_EXPORT __declspec(dllexport)
#else
#   define DLL_EXPORT
#endif

class DatabaseElementWrapper {
public:
    int id_;
    char* name_;
    int nameLen_;
    char* sequence_;
    int sequenceLen_;
};
extern "C" {
DLL_EXPORT Tachyon* Tachyon_init(const char* origPath, const char* reducedPath) {
    return new Tachyon(origPath, reducedPath);
}

DLL_EXPORT void Tachyon_search(const Tachyon* self, const DatabaseElementWrapper* queries, int numbOfQueries, RunParams& params, Alignment** results, int** resultsPerQuery) {

    std::vector<DatabaseElement> queriesV;

    for (int i = 0; i < numbOfQueries; i++) {
        queriesV.emplace_back(DatabaseElement(queries[i].id_, queries[i].name_, queries[i].nameLen_, queries[i].sequence_, queries[i].sequenceLen_));
    }

    std::vector<std::vector<Alignment>>* v = new std::vector<std::vector<Alignment>>;
    self->search(queriesV, params, *v);
    int cn = 0;

    *resultsPerQuery = (int *) malloc(sizeof(int) * numbOfQueries);
    memset(*resultsPerQuery, 0, sizeof(int) * numbOfQueries);

    for (int i = 0; i < numbOfQueries; i++) {
        int size = (*v)[i].size();
        (*resultsPerQuery)[i] = -10;
        cn += size;
    }

    std::cout << *resultsPerQuery[0] << std::endl;
    *results = (Alignment *) malloc(sizeof(Alignment) * cn);
    int p = 0;
    for (int i = 0; i < numbOfQueries; i++) {
        for (int j = 0; j < (*v)[i].size(); j++) {
            results[p++] = &(*v)[i][j];
        }
    }

}
DLL_EXPORT void Tachyon_ispis(const Tachyon* self) {
    std::cout<<self->getBase()->database_size()<<std::endl;
}
}
