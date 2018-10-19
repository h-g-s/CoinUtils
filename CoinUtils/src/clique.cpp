#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cassert>
#include <cstring>
#include <algorithm>
#include <vector>
#include "clique.h"
#include "str_utils.h"

#define INI_CAP 1024
#define HASH_SIZE 8192
#define LINE_SIZE 2048

static const size_t hash_numbers[] = {37, 31, 29, 17, 13, 11, 7, 1};
static const size_t n_hash_numbers = 8;

size_t int_vector_hash_code(const std::vector<int> &clique);

/***
 * returns 1 if clique was already inserted,
 * 0 otherwise
 ***/
int clq_set_clique_already_inserted(const CliqueSet *clqSet, const std::vector<int> &clique);

struct _CliqueSet {
    std::vector<std::vector<int> > cliques;
    std::vector<int> W;

    /* indicates the position of each clique */
    std::vector<size_t> *hash;

    int weightSum;
};

CliqueSet *clq_set_clone(const CliqueSet *clqSet) {
    CliqueSet *clone = new CliqueSet;

    const int numberOfCliques = (int)clqSet->cliques.size();
    clone->weightSum = clqSet->weightSum;
    clone->cliques.resize(numberOfCliques);
    clone->W.resize(numberOfCliques);

    for (int i = 0; i < numberOfCliques; i++) {
        if (!clqSet->cliques[i].empty())
            clone->cliques[i] = clqSet->cliques[i];
        clone->W[i] = clqSet->W[i];
    }

    clone->hash = new std::vector<size_t>[HASH_SIZE];
    for (int i = 0; i < HASH_SIZE; i++) {
        if (!clone->hash[i].empty())
            clone->hash[i] = clqSet->hash[i];
    }

    return clone;
}

CliqueSet *clq_set_load(const char *fileName) {
#define INT_CHARS 10

    char *line = new char[LINE_SIZE * INT_CHARS];
    int *elements = new int[LINE_SIZE];
    int size = 0;
    int weight = 0;
    char **strVector = new char*[LINE_SIZE];

    strVector[0] = new char[LINE_SIZE * INT_CHARS];
    for(int i = 1; i < LINE_SIZE; i++)
        strVector[i] = strVector[i-1] + INT_CHARS;

    CliqueSet *clqSet = clq_set_create();

    FILE *f = fopen(fileName, "r");
    if (!f) {
        fprintf(stderr, "Could not open file %s", &(fileName[0]));
        exit(EXIT_FAILURE);
    }

    while (fgets(line, LINE_SIZE, f)) {
        if (line[strlen(line) - 1] == '\n')
            line[strlen(line) - 1] = '\0';

        if (line[0] == '[') {
            if (size > 0) {
                clq_set_add(clqSet, elements, size, weight);
                size = 0;
                weight = 0;
            }

            sscanf(line + 1, "%d", &weight);
            char *start;
            if ((start = strstr(line + 1, "]")) == nullptr) {
                fprintf(stderr, "Invalid file format.\n");
                exit(EXIT_FAILURE);
            }
            size = splitString(strVector, start + 1, ' ', LINE_SIZE, 1);
            int nextEl = 0;
            for (int i = 0; i < size; i++) {
                if ((strlen(strVector[i]) > 0) && (isdigit(strVector[i][0])))
#ifdef DEBUG
                    assert(nextEl + 1 < LINE_SIZE);
#endif
                    elements[nextEl++] = atoi(strVector[i]) - 1;
            }
            size = nextEl;
        } else {
            if ((strlen(line) == 0) || (digitsInLine(line, LINE_SIZE) == 0))
                continue;
            if (size == 0) {
                fprintf(stderr, "Invalid file format.\n");
                exit(EXIT_FAILURE);
            }
            int newElements = splitString(strVector, line, ' ', LINE_SIZE, 1);
            for (int i = 0; i < newElements; i++) {
#ifdef DEBUG
                assert(size + i < LINE_SIZE);
#endif
                elements[size + i] = atoi(strVector[i]);
            }
            size += newElements;
        }
    }
    if (size > 0)
        clq_set_add(clqSet, elements, size, weight);

    fclose(f);
    delete[] strVector[0];
    delete strVector;
    delete[] line;
    delete[] elements;

    return clqSet;
#undef INT_CHARS
}


CliqueSet *clq_set_create() {
    CliqueSet *clqSet = new CliqueSet;

    clqSet->hash = new std::vector<size_t>[HASH_SIZE];
    clqSet->cliques.reserve(INI_CAP);
    clqSet->W.reserve(INI_CAP);
    clqSet->weightSum = 0;

    return clqSet;
}

int clq_set_add( CliqueSet *clqSet, const int *idxs, const int size, const int w ) {
    std::vector<int> tmpClique(idxs, idxs + size);
    std::sort(tmpClique.begin(), tmpClique.end());
    auto it = std::unique(tmpClique.begin(), tmpClique.end());
    tmpClique.resize(static_cast<unsigned long>(std::distance(tmpClique.begin(), it)));

    if (clq_set_clique_already_inserted(clqSet, tmpClique))
        return 0;

    clqSet->cliques.push_back(tmpClique);
    clqSet->W.push_back(w);
    clqSet->weightSum += w;

    /* inserting into hash table */
    size_t hash_code = int_vector_hash_code(tmpClique);
#ifdef DEBUG
    assert(hash_code < HASH_SIZE);
#endif

    clqSet->hash[hash_code].push_back(clqSet->cliques.size());

    return 1;
}

int clq_set_add( CliqueSet *clqSet, const std::vector<int> &idxs, const int w ) {
    std::vector<int> tmpClique(idxs);
    std::sort(tmpClique.begin(), tmpClique.end());
    auto it = std::unique(tmpClique.begin(), tmpClique.end());
    tmpClique.resize(static_cast<unsigned long>(std::distance(tmpClique.begin(), it)));

    if (clq_set_clique_already_inserted(clqSet, tmpClique))
        return 0;

    clqSet->cliques.push_back(tmpClique);
    clqSet->W.push_back(w);
    clqSet->weightSum += w;

    /* inserting into hash table */
    size_t hash_code = int_vector_hash_code(tmpClique);
#ifdef DEBUG
    assert(hash_code < HASH_SIZE);
#endif

    clqSet->hash[hash_code].push_back(clqSet->cliques.size());

    return 1;
}

int clq_set_weight(const CliqueSet *clqSet, const int clique) {
    return clqSet->W[clique];
}

int clq_set_clique_size(const CliqueSet *clqSet, const int clique) {
    return (int)clqSet->cliques[clique].size();
}

const int* clq_set_clique_elements( const CliqueSet *clqSet, const int clique ) {
    return &clqSet->cliques[clique][0];
}

void clq_set_free(CliqueSet **clqSet) {
    delete[] (*clqSet)->hash;
    delete (*clqSet);
    (*clqSet) = nullptr;
}

int clq_validate(const CGraph *cgraph, const int *idxs, const int size, int *n1, int *n2) {
    *n1 = -1;
    *n2 = -1;
    for (int i = 0; i < size - 1; i++)
        for (int j = i + 1; j < size; j++) {
            if ((!cgraph_conflicting_nodes(cgraph, idxs[i], idxs[j])) || (idxs[i] == idxs[j])) {
                *n1 = idxs[i];
                *n2 = idxs[j];
                return 0;
            }
        }

    return 1;
}

int clq_set_number_of_cliques(const CliqueSet *clqSet) {
    if (!clqSet)
        return 0;
    return (int)clqSet->cliques.size();
}

void clq_set_print(const CliqueSet *clqSet) {
    for (size_t i = 0; i < clqSet->cliques.size(); i++) {
        printf("[%d] ", clqSet->W[i]);
        for (const int j : clqSet->cliques[i])
            printf("%d ", j + 1);
        printf("\n");
    }
}

int clq_set_weight_sum(const CliqueSet *clqSet) {
    return clqSet->weightSum;
}

size_t int_vector_hash_code(const std::vector<int> &clique) {
#ifdef DEBUG
    assert(!clique.empty());
#endif

    size_t code = 0;

    code += (clique.size() * hash_numbers[0]);
    code += (clique[0] * hash_numbers[1]);

    for (size_t i = 1; i < clique.size(); i++)
        code += (hash_numbers[i % n_hash_numbers] * clique[i]);

    code = (code % HASH_SIZE);

#ifdef DEBUG
    assert(code >= 0);
    assert(code < HASH_SIZE);
#endif

    return code;
}

int clq_set_clique_already_inserted(const CliqueSet *clqSet, const std::vector<int> &clique) {
    const size_t hash_code = int_vector_hash_code(clique);

    for(const size_t cliqueIndex : clqSet->hash[hash_code]) {
        if (clique == clqSet->cliques[cliqueIndex])
            return 1;
    }

    return 0;
}

void clq_set_save(const CGraph *cgraph, const CliqueSet *clqSet, const char *fileName) {
    FILE *f = fopen(fileName, "w");

    for (size_t i = 0; i < clqSet->cliques.size(); i++) {
        const int w = clqSet->W[i];
        fprintf(f, "[%d]", w);
        for (int element : clqSet->cliques[i]) {
            fprintf(f, " %d", element + 1);
            if (cgraph_get_node_name(cgraph, element))
                fprintf(f, "(%s)", cgraph_get_node_name(cgraph, element));
            fprintf(f, " ");
        }
        fprintf(f, "\n");
    }

    fclose(f);
}

void clq_set_clear(CliqueSet *clqSet) {
    /* clearing hash contents */
    for (size_t i = 0; i < HASH_SIZE; i++) {
        clqSet->hash[i].clear();
    }

    clqSet->cliques.clear(); clqSet->cliques.reserve(INI_CAP);
    clqSet->W.clear(); clqSet->W.reserve(INI_CAP);
    clqSet->weightSum = 0;
}

int clq_set_add_cliques(CliqueSet *clqs_target, const CliqueSet *clqs_source) {
    int result = 0;
    for (size_t i = 0; i < clqs_source->cliques.size(); i++) {
        result += clq_set_add(clqs_target, &clqs_source->cliques[i][0],
                              (int)clqs_source->cliques[i].size(), clqs_source->W[i]);
    }

    return result;
}

void clq_set_add_using_original_indexes(CliqueSet *target, const CliqueSet *source, const int *orig) {
    for(size_t i = 0; i < source->cliques.size(); i++) {
        const size_t cliqueSize = source->cliques[i].size();
        int *tmp = new int[cliqueSize];
        const int weight = source->W[i];

        for(size_t j = 0; j < cliqueSize; j++) {
            tmp[j] = orig[source->cliques[i][j]];
        }

        clq_set_add(target, tmp, (int)cliqueSize, weight);
        delete[] tmp;
    }
}

int clq_set_clique_has_element( const CliqueSet *clqSet, const int clique, const int element ) {
    return std::binary_search(clqSet->cliques[clique].begin(), clqSet->cliques[clique].end(), element);
}