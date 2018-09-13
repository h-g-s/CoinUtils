#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "clique.h"
#include "vectormgm.h"
#include "memory.h"
#include "macros.h"
#include "str_utils.h"

#define INI_CAP 1024

#define HASH_SIZE 8192

#define ISET_HASH_SIZE 256

static const size_t hash_numbers[] = {37, 31, 29, 17, 13, 11, 7, 1};
static const size_t n_hash_numbers = 8;

size_t int_vector_hash_code(const int n, const int idx[]);

/***
 * returns 1 if clique was already inserted,
 * 0 otherwise
 ***/
int clq_set_clique_already_inserted(const CliqueSet *clqSet, const ISet *iset);

struct _CliqueSet {
    ISet **cliques;
    int *W;

    /* indicates the position of each clique */
    ISet **hash;

    int cliquesCap;
    int numberOfCliques;
    int weightSum;
};

CliqueSet *clq_set_clone(const CliqueSet *clqSet) {
    CliqueSet *clone = xmalloc(sizeof(CliqueSet));

    clone->cliquesCap = clqSet->cliquesCap;
    clone->numberOfCliques = clqSet->numberOfCliques;
    clone->weightSum = clqSet->weightSum;
    clone->cliques = xmalloc(sizeof(ISet *) * clone->cliquesCap);
    clone->W = xmalloc(sizeof(int) * clone->cliquesCap);

    for (int i = 0; i < clone->cliquesCap; i++) {
        clone->cliques[i] = NULL;
    }

    for (int i = 0; i < clone->numberOfCliques; i++) {
        if (clqSet->cliques[i]) {
            clone->cliques[i] = iset_create(ISET_HASH_SIZE);
            iset_cpy(clone->cliques[i], clqSet->cliques[i]);
        }

        clone->W[i] = clqSet->W[i];
    }

    clone->hash = xmalloc(sizeof(ISet *) * HASH_SIZE);
    for (int i = 0; i < HASH_SIZE; i++) {
        if (clqSet->hash[i]) {
            clone->hash[i] = iset_create(ISET_HASH_SIZE);
            iset_cpy(clone->hash[i], clqSet->hash[i]);
        } else {
            clone->hash[i] = NULL;
        }
    }

    return clone;
}

CliqueSet *clq_set_load(const char *fileName) {
#define INT_CHARS 10

    char *line = xmalloc(LINE_SIZE * INT_CHARS);
    int *elements = xmalloc(LINE_SIZE * sizeof(int));
    int size = 0;
    int weight = 0;
    char **strVector;
    CREATE_STRING_VECTOR(strVector, LINE_SIZE, INT_CHARS);

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
                clq_set_add(clqSet, size, elements, weight);
                size = 0;
                weight = 0;
            }

            sscanf(line + 1, "%d", &weight);
            char *start;
            if ((start = strstr(line + 1, "]")) == NULL) {
                fprintf(stderr, "Invalid file format.\n");
                exit(EXIT_FAILURE);
            }
            size = splitString(strVector, start + 1, ' ', LINE_SIZE, INT_CHARS, 1);
            {
                int i, nextEl = 0;
                for (i = 0; (i < size); ++i) {
                    if ((strlen(strVector[i]) > 0) && (isdigit(strVector[i][0])))
                        elements[nextEl++] = atoi(strVector[i]) - 1;
                }
                size = nextEl;
            }
        } else {
            if ((strlen(line) == 0) || (digitsInLine(line, LINE_SIZE) == 0))
                continue;
            if (size == 0) {
                fprintf(stderr, "Invalid file format.\n");
                exit(EXIT_FAILURE);
            }
            int newElements = splitString(strVector, line, ' ', LINE_SIZE, INT_CHARS, 1);
            {
                int i;
                for (i = 0; (i < newElements); ++i)
                    elements[size + i] = atoi(strVector[i]);
            }
            size += newElements;
        }
    }
    if (size > 0)
        clq_set_add(clqSet, size, elements, weight);

    FREE_STRING_VECTOR(strVector);
    fclose(f);
    free(line);
    free(elements);


    return clqSet;
#undef INT_CHARS
}


CliqueSet *clq_set_create() {
    int i;

    CliqueSet *clqSet = xmalloc(sizeof(CliqueSet));

    clqSet->hash = xmalloc(sizeof(ISet*) * HASH_SIZE);
    for (i = 0; (i < HASH_SIZE); ++i)
        clqSet->hash[i] = NULL;

    clqSet->numberOfCliques = 0;
    clqSet->cliquesCap = INI_CAP;
    clqSet->cliques = xmalloc(sizeof(ISet *) * clqSet->cliquesCap);
    for (i = 0; (i < clqSet->cliquesCap); ++i)
        clqSet->cliques[i] = NULL;

    clqSet->W = xmalloc(sizeof(int) * clqSet->cliquesCap);
    clqSet->weightSum = 0;

    return clqSet;
}

int clq_set_add(CliqueSet *clqSet, const int size, const int nodes[], const int w) {
    ISet *tmpClique = iset_create(ISET_HASH_SIZE);

    for (int i = 0; i < size; i++)
        iset_add(tmpClique, nodes[i]);

    if (clq_set_clique_already_inserted(clqSet, tmpClique)) {
        iset_free(&tmpClique);
        return 0;
    }

    if (clqSet->cliquesCap < clqSet->numberOfCliques + 1) {
        int i;
        const int oldCap = clqSet->cliquesCap;
        clqSet->cliquesCap *= 2;
        clqSet->W = xrealloc(clqSet->W, sizeof(int) * clqSet->cliquesCap);
        clqSet->cliques = xrealloc(clqSet->cliques, sizeof(ISet *) * clqSet->cliquesCap);
        for (i = oldCap; (i < clqSet->cliquesCap); ++i)
            clqSet->cliques[i] = NULL;
    }

    if (!clqSet->cliques[clqSet->numberOfCliques])
        clqSet->cliques[clqSet->numberOfCliques] = iset_create(ISET_HASH_SIZE);

    const int *snodes = iset_elements(tmpClique);
    for (int i = 0; i < iset_n_elements(tmpClique); i++)
        iset_add(clqSet->cliques[clqSet->numberOfCliques], snodes[i]);

    clqSet->W[clqSet->numberOfCliques] = w;
    clqSet->weightSum += w;

    /* inserting into hash table */
    size_t hash_code = int_vector_hash_code(size, snodes);
#ifdef DEBUG
    assert( hash_code<HASH_SIZE );
#endif

    if(!clqSet->hash[hash_code])
        clqSet->hash[hash_code] = iset_create(ISET_HASH_SIZE);

    iset_add(clqSet->hash[hash_code], clqSet->numberOfCliques);

    clqSet->numberOfCliques++;

    iset_free(&tmpClique);

    return 1;
}

int clq_set_weight(const CliqueSet *clqSet, const int clique) {
    return clqSet->W[clique];
}


int clq_set_clique_size(const CliqueSet *clqSet, const int clique) {
    if (!clqSet->cliques[clique])
        return 0;
    return (iset_n_elements(clqSet->cliques[clique]));
}

const int *clq_set_clique_elements(const CliqueSet *clqSet, const int clique) {
    if (clqSet->cliques[clique])
        return iset_elements(clqSet->cliques[clique]);
    return NULL;
}

void clq_set_free(CliqueSet **clqSet) {
    int i;
    for (i = 0; (i < (*clqSet)->cliquesCap); ++i)
        if ((*clqSet)->cliques[i])
            iset_free(&(*clqSet)->cliques[i]);

    for (i = 0; (i < HASH_SIZE); ++i)
        if((*clqSet)->hash[i])
            iset_free(&(*clqSet)->hash[i]);

    free((*clqSet)->cliques);
    free((*clqSet)->hash);
    free((*clqSet)->W);
    free((*clqSet));
    (*clqSet) = NULL;
}


int clq_validate(const CGraph *cgraph, const int size, const int nodes[], int *n1, int *n2) {
    int i, j;

    *n1 = -1;
    *n2 = -1;
    const int sizeM1 = size - 1;
    for (i = 0; (i < sizeM1); ++i)
        for (j = i + 1; (j < size); ++j) {
            if ((!cgraph_conflicting_nodes(cgraph, nodes[i], nodes[j])) || (nodes[i] == nodes[j])) {
                *n1 = nodes[i];
                *n2 = nodes[j];
                return 0;
            }
        }

    return 1;
}

int clq_comp_int(const void *v1, const void *v2) {
    return (*((const int *) v1)) - (*((const int *) v2));
}

int clq_conflicts_with_all(const CGraph *cgraph, const int node, const int size, const int nodes[]) {
    int i;
    for (i = 0; (i < size); ++i)
        if (!cgraph_conflicting_nodes(cgraph, node, nodes[i]))
            return 0;

    return 1;
}

int clq_set_number_of_cliques(const CliqueSet *clqSet) {
    if (!clqSet)
        return 0;
    return clqSet->numberOfCliques;
}

void clq_set_print(const CliqueSet *clqSet) {
    int i;
    for (i = 0; (i < clq_set_number_of_cliques(clqSet)); ++i) {
        printf("[%d] ", clq_set_weight(clqSet, i));
        const int *el = clq_set_clique_elements(clqSet, i);
        int j;
        for (j = 0; (j < clq_set_clique_size(clqSet, i)); ++j)
            printf("%d ", el[j] + 1);
        printf("\n");
    }
}

int clq_set_weight_sum(CliqueSet *clqSet) {
    return clqSet->weightSum;
}

size_t int_vector_hash_code(const int n, const int idx[]) {
    size_t code = 0;

    code += n * hash_numbers[0];

    code += idx[0] * hash_numbers[1];

    size_t i;
    for (i = 1; (i < n); ++i)
        code += hash_numbers[i % n_hash_numbers] * idx[i];

    code = code % HASH_SIZE;

#ifdef DEBUG
    assert( code >=0 ); assert( code < HASH_SIZE );
#endif

    return code;
}

int clq_set_clique_already_inserted(const CliqueSet *clqSet, const ISet *iset) {
    const size_t hash_code = int_vector_hash_code(iset_n_elements(iset), iset_elements(iset));
    const ISet *hashEntry = clqSet->hash[hash_code];

    if(!hashEntry)
        return 0;

    int i;
    const int *hashEntryElements = iset_elements(hashEntry);
    for(i = 0; i < iset_n_elements(hashEntry); i++) {
        const int cliqueIndex = hashEntryElements[i];
        const ISet *otherClique = clqSet->cliques[cliqueIndex];

        if(iset_equals(iset, otherClique))
        	return 1;
    }

    return 0;
}

int clq_set_clique_has_element(const CliqueSet *clqSet, const int clique, const int element) {
    if (!clqSet->cliques[clique])
        return 0;

    if (iset_has(clqSet->cliques[clique], element))
        return 1;

    return 0;
}

void clq_set_save(const CGraph *cgraph, const CliqueSet *clqSet, const char *fileName) {
    FILE *f = fopen(fileName, "w");

    int i;
    for (i = 0; (i < clqSet->numberOfCliques); ++i) {
        const int size = clq_set_clique_size(clqSet, i);
        const int *elements = clq_set_clique_elements(clqSet, i);
        const int w = clq_set_weight(clqSet, i);
        fprintf(f, "[%d]", w);
        int j;
        for (j = 0; (j < size); ++j) {
            fprintf(f, " %d", elements[j] + 1);
            if (cgraph_get_node_name(cgraph, elements[j]))
                fprintf(f, "(%s)", cgraph_get_node_name(cgraph, elements[j]));
            fprintf(f, " ");
        }
        fprintf(f, "\n");
    }

    fclose(f);
}

void clq_set_clear(CliqueSet *clqSet) {
    /* clearing hash contents */
    {
        int i;
        for (i = 0; (i < clq_set_number_of_cliques(clqSet)); ++i) {
            const int n = clq_set_clique_size(clqSet, i);
            const int *el = clq_set_clique_elements(clqSet, i);
            size_t entry = int_vector_hash_code(n, el);

            if(clqSet->hash[entry])
                iset_clear(clqSet->hash[entry]);
        }
    }
    /* clearing cliques */
    {
        int i;
        for (i = 0; (i < clqSet->cliquesCap); ++i)
            if (clqSet->cliques[i])
                iset_clear(clqSet->cliques[i]);
    }

    clqSet->numberOfCliques = 0;
    clqSet->weightSum = 0;
    memset(clqSet->W, 0, sizeof(int) * clqSet->cliquesCap);
}

void clq_set_cpy(CliqueSet *clqs_target, const CliqueSet *clqs_source) {
    clq_set_clear(clqs_target);
    clq_set_add_cliques(clqs_target, clqs_source);
}

int clq_set_add_cliques(CliqueSet *clqs_target, const CliqueSet *clqs_source) {
    int result = 0;
    {
        int i;
        for (i = 0; (i < clq_set_number_of_cliques(clqs_source)); ++i) {
            const int *el = clq_set_clique_elements(clqs_source, i);
            const int weight = clq_set_weight(clqs_source, i);

            result += clq_set_add(clqs_target, clq_set_clique_size(clqs_source, i), el, weight);
        }
    }
    return result;
}

void clq_set_add_using_original_indexes(CliqueSet *target, const CliqueSet *source, const int orig[]) {
    int i;
    for (i = 0; (i < clq_set_number_of_cliques(source)); ++i) {
        ISet *tmps = iset_create(ISET_HASH_SIZE);
        const ISet *clique = clq_set_get_clique(source, i);
        const int weight = clq_set_weight(source, i);
        const int size = iset_n_elements(clique);
        const int *el = iset_elements(clique);
        iset_add_using_original_indexes(tmps, el, size, orig);
        clq_set_add(target, iset_n_elements(tmps), iset_elements(tmps), weight);
        iset_free(&tmps);
    }
}

const ISet *clq_set_get_clique(const CliqueSet *clqSet, const int idx) {
    if (!clqSet->cliques[idx])
        return NULL;
    return clqSet->cliques[idx];
}

int clq_dominates(const ISet *a, const ISet *b) {
    const int sizeA = iset_n_elements(a);
    const int sizeB = iset_n_elements(b);

    if (sizeA <= sizeB)
        return 0;

    /* a is bigger, checking if it has all elements of b */
    const int *elB = iset_elements(b);
    int i;
    for (i = 0; (i < sizeB); ++i)
        if (!iset_has(a, elB[i]))
            return 0;

    return 1;
}


