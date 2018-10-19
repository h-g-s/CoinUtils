#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <limits>
#include <algorithm>
#include "cgraph.h"
#include "clique.h"
#include "clique_extender.h"
#include "bron_kerbosch.h"

#define CLQE_DEF_MAX_CANDIDATES   256
#define CLQE_DEF_MAX_COST         100.0
#define CLQE_DEF_MAX_GEN          5
#define CLQE_DEF_RC_PERCENTAGE    0.6
#define CLQE_DEF_MAX_IT_BK        (std::numeric_limits<int>::max()/10)

struct _CliqueExtender {
    const CGraph *cgraph;
    CliqueSet *clqSet;
    int *candidates;

    int *newClique;
    int newCliqueSize;

    int maxCandidates;
    double maxCost;
    int maxClqGen;
    int maxItBK;
    double rcPercentage;

    double *costs;
};

CliqueExtender *clqe_create(const CGraph *cgraph) {
    CliqueExtender *clqe = new CliqueExtender;
    const int cgSize = cgraph_size(cgraph);


    clqe->cgraph = cgraph;
    clqe->clqSet = clq_set_create();
    clqe->candidates = new int[cgSize];
    clqe->newClique = new int[cgSize];
    clqe->costs = nullptr;
    clqe->maxCandidates = CLQE_DEF_MAX_CANDIDATES;
    clqe->maxCost = CLQE_DEF_MAX_COST;
    clqe->maxClqGen = CLQE_DEF_MAX_GEN;
    clqe->rcPercentage = CLQE_DEF_RC_PERCENTAGE;
    clqe->maxItBK = CLQE_DEF_MAX_IT_BK;

    return clqe;
}

int clqe_insert_best_candidates(CliqueExtender *clqe, const int *clique, const int size, const int weight,
                                const CliqueExtendingMethod clqem) {
    const CGraph *cgraph = clqe->cgraph;
    const int cgSize = cgraph_size(cgraph);
    const int nCols = cgSize / 2;
#ifdef DEBUG
    assert((cgSize % 2) == 0);
#endif
    char *iv = new char[cgSize]();
    int nodeSD = -1, minDegree = std::numeric_limits<int>::max();

    clqe->newCliqueSize = 0;

    /* picking node with the smallest degree */
    /* adding clique elements in newClique */
    for (int i = 0; i < size; i++) {
        const int degree = cgraph_degree(cgraph, clique[i]);
        if (degree < minDegree) {
            minDegree = degree;
            nodeSD = clique[i];
        }
        clqe->newClique[clqe->newCliqueSize++] = clique[i];
    }

    const double *costs = clqe->costs;

    if (clqem == CLQEM_PRIORITY_GREEDY) { //clique extender method uses greedy selection (reduced cost)
        int nNeighs = cgraph_get_best_n_neighbors(cgraph, nodeSD, costs, clqe->candidates, clqe->maxCandidates);
#ifdef DEBUG
        assert(nNeighs >= 0 && nNeighs <= clqe->maxCandidates);
#endif
        for (int i = 0; i < nNeighs; i++) {
            const int selected = clqe->candidates[i];

            if ((costs[selected] > clqe->maxCost) || (iv[selected]))
                continue;

            /* need to have conflict with all nodes in clique and all others inserted */
            bool insert = true;
            for (int j = 0; j < clqe->newCliqueSize; j++) {
                int complement = (clqe->newClique[j] < nCols) ? (clqe->newClique[j] + nCols) : (clqe->newClique[j] -
                                                                                                nCols);
                if ((!cgraph_conflicting_nodes(cgraph, clqe->newClique[j], selected)) ||
                    (selected == clqe->newClique[j]) || (selected == complement)) {
                    insert = false;
                    break;
                }
            }
            if (insert) {
                clqe->newClique[clqe->newCliqueSize++] = selected;
                int complement = (selected < nCols) ? (selected + nCols) : (selected - nCols);
                iv[selected] = 1;
                iv[complement] = 1;
            }
        }
    } else if (clqem == CLQEM_MAX_DEGREE) {
        double *degree = new double[cgSize];
        const int maxDegree = cgraph_max_degree(cgraph);

        for (int i = 0; i < cgSize; i++)
            degree[i] = (double) (maxDegree - cgraph_degree(cgraph, i));

        int nNeighs = cgraph_get_best_n_neighbors(cgraph, nodeSD, degree, clqe->candidates, clqe->maxCandidates);
#ifdef DEBUG
        assert(nNeighs >= 0 && nNeighs <= clqe->maxCandidates);
#endif
        for (int i = 0; i < nNeighs; i++) {
            const int selected = clqe->candidates[i];

            if ((costs[selected] > clqe->maxCost) || (iv[selected]))
                continue;

            /* need to have conflict with all nodes in clique and all others inserted */
            bool insert = true;
            for (int j = 0; j < clqe->newCliqueSize; j++) {
                int complement = (clqe->newClique[j] < nCols) ? (clqe->newClique[j] + nCols) : (clqe->newClique[j] -
                                                                                                nCols);
                if ((!cgraph_conflicting_nodes(cgraph, clqe->newClique[j], selected)) ||
                    (selected == clqe->newClique[j]) || (selected == complement)) {
                    insert = false;
                    break;
                }
            }
            if (insert) {
                clqe->newClique[clqe->newCliqueSize++] = selected;
                int complement = (selected < nCols) ? (selected + nCols) : (selected - nCols);
                iv[selected] = 1;
                iv[complement] = 1;
            }
        }
        delete[] degree;
    } else { //clique extender method uses random selection
        int nConflicts = cgraph_get_all_conflicting(cgraph, nodeSD, clqe->candidates, cgSize);

        if (nConflicts < clqe->maxCandidates) {
            for (int i = 0; i < nConflicts; i++) {
                const int selected = clqe->candidates[i];

                if ((costs[selected] > clqe->maxCost) || (iv[selected]))
                    continue;

                bool insert = true;
                for (int j = 0; j < clqe->newCliqueSize; j++) {
                    int complement = (clqe->newClique[j] < nCols) ? (clqe->newClique[j] + nCols) : (clqe->newClique[j] -
                                                                                                    nCols);
                    if ((!cgraph_conflicting_nodes(cgraph, clqe->newClique[j], selected)) ||
                        (selected == clqe->newClique[j]) || (selected == complement)) {
                        insert = false;
                        break;
                    }
                }
                if (insert) {
                    clqe->newClique[clqe->newCliqueSize++] = selected;
                    int complement = (selected < nCols) ? (selected + nCols) : (selected - nCols);
                    iv[selected] = 1;
                    iv[complement] = 1;
                }
            }
        } else {
            int r;
            char *isSelected = new char[nConflicts]();
            int selected;

            for (int i = 0; i < clqe->maxCandidates; i++) {
                do {
                    r = rand() % nConflicts;
                    selected = clqe->candidates[r];
                } while (isSelected[r]);

                isSelected[r] = 1;

                if ((costs[selected] > clqe->maxCost) || (iv[selected]))
                    continue;

                bool insert = true;
                for (int j = 0; j < clqe->newCliqueSize; j++) {
                    int complement = (clqe->newClique[j] < nCols) ? (clqe->newClique[j] + nCols) : (clqe->newClique[j] -
                                                                                                    nCols);
                    if ((!cgraph_conflicting_nodes(cgraph, clqe->newClique[j], selected)) ||
                        (selected == clqe->newClique[j]) || (selected == complement)) {
                        insert = false;
                        break;
                    }
                }
                if (insert) {
                    clqe->newClique[clqe->newCliqueSize++] = selected;
                    int complement = (selected < nCols) ? (selected + nCols) : (selected - nCols);
                    iv[selected] = 1;
                    iv[complement] = 1;
                }
            }
            delete[] isSelected;
        }
    }

#ifdef DEBUG
    assert(clqe->newCliqueSize <= cgSize);
#endif

    delete[] iv;

    if (clqe->newCliqueSize == size)
        return 0;

#ifdef DEBUG
    int en1, en2;
    if (!clq_validate(clqe->cgraph, clqe->newClique, clqe->newCliqueSize, &en1, &en2)) {
        fprintf(stderr, "ERROR clqe_extend : Nodes %d and %d are not in conflict.\n", en1, en2);
        exit(EXIT_FAILURE);
    }
#endif

    return clq_set_add(clqe->clqSet, clqe->newClique, clqe->newCliqueSize, weight);
}

typedef struct {
    int idx;
    int weight;
} CliqueWeight;

bool cmp_clq_weight(const CliqueWeight &e1, const CliqueWeight &e2) {
    if (e1.weight != e2.weight)
        return e1.weight < e2.weight;
    return e1.idx < e2.idx;
}

int exact_clique_extension(CliqueExtender *clqe, const int *clique, const int size, const int weight) {
    int i, j, nCandidates = 0;
    int cgSize = cgraph_size(clqe->cgraph), newCliques = 0;
    int *nindexes = new int[cgSize];
    int *candidates = new int[cgSize];
    int nodeSD = -1, degree = std::numeric_limits<int>::max() / 10;
    double minRC = std::numeric_limits<double>::max() / 10.0, maxRC = -(std::numeric_limits<double>::max() / 10.0);
    const int nCols = cgSize / 2;
#ifdef DEBUG
    assert((cgSize % 2) == 0);
#endif
    char *iv = new char[cgSize]();

    std::fill(nindexes, nindexes + cgSize, -1);

    for (i = 0; i < size; i++) {
        clqe->newClique[i] = clique[i];
        if (cgraph_degree(clqe->cgraph, clique[i]) < degree) {
            nodeSD = clique[i];
            degree = cgraph_degree(clqe->cgraph, clique[i]);
        }
    }

    int *neighs = new int[cgSize];
    int n = cgraph_get_all_conflicting(clqe->cgraph, nodeSD, neighs, cgSize);

    for (i = 0; i < n; i++) {
        int neigh = neighs[i];

        if (clqe->costs[neigh] > clqe->maxCost || iv[neigh] == 1)
            continue;

        for (j = 0; j < size; j++) {
            int complement = (clique[j] < nCols) ? (clique[j] + nCols) : (clique[j] - nCols);
            if ((neigh == clique[j]) || (neigh == complement) ||
                (!cgraph_conflicting_nodes(clqe->cgraph, clique[j], neigh)))
                break;
        }
        if (j >= size) {
            candidates[nCandidates++] = neigh;
            minRC = std::min(minRC, clqe->costs[neigh]);
            maxRC = std::max(maxRC, clqe->costs[neigh]);

            int complement = (neigh < nCols) ? (neigh + nCols) : (neigh - nCols);
            iv[neigh] = 1;
            iv[complement] = 1;
        }
    }

    double maxValue = minRC + floor(clqe->rcPercentage * (maxRC - minRC));
    j = nCandidates;
    nCandidates = 0;
    for (i = 0; i < j; i++)
        if (clqe->costs[candidates[i]] <= maxValue)
            nindexes[nCandidates++] = candidates[i];

    delete[] neighs;
    delete[] candidates;

    if (nCandidates == 0) {
        delete[] nindexes;
        delete[] iv;
        return 0;
    }

    CGraph *cg = cgraph_create_induced_subgraph(clqe->cgraph, nindexes, nCandidates);
    delete[] nindexes;

    if (fabs(minRC - maxRC) < 0.00001) {
        for (i = 0; i < nCandidates; i++) {
            cgraph_set_node_weight(cg, i, 1);
        }
    } else {
        for (i = 0; i < nCandidates; i++) {
            int origIdx = cgraph_get_original_node_index(cg, i);
            double normRC = 1.0 - ((clqe->costs[origIdx] - minRC) / (maxRC - minRC));
#ifdef DEBUG
            assert(normRC > -0.00001 && normRC < 1.00001);
#endif
            cgraph_set_node_weight(cg, i, cgraph_weight(normRC));
        }
    }

    BronKerbosch *bk = bk_create(cg);
    bk_set_max_it(bk, clqe->maxItBK);
    bk_set_min_weight(bk, 0);
    bk_run(bk);
    const CliqueSet *clqSet = bk_get_clq_set(bk);
    const int numClqs = clq_set_number_of_cliques(clqSet);
    if (clqSet && numClqs) {
        if (numClqs <= clqe->maxClqGen) {
            for (i = 0; i < numClqs; i++) {
                const int *extClqEl = clq_set_clique_elements(clqSet, i);
                const int extClqSize = clq_set_clique_size(clqSet, i);
                clqe->newCliqueSize = size;
#ifdef DEBUG
                assert(extClqSize > 0);
#endif
                for (j = 0; j < extClqSize; j++) {
                    int origIdx = cgraph_get_original_node_index(cg, extClqEl[j]);
                    clqe->newClique[size + j] = origIdx;
                    clqe->newCliqueSize++;
                }
#ifdef DEBUG
                int en1, en2;
                if (!clq_validate(clqe->cgraph, clqe->newClique, clqe->newCliqueSize, &en1, &en2)) {
                    fprintf(stderr, "ERROR clqe_extend : Nodes %d and %d are not in conflict.\n", en1, en2);
                    exit(EXIT_FAILURE);
                }
#endif
                newCliques += clq_set_add(clqe->clqSet, clqe->newClique, clqe->newCliqueSize, weight);
            }
        } else {
            CliqueWeight *clqw = new CliqueWeight[numClqs];
            for (i = 0; i < numClqs; i++) {
                clqw[i].idx = i;
                clqw[i].weight = clq_set_clique_size(clqSet, i) * clq_set_weight(clqSet, i);
            }
            std::sort(clqw, clqw + numClqs, cmp_clq_weight);
            for (i = 0; i < clqe->maxClqGen; i++) {
                const int *extClqEl = clq_set_clique_elements(clqSet, clqw[i].idx);
                const int extClqSize = clq_set_clique_size(clqSet, clqw[i].idx);
                clqe->newCliqueSize = size;
#ifdef DEBUG
                assert(extClqSize > 0);
#endif
                for (j = 0; j < extClqSize; j++) {
                    int origIdx = cgraph_get_original_node_index(cg, extClqEl[j]);
                    clqe->newClique[size + j] = origIdx;
                    clqe->newCliqueSize++;
                }
#ifdef DEBUG
                int en1, en2;
                if (!clq_validate(clqe->cgraph, clqe->newClique, clqe->newCliqueSize, &en1, &en2)) {
                    fprintf(stderr, "ERROR clqe_extend : Nodes %d and %d are not in conflict.\n", en1, en2);
                    exit(EXIT_FAILURE);
                }
#endif
                newCliques += clq_set_add(clqe->clqSet, clqe->newClique, clqe->newCliqueSize, weight);
            }
            delete[] clqw;
        }
    }

    bk_free(&bk);
    cgraph_free(&cg);
    delete[] iv;

    return newCliques;
}

int clqe_extend(CliqueExtender *clqe, const int *clique, const int size, const int weight,
                CliqueExtendingMethod clqem) {
#ifdef DEBUG
    assert(clqem != CLQEM_NO_EXTENSION);
#endif

    if ((clqem == CLQEM_EXACT || clqem == CLQEM_PRIORITY_GREEDY) && !clqe->costs) {
        fprintf(stderr, "Warning: using random selection for extension since no costs were informed.\n");
        clqem = CLQEM_RANDOM;
    }

    int result = 0;

    if (clqem == CLQEM_EXACT) {
        result += exact_clique_extension(clqe, clique, size, weight);
    } else {
        result += clqe_insert_best_candidates(clqe, clique, size, weight, clqem);
    }

    return result;
}

const CliqueSet *clqe_get_cliques(CliqueExtender *clqe) {
    return clqe->clqSet;
}

void clqe_set_costs(CliqueExtender *clqe, const double *costs, const int n) {
#ifdef DEBUG
    assert(!clqe->costs);
    assert(cgraph_size(clqe->cgraph) == n);
#endif
    clqe->costs = new double[n];
    std::copy(costs, costs + n, clqe->costs);
}

void clqe_free(CliqueExtender **clqe) {
    delete[] (*clqe)->newClique;
    delete[] (*clqe)->candidates;

    clq_set_free(&(((*clqe)->clqSet)));

    delete[] (*clqe)->costs;

    delete (*clqe);
    (*clqe) = nullptr;
}

void clqe_set_max_candidates(CliqueExtender *clqe, const int max_size) {
#ifdef DEBUG
    assert(max_size > 0);
#endif
    clqe->maxCandidates = max_size;
}

void clqe_set_max_cost(CliqueExtender *clqe, const double maxCost) {
    clqe->maxCost = maxCost;
}

void clqe_set_max_clq_gen(CliqueExtender *clqe, const int maxClqGen) {
#ifdef DEBUG
    assert(maxClqGen > 0);
#endif
    clqe->maxClqGen = maxClqGen;
}

int clqe_get_max_clq_gen(CliqueExtender *clqe) {
    return clqe->maxClqGen;
}

void clqe_set_rc_percentage(CliqueExtender *clqe, const double rcPercentage) {
    clqe->rcPercentage = rcPercentage;
}

double clqe_get_rc_percentage(CliqueExtender *clqe) {
    return clqe->rcPercentage;
}

int clqe_get_max_it_bk(CliqueExtender *clqe) {
    return clqe->maxItBK;
}

void clqe_set_max_it_bk(CliqueExtender *clqe, int maxItBK) {
    clqe->maxItBK = maxItBK;
}
