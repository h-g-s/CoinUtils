#include <cstdio>
#include <cstring>
#include <ctime>
#include <cassert>
#include <limits>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <vector>
#include "clique_separation.h"
#include "bron_kerbosch.h"
#include "clique_extender.h"

/* default values */
#define CLQ_SEP_DEF_VERBOSE             0
#define CLQ_SEP_DEF_MIN_VIOL            0.02
#define CLQ_SEP_DEF_MIN_FRAC            0.001
#define CLQ_SEP_DEF_MAX_IT_BK           (std::numeric_limits<int>::max()/10)
#define CLQ_SEP_DEF_CLQE_EXTEND         CLQEM_EXACT

double fracPart(const double x) {
    double nextInteger = ceil(x);
    double previousInteger = floor(x);
    return std::min(nextInteger - x, x - previousInteger);
}

struct _CliqueSeparation {
    /* original conflict graph */
    const CGraph *cgraph;

    /* indicates the variables that will be considered in the pre-processed separation graph */
    int *nindexes;

    /* clique extender */
    CliqueExtender *clqe;

    /* clique set with original nodes, only translated cliques */
    CliqueSet *clqSetOrig;

    /* final clique set, possibly with additional extended cliques */
    CliqueSet *clqSet;

    /* minViol: default 0.02 */
    double minViol;

    /* extendCliques: 0: off - 1: random (default) - 2: max degree selection - 3: greedy extension - 4: exact extension */
    int extendCliques;

    /* costs based in reduced cost for extending cliques */
    double *costs;
    char hasCosts;

    /* minimum fractional value that a variable must have to be considered for separation */
    double minFrac;

    BronKerbosch *bk;

    /* max iterations for bron-kerbosch */
    int maxItBK;

    int verbose;
};

CliqueSeparation *clq_sep_create(const CGraph *origGraph) {
    CliqueSeparation *clqSep = new CliqueSeparation;
    const int cgSize = cgraph_size(origGraph);

    clqSep->minViol = CLQ_SEP_DEF_MIN_VIOL;
    clqSep->minFrac = CLQ_SEP_DEF_MIN_FRAC;
    clqSep->extendCliques = CLQ_SEP_DEF_CLQE_EXTEND;
    clqSep->verbose = CLQ_SEP_DEF_VERBOSE;
    clqSep->maxItBK = CLQ_SEP_DEF_MAX_IT_BK;
    clqSep->hasCosts = 0;
    clqSep->nindexes = new int[cgSize];
    clqSep->costs = new double[cgSize];
    clqSep->cgraph = origGraph;
    clqSep->clqe = clqe_create(origGraph);
    clqSep->clqSetOrig = clq_set_create();
    clqSep->clqSet = clq_set_create();

    return clqSep;
}

void clq_sep_set_rc(CliqueSeparation *sep, const double rc[]) {
    const int cgSize = cgraph_size(sep->cgraph);
    std::copy(rc, rc + cgSize, sep->costs);
    sep->hasCosts = 1;
}

void clq_sep_update_ppgraph_weights(CGraph *ppcg, const double x[]) {
    /* weights for fractional variables */
    for (int i = 0; i < cgraph_size(ppcg); i++) {
        const int origIdx = cgraph_get_original_node_index(ppcg, i);
        const double f = x[origIdx];
        cgraph_set_node_weight(ppcg, i, cgraph_weight(f));
    }
}

void clq_sep_separate(CliqueSeparation *sep, const double x[]) {
    const CGraph *cgraph = sep->cgraph;
    const int csize = cgraph_size(cgraph);
    const double minFrac = sep->minFrac;

    CliqueSet *clqSetOrig = sep->clqSetOrig;
    clq_set_clear(clqSetOrig); /* before extension, orig indexes */
    CliqueSet *clqSet = sep->clqSet;
    clq_set_clear(clqSet); /* final clique set */

    int *nindexes = sep->nindexes;
    int n = 0;
    for (int i = 0; i < csize; i++) {
        const int degree = cgraph_degree(cgraph, i);
        //disconsidering variables that are not binary
        //and variables that have conflict only with their complements
        if (degree < 2) {
            continue;
        } else if (fracPart(x[i]) < minFrac && x[i] < 0.98) {
            //variables at zero in x are disconsidered
            continue;
        } else {
            nindexes[n++] = i;
        }
    }

    CGraph *ppcg = cgraph_create_induced_subgraph(cgraph, nindexes, n);
    clq_sep_update_ppgraph_weights(ppcg, x);

    if (sep->verbose)
        cgraph_print_summary(ppcg, "pre-processed graph - part 1");

    /* separation works with integer weights */
    const int minW = (int) (1000.0 + (sep->minViol * 1000.0));

    if (cgraph_size(ppcg) >= 2) {
        sep->bk = bk_create(ppcg);
        clock_t startBK = clock();
        bk_set_max_it(sep->bk, sep->maxItBK);
        bk_set_min_weight(sep->bk, minW);
        bk_run(sep->bk);
        clock_t endBK = clock();
        if (sep->verbose) {
            printf("bk took %.3g seconds\n", ((double) (endBK - startBK)) / ((double) CLOCKS_PER_SEC));
        }

        const CliqueSet *bkClqSet = bk_get_clq_set(sep->bk);

        if (bkClqSet) {
            if (clq_set_number_of_cliques(bkClqSet)) {
#ifdef DEBUG
                for (int ic = 0; ic < clq_set_number_of_cliques(bkClqSet); ic++) {
                    const int size = clq_set_clique_size(bkClqSet, ic);
                    const int *el = clq_set_clique_elements(bkClqSet, ic);
                    int n1, n2;
                    if (!clq_validate(ppcg, el, size, &n1, &n2)) {
                        fprintf(stderr, "Nodes %d and %d are not in conflict in ppcg.\n", n1, n2);
                        exit(EXIT_FAILURE);
                    }
                    for (int j = 0; j < size; j++) {
                        const int vidx = el[j];
                        assert(vidx >= 0);
                        assert(vidx < cgraph_size(ppcg));
                    }
                }
#endif
                clq_set_add_using_original_indexes(clqSetOrig, bkClqSet, cgraph_get_original_node_indexes(ppcg));
            }
        }

        bk_free(&(sep->bk));
        sep->bk = nullptr;

        /* extending cliques */;
        if (sep->extendCliques != CLQEM_NO_EXTENSION) {
            clock_t startExtend = clock();
            CliqueExtender *clqe = sep->clqe;

            if (sep->hasCosts)
                clqe_set_costs(clqe, sep->costs, cgraph_size(cgraph));

            for (int i = 0; i < clq_set_number_of_cliques(clqSetOrig); i++) {
                const int *clqEl = clq_set_clique_elements(clqSetOrig, i);
                const int nClqEl = clq_set_clique_size(clqSetOrig, i);
                const int weight = clq_set_weight(clqSetOrig, i);
#ifdef DEBUG
                assert(clqEl && nClqEl);
                int en1, en2;
                if (!clq_validate(sep->cgraph, clqEl, nClqEl, &en1, &en2)) {
                    fprintf(stderr, "ERROR clqe_extend : Nodes %d and %d are not in conflict.\n", en1, en2);
                    exit(EXIT_FAILURE);
                }
#endif
                int nNewClqs = clqe_extend(clqe, clqEl, nClqEl, weight,
                                           static_cast<CliqueExtendingMethod>(sep->extendCliques));

                /* adding cliques which were not extended */
                if (!nNewClqs) {
                    clq_set_add(clqSet, clqEl, nClqEl, weight);
                }
            }

            /* adding all extended cliques */
            clq_set_add_cliques(clqSet, clqe_get_cliques(clqe));
            clock_t endExtend = clock();
            const double timeExtend = ((double) (endExtend - startExtend)) / ((double) CLOCKS_PER_SEC);
            if (sep->verbose) {
                printf("clique extension took %.3f seconds.\n", timeExtend);
            }
        } else {//no extension
            for (int i = 0; i < clq_set_number_of_cliques(clqSetOrig); i++) {
                const int *clqEl = clq_set_clique_elements(clqSetOrig, i);
                const int nClqEl = clq_set_clique_size(clqSetOrig, i);
                const int weight = clq_set_weight(clqSetOrig, i);
#ifdef DEBUG
                assert(clqEl && nClqEl);
                int en1, en2;
                if (!clq_validate(sep->cgraph, clqEl, nClqEl, &en1, &en2)) {
                    fprintf(stderr, "ERROR clqe_extend : Nodes %d and %d are not in conflict.\n", en1, en2);
                    exit(EXIT_FAILURE);
                }
#endif
                clq_set_add(clqSet, clqEl, nClqEl, weight);
            }
        }

#ifdef DEBUG
        for (int i = 0; i < clq_set_number_of_cliques(clqSet); i++) {
            assert(clq_set_clique_size(clqSet, i) && clq_set_clique_elements(clqSet, i));
            int en1, en2;
            if (!clq_validate(sep->cgraph, clq_set_clique_elements(clqSet, i), clq_set_clique_size(clqSet, i), &en1,
                              &en2)) {
                fprintf(stderr, "ERROR clqe_extend : Nodes %d and %d are not in conflict.\n", en1, en2);
                exit(EXIT_FAILURE);
            }
        }
#endif
    }

    /* need to be informed again next call */
    sep->hasCosts = 0;

    cgraph_free(&ppcg);
}

const CliqueSet *clq_sep_get_cliques(CliqueSeparation *sep) {
    return sep->clqSet;
}

void clq_sep_free(CliqueSeparation **clqSep) {
    clqe_free(&((*clqSep)->clqe));
    clq_set_free(&((*clqSep)->clqSetOrig));
    clq_set_free(&((*clqSep)->clqSet));

    delete[] (*clqSep)->nindexes;
    delete[] (*clqSep)->costs;

    delete (*clqSep);

    *clqSep = nullptr;
}

void clq_sep_set_max_it_bk(CliqueSeparation *clqSep, int maxItBK) {
    clqSep->maxItBK = maxItBK;
}

void clq_sep_set_max_it_bk_ext(CliqueSeparation *clqSep, int maxItBK) {
    clqe_set_max_it_bk(clqSep->clqe, maxItBK);
}

void clq_sep_set_extend_method(CliqueSeparation *sep, const int extendC) {
#ifdef DEBUG
    assert(extendC >= 0 && extendC <= 4);
#endif
    sep->extendCliques = extendC;
}

