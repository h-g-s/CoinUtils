#ifndef CLIQUE_EXTENDER_H
#define CLIQUE_EXTENDER_H

/**
 * starting from a clique found in a preprocessed graph,
 * extends it to include more nodes from the original graph
 **/

#include "cgraph.h"

typedef struct _CliqueExtender CliqueExtender;

typedef enum {
    CLQEM_NO_EXTENSION = 0,
    CLQEM_RANDOM = 1,
    CLQEM_MAX_DEGREE = 2,
    CLQEM_PRIORITY_GREEDY = 3,
    CLQEM_EXACT = 4
} CliqueExtendingMethod;

CliqueExtender *clqe_create(const CGraph *cgraph);

/**
 * extends clique "clique". extended
 * cliques are stores in internal cliqueSet
 **/

int clqe_extend(CliqueExtender *clqe, const int *clique, const int size, const int weight, CliqueExtendingMethod clqem);

const CliqueSet *clqe_get_cliques(CliqueExtender *clqe);

/* sets up costs for n variables */
void clqe_set_costs(CliqueExtender *clqe, const double *costs, const int n);

void clqe_free(CliqueExtender **clqe);


/* parameters */
void clqe_set_max_candidates(CliqueExtender *clqe, const int max_size);

void clqe_set_max_cost(CliqueExtender *clqe, const int maxCost);

void clqe_set_max_clq_gen(CliqueExtender *clqe, const int maxClqGen);

int clqe_get_max_clq_gen(CliqueExtender *clqe);

void clqe_set_rc_percentage(CliqueExtender *clqe, const double rcPercentage);

double clqe_get_rc_percentage(CliqueExtender *clqe);

int clqe_get_max_it_bk(CliqueExtender *clqe);

void clqe_set_max_it_bk(CliqueExtender *clqe, int maxItBK);

#endif
