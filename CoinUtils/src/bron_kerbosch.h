#ifndef BRON_KERBOSCH_H_INCLUDED
#define BRON_KERBOSCH_H_INCLUDED

#include "cgraph.h"
#include "clique.h"

typedef struct _BronKerbosch BronKerbosch;

BronKerbosch* bk_create(const CGraph *cgraph);
void bk_free(BronKerbosch **_bk);

int bk_run(BronKerbosch *bk);

const CliqueSet* bk_get_clq_set(const BronKerbosch *bk);
int bk_get_max_weight(const BronKerbosch *bk);

void bk_set_min_weight(BronKerbosch *bk, int minWeight);
void bk_set_max_it(BronKerbosch *bk, int maxIt);


#endif
