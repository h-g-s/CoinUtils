#include <vector>
#include <list>
#include <cassert>
#include <algorithm>
#include <climits>
#include <cstdio>

extern "C" {
    #include "bron_kerbosch.h"
    #include "memory.h"
}

#define INT_SIZE (8*sizeof(int))

using namespace std;

typedef struct {
    size_t id;
    size_t weight;
    size_t degree;
    size_t mdegree; //modified degree: degree of vertex + sum of degrees of adjacent vertices
} BKVertex;

typedef struct {
    list<size_t> vertices;
    size_t totalWeight;
} ListOfVertices;

typedef struct {
    size_t *vertices;
    size_t numVertices; //size of array "vertices"
    size_t totalWeight; //sum of weights of active vertices
    size_t filledPositions; //controls the number of active vertices
} ArrayOfVertices;

struct _BronKerbosch {
    //graph data
    const CGraph *cgraph;
    vector<BKVertex> vertices;

    //bk data
    int **bit;
    size_t *mask;

    //bk parameters
    size_t minWeight;
    size_t maxIt;

    //bk statistics
    size_t it;
    size_t maxWeight;
    int status;

    //bk result
    CliqueSet *clqSet;
};

ArrayOfVertices* array_of_vertices_create(size_t size) {
    ArrayOfVertices* av = new ArrayOfVertices;

    av->vertices = new size_t[size]();
    av->numVertices = size;
    av->totalWeight = 0L;
    av->filledPositions = 0L;

    return av;
}

void array_of_vertices_free(ArrayOfVertices *av) {
    delete[] av->vertices;
    delete av;
    av = nullptr;
}

ListOfVertices* list_of_vertices_create() {
    ListOfVertices* lv = new ListOfVertices;
    lv->totalWeight = 0L;
    return lv;
}

void list_of_vertices_free(ListOfVertices* lv) {
    delete lv;
    lv = nullptr;
}

BronKerbosch* bk_create(const CGraph *cgraph) {
    BronKerbosch *bk = new BronKerbosch;
    BKVertex aux;

    const int* weights = cgraph_get_node_weights(cgraph);
    const size_t nVertices = (const size_t) cgraph_size(cgraph);
    int *neighs = new int[nVertices*2];

    bk->vertices.reserve(nVertices);

    for(size_t i = 0; i < nVertices; i++) {
        int realDegree = cgraph_degree(cgraph, (int)i);
        int check = cgraph_get_all_conflicting( cgraph, (int)i, neighs, (int)nVertices*2 );
        assert(check == realDegree);

        size_t mdegree = (size_t)realDegree;
        for(int j = 0; j < realDegree; j++) {
            assert(neighs[j] != i && neighs[j] >= 0 && neighs[j] < nVertices);
            mdegree += cgraph_degree(cgraph, neighs[j]);
        }

        aux.id = i;
        aux.weight = (size_t)weights[i];
        aux.degree = (size_t)realDegree;
        aux.mdegree = mdegree;
        bk->vertices.push_back(aux);
    }

    assert(nVertices == bk->vertices.size());

    bk->cgraph = cgraph;
    bk->clqSet = clq_set_create();
    bk->status = 0L;
    bk->maxWeight = 0L;
    bk->it = 0L;
    bk->maxIt = (INT_MAX/100);
    bk->minWeight = 0L;

    bk->mask = new size_t[INT_SIZE];
    bk->mask[0] = 1;
    for(size_t h = 1; h < INT_SIZE; h++)
        bk->mask[h] = bk->mask[h-1] << 1U;

    bk->bit = new int*[nVertices]();
    for(int i = 0; i < nVertices; ++i)
        bk->bit[i] = new int[nVertices/INT_SIZE + 1]();

    for(size_t v = 0; v < nVertices; ++v) {
        for(size_t y = (v+1); y < nVertices; ++y) {
            if(cgraph_conflicting_nodes(bk->cgraph, (int)v, (int)y) == 1) {
                bk->bit[y][v/INT_SIZE] |= bk->mask[v%INT_SIZE];
                bk->bit[v][y/INT_SIZE] |= bk->mask[y%INT_SIZE];
            }
        }
    }

    delete[] neighs;

    return bk;
}

void bk_free(BronKerbosch **_bk) {
    BronKerbosch *bk = *_bk;

    clq_set_free(&bk->clqSet);

    for(size_t i = 0; i < bk->vertices.size(); i++)
        delete[] bk->bit[i];
    delete[] bk->bit;

    delete[] bk->mask;

    delete bk;
    *_bk = nullptr;
}

void exclude_neighbors_u(const BronKerbosch *bk, vector<size_t> &P_excluding_N_u, const ListOfVertices *P, size_t u) {
    for(const size_t vertex : P->vertices) {
        if(!cgraph_conflicting_nodes(bk->cgraph, (int)u, (int)vertex)) {
            P_excluding_N_u.push_back(vertex);
        }
    }
}

vector<int> generate_clique(const BronKerbosch *bk, const ArrayOfVertices *C) {
    size_t node, w = 0L, value;
    vector<int> nodes;

    nodes.reserve(C->numVertices);

    for(size_t t = 0; t < C->numVertices; t++) {
        node = INT_SIZE * t;
        value =  C->vertices[t];

        while(value > 1) {
            if(value % 2 == 1) {
                nodes.push_back((int)node);
                w += bk->vertices[node].weight;
            }

            value = value / 2;
            node++;
        }

        if(value == 1) {
            nodes.push_back((int)node);
            w += bk->vertices[node].weight;
        }
    }

    assert(w == C->totalWeight);

    /*int n1[nodes.size()], n2[nodes.size()];
    if(!clq_validate(bk->cgraph, (int)nodes.size(), &nodes[0], n1, n2)) {
        fprintf(stderr, "Clique invalido!\n");
        for(const int n : nodes)
            fprintf(stderr, "%d ", n);
        fprintf(stderr, "\n");
        exit(EXIT_FAILURE);
    }*/

    return nodes;
}

ArrayOfVertices* create_new_C(const BronKerbosch *bk, const ArrayOfVertices *C, const size_t v) {
    ArrayOfVertices* newC = array_of_vertices_create((int)bk->vertices.size()/INT_SIZE + 1);

    //copying elements of C
    copy(C->vertices, C->vertices + C->numVertices, newC->vertices);
    newC->numVertices = C->numVertices;
    newC->totalWeight = C->totalWeight;
    newC->filledPositions = C->filledPositions;

    //adding v in C
    newC->vertices[v/INT_SIZE] |= bk->mask[v%INT_SIZE];
    newC->filledPositions++;
    newC->totalWeight += bk->vertices[v].weight;

    return newC;
}

ArrayOfVertices* create_new_S(const BronKerbosch *bk, const ArrayOfVertices *S, const size_t v) {
    ArrayOfVertices *newS = array_of_vertices_create((int)bk->vertices.size()/INT_SIZE + 1);

    //newS = S intersection N(v)
    newS->filledPositions = 0L;
    for(size_t i = 0; i < S->numVertices; ++i) {
        newS->vertices[i] = S->vertices[i] & bk->bit[v][i];
        if(newS->vertices[i])
            newS->filledPositions++;
    }

    return newS;
}

ListOfVertices* create_new_P(const BronKerbosch *bk, const ListOfVertices *P, const size_t v) {
    ListOfVertices *newP = list_of_vertices_create();

    //newP = P intersection N(v)
    for(const size_t vertex : P->vertices) {
        if(cgraph_conflicting_nodes(bk->cgraph, (int)v, (int)vertex)) {
            newP->vertices.push_back(vertex);
            newP->totalWeight += bk->vertices[vertex].weight;
        }
    }

    return newP;
}

void bron_kerbosch_algorithm(BronKerbosch *bk, const ArrayOfVertices *C, ListOfVertices *P, ArrayOfVertices *S) {
    //P and S are empty
    //maximal clique above a threshold found
    if( (P->vertices.empty()) && (!S->filledPositions) && (C->totalWeight >= bk->minWeight) ) {
        vector<int> clique = generate_clique(bk, C);
        clq_set_add(bk->clqSet, (int)clique.size(), &clique[0], (int)C->totalWeight);
        if(C->totalWeight > bk->maxWeight)
            bk->maxWeight = C->totalWeight;
    }

    bk->it++;

    if(bk->it > bk->maxIt || P->vertices.empty())
        return;

    if(C->totalWeight + P->totalWeight >= bk->minWeight) {
        const size_t u = *(P->vertices.begin());
        vector<size_t> P_excluding_N_u;

        P_excluding_N_u.reserve(bk->vertices.size());
        exclude_neighbors_u(bk, P_excluding_N_u, P, u);

        for(const size_t v : P_excluding_N_u) {
            ArrayOfVertices *newC = create_new_C(bk, C, v);
            ArrayOfVertices *newS = create_new_S(bk, S, v);
            ListOfVertices *newP = create_new_P(bk, P, v);

            bron_kerbosch_algorithm(bk, newC, newP, newS);

            //freeing memory
            array_of_vertices_free(newC);
            array_of_vertices_free(newS);
            list_of_vertices_free(newP);

            //P = P \ {v}
            P->vertices.remove(v);
            P->totalWeight -= bk->vertices[v].weight;

            //S = S U {v}
            S->vertices[v/INT_SIZE] |= bk->mask[v%INT_SIZE];//adding v in S
            S->filledPositions++;
        }
    }
}

int bk_run(BronKerbosch *bk) {
    ArrayOfVertices *C = array_of_vertices_create((int)bk->vertices.size()/INT_SIZE + 1);
    ArrayOfVertices *S = array_of_vertices_create((int)bk->vertices.size()/INT_SIZE + 1);
    ListOfVertices *P = list_of_vertices_create();

    for(size_t i = 0; i < bk->vertices.size(); i++) {
        P->vertices.push_back(i);
        P->totalWeight += bk->vertices[i].weight;
    }

    P->vertices.sort(
            [bk](const size_t &a, const size_t &b) -> bool {
                if(bk->vertices[a].mdegree != bk->vertices[b].mdegree)
                    return bk->vertices[a].mdegree > bk->vertices[b].mdegree;
                return bk->vertices[a].degree > bk->vertices[b].degree;
            }
    );

    bron_kerbosch_algorithm(bk, C, P, S);

    array_of_vertices_free(C);
    list_of_vertices_free(P);
    array_of_vertices_free(S);

    return (bk->it > bk->maxIt);
}

const CliqueSet* bk_get_clq_set(const BronKerbosch *bk) {
    return bk->clqSet;
}

size_t bk_get_max_weight(const BronKerbosch *bk) {
    return bk->maxWeight;
}

void bk_set_min_weight(BronKerbosch *bk, size_t minWeight) {
    assert(minWeight >= 0);
    bk->minWeight = minWeight;
}

void bk_set_max_it(BronKerbosch *bk, size_t maxIt) {
    assert(maxIt > 0);
    bk->maxIt = maxIt;
}