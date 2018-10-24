#include <vector>
#include <list>
#include <limits>
#include <cassert>
#include <algorithm>
#include <cstdio>
#include "bron_kerbosch.h"

#define INT_SIZE (8*sizeof(int))

typedef struct {
    int id;
    int weight;
    int degree;
    int mdegree; //modified degree: degree of vertex + sum of degrees of adjacent vertices
} BKVertex;

typedef struct {
    std::list<int> vertices;
    int totalWeight;
} ListOfVertices;

typedef struct {
    size_t *vertices;
    int numVertices; //size of array "vertices"
    int totalWeight; //sum of weights of active vertices
    int filledPositions; //controls the number of active vertices
} ArrayOfVertices;

struct _BronKerbosch {
    //graph data
    const CGraph *cgraph;
    std::vector<BKVertex> vertices;
    int nVertices;

    //bk data
    size_t **bit;
    size_t *mask;

    //bk parameters
    int minWeight;
    int maxIt;

    //bk statistics
    int it;
    int maxWeight;
    int status;

    //bk result
    CliqueSet *clqSet;
};

ArrayOfVertices* array_of_vertices_create(int size) {
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
    const int cgSize = cgraph_size(cgraph);
    int *neighs = new int[cgSize];

    bk->vertices.reserve(cgSize);
    bk->bit = nullptr;
    bk->mask = nullptr;
    bk->clqSet = nullptr;

    for(int i = 0; i < cgSize; i++) {
        int realDegree = cgraph_degree(cgraph, i);

        if(realDegree == 0)
            continue;

        int check = cgraph_get_all_conflicting(cgraph, i, neighs, cgSize);
#ifdef DEBUG
        assert(check == realDegree);
#endif

        int mdegree = realDegree;
        for(int j = 0; j < realDegree; j++) {
#ifdef DEBUG
            assert(neighs[j] != i && neighs[j] >= 0 && neighs[j] < cgSize);
#endif
            mdegree += cgraph_degree(cgraph, neighs[j]);
        }

        aux.id = i;
        aux.weight = weights[i];
        aux.degree = realDegree;
        aux.mdegree = mdegree;
        bk->vertices.push_back(aux);
    }

    bk->nVertices = (int)bk->vertices.size();

    if(bk->nVertices > 0) {
        bk->cgraph = cgraph;
        bk->clqSet = clq_set_create();
        bk->status = 0L;
        bk->maxWeight = 0L;
        bk->it = 0L;
        bk->maxIt = (std::numeric_limits<int>::max() / 100);
        bk->minWeight = 0L;

        bk->mask = new size_t[INT_SIZE];
        bk->mask[0] = 1;
        for (size_t h = 1; h < INT_SIZE; h++)
            bk->mask[h] = bk->mask[h - 1] << 1U;

        bk->bit = new size_t*[bk->nVertices]();
        for (int i = 0; i < bk->nVertices; ++i)
            bk->bit[i] = new size_t[bk->nVertices / INT_SIZE + 1]();

        for (int v = 0; v < bk->nVertices; ++v) {
            for (int y = (v + 1); y < bk->nVertices; ++y) {
                if (cgraph_conflicting_nodes(bk->cgraph, bk->vertices[v].id, bk->vertices[y].id) == 1) {
                    bk->bit[y][v / INT_SIZE] |= bk->mask[v % INT_SIZE];
                    bk->bit[v][y / INT_SIZE] |= bk->mask[y % INT_SIZE];
                }
            }
        }
    }

    delete[] neighs;

    return bk;
}

void bk_free(BronKerbosch **_bk) {
    BronKerbosch *bk = *_bk;

    if(bk->clqSet)
        clq_set_free(&bk->clqSet);

    if(bk->bit) {
        for (int i = 0; i < bk->nVertices; i++)
            delete[] bk->bit[i];
        delete[] bk->bit;
    }

    if(bk->mask)
        delete[] bk->mask;

    delete bk;
    *_bk = nullptr;
}

std::vector<int> exclude_neighbors_u(const BronKerbosch *bk, const ListOfVertices *P, int u) {
    std::vector<int> P_excluding_N_u;
    P_excluding_N_u.reserve(bk->nVertices);
    for(const int &vertex : P->vertices) {
        if(!cgraph_conflicting_nodes(bk->cgraph, bk->vertices[u].id, bk->vertices[vertex].id)) {
            P_excluding_N_u.push_back(vertex);
        }
    }
    return P_excluding_N_u;
}

std::vector<int> generate_clique(const BronKerbosch *bk, const ArrayOfVertices *C) {
    size_t node, w = 0L, value;
    std::vector<int> nodes;

    nodes.reserve(C->numVertices);

    for(int t = 0; t < C->numVertices; t++) {
        node = INT_SIZE * t;
        value =  C->vertices[t];

        while(value > 1) {
            if(value % 2 == 1) {
                nodes.push_back(bk->vertices[node].id);
                w += bk->vertices[node].weight;
            }

            value = value / 2;
            node++;
        }

        if(value == 1) {
            nodes.push_back(bk->vertices[node].id);
            w += bk->vertices[node].weight;
        }
    }

#ifdef DEBUG
    assert(w == C->totalWeight);
#endif

    /*int n1[nodes.size()], n2[nodes.size()];
    if(!clq_validate(bk->cgraph, (int)nodes.size(), &nodes[0], n1, n2)) {
        fprintf(stderr, "Clique invalido!\n");
        for(const int &n : nodes)
            fprintf(stderr, "%d ", n);
        fprintf(stderr, "\n");
        exit(EXIT_FAILURE);
    }*/

    return nodes;
}

ArrayOfVertices* create_new_C(const BronKerbosch *bk, const ArrayOfVertices *C, const int v) {
    ArrayOfVertices* newC = array_of_vertices_create(bk->nVertices/INT_SIZE + 1);

    //copying elements of C
    std::copy(C->vertices, C->vertices + C->numVertices, newC->vertices);
    newC->numVertices = C->numVertices;
    newC->totalWeight = C->totalWeight;
    newC->filledPositions = C->filledPositions;

    //adding v in C
    newC->vertices[v/INT_SIZE] |= bk->mask[v%INT_SIZE];
    newC->filledPositions++;
    newC->totalWeight += bk->vertices[v].weight;

    return newC;
}

ArrayOfVertices* create_new_S(const BronKerbosch *bk, const ArrayOfVertices *S, const int v) {
    ArrayOfVertices *newS = array_of_vertices_create(bk->nVertices/INT_SIZE + 1);

    //newS = S intersection N(v)
    newS->filledPositions = 0L;
    for(int i = 0; i < S->numVertices; ++i) {
        newS->vertices[i] = S->vertices[i] & bk->bit[v][i];
        if(newS->vertices[i])
            newS->filledPositions++;
    }

    return newS;
}

ListOfVertices* create_new_P(const BronKerbosch *bk, const ListOfVertices *P, const int v) {
    ListOfVertices *newP = list_of_vertices_create();

    //newP = P intersection N(v)
    for(const int &vertex : P->vertices) {
        if(cgraph_conflicting_nodes(bk->cgraph, bk->vertices[v].id, bk->vertices[vertex].id)) {
            newP->vertices.push_back(vertex);
            newP->totalWeight += bk->vertices[vertex].weight;
        }
    }

    return newP;
}

void bron_kerbosch_algorithm(BronKerbosch *bk, const ArrayOfVertices *C, ListOfVertices *P, ArrayOfVertices *S) {
    //P and S are empty
    //maximal clique above a threshold found
    if( (P->vertices.empty()) && (S->filledPositions == 0) && (C->filledPositions > 0) &&
        (C->totalWeight >= bk->minWeight) ) {
        std::vector<int> clique = generate_clique(bk, C);
        clq_set_add(bk->clqSet, clique, C->totalWeight);
        if(C->totalWeight > bk->maxWeight)
            bk->maxWeight = C->totalWeight;
    }

    bk->it++;

    if(bk->it > bk->maxIt || P->vertices.empty())
        return;

    if(C->totalWeight + P->totalWeight >= bk->minWeight) {
        const int u = *(P->vertices.begin());
        const std::vector<int> P_excluding_N_u = exclude_neighbors_u(bk, P, u);

        for(const int &v : P_excluding_N_u) {
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

struct CompareMdegree {
    CompareMdegree(const BronKerbosch *bk) { this->bk = bk; }
    bool operator () (const int &i, const int &j) {
        if(bk->vertices[i].mdegree != bk->vertices[j].mdegree)
            return bk->vertices[i].mdegree >= bk->vertices[j].mdegree;
        return bk->vertices[i].degree > bk->vertices[j].degree;
    }
    const BronKerbosch *bk;
};

int bk_run(BronKerbosch *bk) {
    ArrayOfVertices *C = array_of_vertices_create(bk->nVertices/INT_SIZE + 1);
    ArrayOfVertices *S = array_of_vertices_create(bk->nVertices/INT_SIZE + 1);
    ListOfVertices *P = list_of_vertices_create();

    for(int i = 0; i < bk->nVertices; i++) {
        P->vertices.push_back(i);
        P->totalWeight += bk->vertices[i].weight;
    }

    P->vertices.sort(CompareMdegree(bk));

    bron_kerbosch_algorithm(bk, C, P, S);

    array_of_vertices_free(C);
    list_of_vertices_free(P);
    array_of_vertices_free(S);

    return (bk->it > bk->maxIt);
}

const CliqueSet* bk_get_clq_set(const BronKerbosch *bk) {
    return bk->clqSet;
}

int bk_get_max_weight(const BronKerbosch *bk) {
    return bk->maxWeight;
}

void bk_set_min_weight(BronKerbosch *bk, int minWeight) {
#ifdef DEBUG
    assert(minWeight >= 0);
#endif
    bk->minWeight = minWeight;
}

void bk_set_max_it(BronKerbosch *bk, int maxIt) {
#ifdef DEBUG
    assert(maxIt > 0);
#endif
    bk->maxIt = maxIt;
}