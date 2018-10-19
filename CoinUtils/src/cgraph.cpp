#include <cctype>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <vector>
#include "cgraph.h"
#include "clique.h"
#include "vint_set.h"

#define LINE_SIZE 2048
#define MAX_NAME_SIZE 64
#define MIN_CLIQUE_ROW 256 /* mininum size for a row to be considered a clique row */

struct _CGraph {
    /* per node stored conflicts
       at position i there is a set with
       all node conflicts */
    IntSet **nodeConflicts;
    std::vector<int> *nodeCliques;  /* all cliques in which a node appears */
    int nodeSize; /* number of nodes considered */

    CliqueSet *clqSet;

    /* degree of each node */
    int *degree;

    int minDegree;
    int maxDegree;
    int lowDegree;

    int *origIdx; /* if it is a preprocessed graph,
                    indicates for each node i its original node */

    /* node names and weights are optional information */
    char **nodeNames;

    int *w;
};

typedef struct CGArc {
    int tail;
    int head;
} CGArc;

/** for every tail add all heads as neighbors. vector arcs should have capacit of nArcs+1 **/
void cgraph_add_vector_arcs(CGraph *cgraph, const std::vector<CGArc> &arcs, int *neighs);

CGraph *cgraph_clone(const CGraph *cg) {
    CGraph *clone = new CGraph;
    clone->nodeSize = cg->nodeSize;
    clone->nodeConflicts = new IntSet*[clone->nodeSize];
    clone->nodeCliques = new std::vector<int>[clone->nodeSize];

    for (int i = 0; i < clone->nodeSize; i++) {
        clone->nodeConflicts[i] = vint_set_clone(cg->nodeConflicts[i]);

        if (!cg->nodeCliques[i].empty())
            clone->nodeCliques[i] = cg->nodeCliques[i];
    }

    clone->clqSet = clq_set_clone(cg->clqSet);

    clone->degree = new int[clone->nodeSize];
    std::copy(cg->degree, cg->degree + cg->nodeSize, clone->degree);

    clone->minDegree = cg->minDegree;
    clone->maxDegree = cg->maxDegree;
    clone->lowDegree = cg->lowDegree;

    if (cg->origIdx) {
        clone->origIdx = new int[clone->nodeSize];
        std::copy(cg->origIdx, cg->origIdx + cg->nodeSize, clone->origIdx);
    } else
        clone->origIdx = nullptr;

    if (cg->nodeNames) {
        clone->nodeNames = new char *[clone->nodeSize];
        clone->nodeNames[0] = new char[MAX_NAME_SIZE * clone->nodeSize];
        for (int i = 1; i < clone->nodeSize; i++)
            clone->nodeNames[i] = clone->nodeNames[i - 1] + MAX_NAME_SIZE;
        for (int i = 0; i < clone->nodeSize; i++)
            strncpy(clone->nodeNames[i], cg->nodeNames[i], MAX_NAME_SIZE);
    } else
        clone->nodeNames = nullptr;

    if (cg->w) {
        clone->w = new int[cg->nodeSize];
        std::copy(cg->w, cg->w + cg->nodeSize, clone->w);
    } else
        clone->w = nullptr;

    return clone;
}

CGraph *cgraph_create(int numColumns) {
    CGraph *result = new CGraph;

    result->nodeSize = numColumns;
    result->nodeConflicts = new IntSet*[result->nodeSize];
    result->nodeCliques = new std::vector<int>[result->nodeSize];

    for (int i = 0; (i < result->nodeSize); i++) {
        result->nodeConflicts[i] = vint_set_create();
    }

    result->degree = new int[numColumns]();
    result->clqSet = clq_set_create();
    result->nodeNames = nullptr;
    result->w = nullptr;
    result->origIdx = nullptr;
    result->lowDegree = 15;

    return result;
}

void cgraph_add_node_conflicts(CGraph *cgraph, const int node, const int *conflicts, const int size) {
#ifdef DEBUG
    assert(node >= 0 && node < cgraph->nodeSize);
#endif

    if(clq_set_number_of_cliques(cgraph->clqSet) > 0) {
        /* adding conflicts (node, conflicts[i]) and computing a estimated degree */
        for (int i = 0; i < size; i++) {
#ifdef DEBUG
            assert(conflicts[i] >= 0 && conflicts[i] < cgraph->nodeSize);
#endif
            vint_set_add(cgraph->nodeConflicts[conflicts[i]], &node, 1);
            cgraph->degree[conflicts[i]]++;
        }
        cgraph->degree[node] += size;
    } else {
        /* adding conflicts and computing the correct degree */
        for (int i = 0; i < size; i++) {
#ifdef DEBUG
            assert(conflicts[i] >= 0 && conflicts[i] < cgraph->nodeSize);
#endif
            if (!cgraph_conflicting_nodes(cgraph, node, conflicts[i]))
                cgraph->degree[node]++;
            if (!cgraph_conflicting_nodes(cgraph, conflicts[i], node))
                cgraph->degree[conflicts[i]]++;
            vint_set_add(cgraph->nodeConflicts[conflicts[i]], &node, 1);
        }
    }

    vint_set_add(cgraph->nodeConflicts[node], conflicts, size);
}

void cgraph_add_node_conflicts_no_sim(CGraph *cgraph, const int node, const int *conflicts, const int size) {
    if (size <= 0)
        return;

#ifdef DEBUG
    assert(node >= 0 && node < cgraph->nodeSize);
#endif

    if(clq_set_number_of_cliques(cgraph->clqSet) > 0) {
        /* computing a estimated degree */
        cgraph->degree[node] += size;
    } else {
        /* computing the correct degree */
        for (int i = 0; i < size; i++) {
#ifdef DEBUG
            assert(conflicts[i] >= 0 && conflicts[i] < cgraph->nodeSize);
#endif
            if (!cgraph_conflicting_nodes(cgraph, node, conflicts[i]))
                cgraph->degree[node]++;
        }
    }

    vint_set_add(cgraph->nodeConflicts[node], conflicts, size);
}

void cgraph_add_clique(CGraph *cgraph, const int *idxs, const int size) {
    if (!(clq_set_add(cgraph->clqSet, idxs, size, 0)))
        return;

    const int idxClique = clq_set_number_of_cliques(cgraph->clqSet) - 1;
    /* computes a estimated degree */
    for (int i = 0; i < size; i++) {
        cgraph->degree[idxs[i]] += (size - 1);
        /* this clique will be related to all nodes inside it */
        cgraph->nodeCliques[idxs[i]].push_back(idxClique);
    }
}

int cgraph_size(const CGraph *cgraph) {
    return cgraph->nodeSize;
}

int cgraph_conflicting_nodes(const CGraph *cgraph, const int i, const int j) {
    if (i == j)
        return 0;

#ifdef DEBUG
    assert(i >= 0);
    assert(j >= 0);
    assert(i < cgraph_size(cgraph));
#endif
    if (vint_set_size(cgraph->nodeConflicts[i]) == 0)
        goto FIND_IN_CLIQUES;

    if (vint_set_find(cgraph->nodeConflicts[i], j) != -1)
        return 1;

    FIND_IN_CLIQUES:
    for (const int idx : cgraph->nodeCliques[i]) {
        if (clq_set_clique_has_element(cgraph->clqSet, idx, j)) {
#ifdef DEBUG
            if (!clq_set_clique_has_element(cgraph->clqSet, idx, i)) {
                fprintf(stderr, "Error:  node %d does not appears in clique %d\n", i, idx);
            }
            assert(clq_set_clique_has_element(cgraph->clqSet, idx, i));
            assert(clq_set_clique_has_element(cgraph->clqSet, idx, j));
#endif
            return 1;
        }
    }

    return 0;
}

int cgraph_get_all_conflicting(const CGraph *cgraph, int node, int *neighs, int maxSize) {
    const int size = cgraph->nodeSize;
    char *iv = new char[size]();
    IntSet *isnc = cgraph->nodeConflicts[node];
    const std::vector<int> &el = vint_set_get_elements(isnc);
    int nConfs = (int)el.size();
    int clqSize = 0;

    if (nConfs > maxSize) {
        fprintf(stderr, "ERROR: cgraph_get_all_conflicting:: Not enough space specified in maxSize.\n");
        fprintf(stderr,
                "Working with node %d, which appears in %lu cliques.\n", node, cgraph->nodeCliques[node].size());
        fprintf(stderr, "at: %s:%d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    iv[node] = 1;
    for (int i = 0; i < nConfs; i++) {
        iv[el[i]] = 1;
        neighs[i] = el[i];
    }

    /* now filling information from cliques, i.e., implicitly stored conflicts */
    for (const int idxClique : cgraph->nodeCliques[node]) {
        const int *clqEl = clq_set_clique_elements(cgraph->clqSet, idxClique);
        for (int j = 0; j < clq_set_clique_size(cgraph->clqSet, idxClique); j++) {
            if (!iv[clqEl[j]]) {
                iv[clqEl[j]] = 1;
                neighs[nConfs++] = clqEl[j];

                if (nConfs > maxSize) {
                    fprintf(stderr, "ERROR: cgraph_get_all_conflicting:: Not enough space specified in maxSize.\n");
                    fprintf(stderr,
                            "Working with node %d, which appears in %lu cliques. When adding clique %d size %d. Result %d. MaxSize %d.\n",
                            node, cgraph->nodeCliques[node].size(), idxClique, clqSize, nConfs, maxSize);
                    fprintf(stderr, "at: %s:%d\n", __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    delete[] iv;

    return nConfs;
}

void cgraph_save(CGraph *cgraph, const char *fileName) {
    int *w = cgraph->w;

    FILE *f = fopen(fileName, "w");
    if (!f) {
        fprintf(stderr, "Could not open file %s.", &(fileName[0]));
        exit(EXIT_FAILURE);
    }

    /* counting edges */
    int nEdges = 0;
    int nodes = cgraph_size(cgraph);

    for (int i = 0; i < cgraph_size(cgraph); i++) {
        IntSet *is = cgraph->nodeConflicts[i];
        const std::vector<int> &el = vint_set_get_elements(is);
        const size_t nEl = el.size();
        for (size_t j = 0; j < nEl; j++)
            if (el[j] > i)
                nEdges++;
    }

    fprintf(f, "p edges %d %d\n", cgraph_size(cgraph), nEdges);

    for (int i = 0; i < cgraph_size(cgraph); i++) {
        IntSet *is = cgraph->nodeConflicts[i];
        const std::vector<int> &el = vint_set_get_elements(is);
        const int nEl = (int)el.size();
        for (int j = 0; j < nEl; j++)
            if (el[j] > i)
                fprintf(f, "e %d %d\n", i + 1, el[j] + 1);
    }

    if (w) {
        int i;
        for (i = 0; (i < nodes); ++i)
            fprintf(f, "w %d %d\n", i + 1, w[i]);
    }

    for (int i = 0; i < clq_set_number_of_cliques(cgraph->clqSet); i++) {
        fprintf(f, "c %d\n", clq_set_clique_size(cgraph->clqSet, i));
        const int *element = clq_set_clique_elements(cgraph->clqSet, i);
        for (int j = 0; j < clq_set_clique_size(cgraph->clqSet, i); j++)
            fprintf(f, "%d\n", element[j] + 1);
    }

    if (cgraph->nodeNames) {
        int i;
        for (i = 0; (i < cgraph_size(cgraph)); ++i) {
            if (cgraph_get_node_name(cgraph, i)) {
                char name[256];
                strcpy(name, cgraph_get_node_name(cgraph, i));
                if (name[strlen(name) - 1] == '\n')
                    name[strlen(name) - 1] = '\0';
                fprintf(f, "n %d %s\n", i + 1, name);
            }
        }
    }

    if (cgraph->origIdx) {
        int i;
        for (i = 0; (i < cgraph_size(cgraph)); ++i)
            fprintf(f, "o %d %d\n", i + 1, cgraph->origIdx[i] + 1);
    }

    fclose(f);
}

bool compArc(const CGArc &a1, const CGArc &a2) {
    if(a1.tail != a2.tail)
        return a1.tail < a2.tail;
    return a1.head < a2.head;
}

void cgraph_add_vector_arcs(CGraph *cgraph, const std::vector<CGArc> &arcs, int *neighs) {
    if (arcs.empty())
        return;

    neighs[0] = arcs[0].head;
    int nNeighs = 1;
    int lastNode = arcs[0].tail;

    for (size_t i = 1; i < arcs.size() + 1; i++) {
        if (lastNode != arcs[i].tail) {
#ifdef DEBUG
            assert(lastNode >= 0);
            assert(lastNode < cgraph_size(cgraph));
#endif
            cgraph_add_node_conflicts_no_sim(cgraph, lastNode, neighs, nNeighs);
            lastNode = arcs[i].tail;
            nNeighs = 0;
            if (i == arcs.size())
                break;
        }
#ifdef DEBUG
        assert(arcs[i].head >= 0);
        assert(arcs[i].head < cgraph_size(cgraph));
#endif
        neighs[nNeighs++] = arcs[i].head;
    }
}

CGraph *cgraph_load(const char *fileName) {
    int nodes = -1;
    int edges = -1;
    char recomputeDegree = 0;

    FILE *f = fopen(fileName, "r");
    if (!f) {
        fprintf(stderr, "Could not open file %s.\n", fileName);
        exit(EXIT_FAILURE);
    }

    char line[LINE_SIZE];
    while (fgets(line, LINE_SIZE, f)) {
        if ((line[0] == 'p') || (line[0] == 'P')) {
            char *ptr = line;
            char *end = ptr + LINE_SIZE;
            while ((!isdigit(*ptr)) && (ptr < end))
                ++ptr;
            sscanf(ptr, "%d %d", &nodes, &edges);
            break;
        }
    }

    if (nodes == -1) {
        fprintf(stderr, "Invalid file format.\n");
        exit(EXIT_FAILURE);
    }

    std::vector<CGArc> arcs; arcs.reserve(8192);
    CGArc tmpArc;
    int *clqElement = new int [nodes];
    CliqueSet *clqSet = clq_set_create();
    int weightsLoaded = 0; /* lineNumber = 1; */

    char **names = nullptr;
    int *tw = nullptr;

    while (fgets(line, LINE_SIZE, f)) {
        if ((line[0] == 'e') || (line[0] == 'E')) {
            int tail, head;
            sscanf(line + 2, "%d %d", &tail, &head);
            --tail;
            --head;
#ifdef DEBUG
            assert(tail >= 0);
            assert(tail < nodes);
            assert(head >= 0);
            assert(head < nodes);
#endif
            if (tail > head) {
                int t = tail;
                tail = head;
                head = t;
            }

            tmpArc.tail = tail;
            tmpArc.head = head;
            arcs.push_back(tmpArc);
        } else {
            if (((line[0] == 'w') || (line[0] == 'W'))) {
                if (!tw) {
                    tw = new int[nodes];
                }

                int nodeIdx, nodeWeight;
                sscanf(line + 2, "%d %d", &nodeIdx, &nodeWeight);
                nodeIdx--;
#ifdef DEBUG
                assert((nodeIdx >= 0) && (nodeIdx < nodes));
#endif
                tw[nodeIdx] = nodeWeight;

                weightsLoaded = 1;
            } else {
                if (((line[0] == 'c') || (line[0] == 'C'))) {
#ifdef DEBUG
                    assert(weightsLoaded);
#endif
                    int clqSize;
                    int totalWeight = 0;
                    sscanf(line + 2, "%d", &clqSize);
                    for (int i = 0; i < clqSize; i++) {
                        if (fscanf(f, "%d", &clqElement[i]) != 1) {
                            fprintf(stderr, "\nERROR: Expected clique element not found.\n");
                            exit(EXIT_FAILURE);
                        }
                        clqElement[i]--;
#ifdef DEBUG
                        assert(clqElement[i] >= 0);
                        assert(clqElement[i] < nodes);
#endif
                        totalWeight += tw[clqElement[i]];
                    }
                    clq_set_add(clqSet,clqElement, clqSize, totalWeight);
                    recomputeDegree = 1;
                } /* reading clique */
                else {
                    /* reading column names */
                    if (((line[0] == 'n') || (line[0] == 'N'))) {
                        if (names == nullptr) {
                            names = new char *[nodes];
                            names[0] = new char[MAX_NAME_SIZE * nodes];
                            int i;
                            for (i = 1; (i < nodes); ++i)
                                names[i] = names[i - 1] + MAX_NAME_SIZE;
                        }

                        int idxNode;
                        sscanf(line + 2, "%d", &idxNode);
                        int digits = 1;
                        idxNode--;
                        if (idxNode != 0)
                            digits = (int) (log10(((double) idxNode)) + 1.0);

                        strncpy(names[idxNode], line + 2 + digits + 1, MAX_NAME_SIZE);
                        size_t len = strlen(names[idxNode]);
                        if (len > 0)
                            if (names[idxNode][len - 1] == '\n')
                                names[idxNode][len - 1] = '\0';
                    }
                } /* not clique */
            } /* not weight */
        } /* not edge */
    }

    std::sort(arcs.begin(), arcs.end(), compArc);

    /* removing repetitions */
    int nRep = 0;
    CGArc lastArc = {-1, -1};
    for (auto &arc : arcs)
        if (arc.head == lastArc.head && arc.tail == lastArc.tail) {
            arc.tail = std::numeric_limits<int>::max();
            arc.head = std::numeric_limits<int>::max();
            nRep++;
        } else
            lastArc = arc;

    if (nRep)
        std::sort(arcs.begin(), arcs.end(), compArc);

    arcs.resize(arcs.size() - nRep);

    CGraph *cgraph = cgraph_create(nodes);

    if (names)
        cgraph->nodeNames = names;
    if (tw) {
        cgraph->w = new int[nodes];
        memcpy(cgraph->w, tw, sizeof(int) * nodes);
        delete[] tw;
    }

    int *neighs = new int[arcs.size()];

    arcs.rbegin()->tail = std::numeric_limits<int>::max();
    arcs.rbegin()->head = std::numeric_limits<int>::max();

    cgraph_add_vector_arcs(cgraph, arcs, neighs);

    /* inverting head and tail */
    for (auto &arc : arcs) {
        const int t = arc.tail;
        arc.tail = arc.head;
        arc.head = t;
    }

    std::sort(arcs.begin(), arcs.end(), compArc);

    cgraph_add_vector_arcs(cgraph, arcs, neighs);

    for (int i = 0; i < clq_set_number_of_cliques(clqSet); i++) {
        cgraph_add_clique(cgraph, clq_set_clique_elements(clqSet, i), clq_set_clique_size(clqSet, i));
    }

    fclose(f);
    delete[] neighs;
    delete[] clqElement;
    clq_set_free(&clqSet);

    if (recomputeDegree)
        cgraph_recompute_degree(cgraph);
    else
        cgraph_update_min_max_degree(cgraph);

    return cgraph;
}

void cgraph_print(CGraph *cgraph, const int *w) {
    int *neighs = new int[cgraph_size(cgraph)];

    for (int i = 0; i < cgraph_size(cgraph); i++) {
        printf("[%d] ", i + 1);
        int n = cgraph_get_all_conflicting(cgraph, i, neighs, cgraph_size(cgraph));
        for (int j = 0; j < n; j++)
            printf("%d ", neighs[j] + 1);
        printf("\n");
    }

    if (w) {
        for (int i = 0; i < cgraph_size(cgraph); i++)
            printf("w[%d] %d\n", i + 1, w[i]);
    }

    delete[] neighs;
}

void cgraph_update_min_max_degree(CGraph *cgraph) {
    cgraph->minDegree = std::numeric_limits<int>::max();
    cgraph->maxDegree = 0;
    for (int i = 0; (i < cgraph_size(cgraph)); ++i) {
        cgraph->minDegree = std::min(cgraph->minDegree, cgraph->degree[i]);
        cgraph->maxDegree = std::max(cgraph->maxDegree, cgraph->degree[i]);
    }
}

int cgraph_degree(const CGraph *cgraph, const int node) {
    return cgraph->degree[node];
}

int cgraph_min_degree(const CGraph *cgraph) {
    return cgraph->minDegree;
}

int cgraph_max_degree(const CGraph *cgraph) {
    return cgraph->maxDegree;
}

void cgraph_free(CGraph **cgraph) {
    for(int i = 0; i < (*cgraph)->nodeSize; i++) {
        vint_set_free(&(*cgraph)->nodeConflicts[i]);
    }


    delete[] (*cgraph)->nodeConflicts;
    delete[] (*cgraph)->nodeCliques;
    delete[] (*cgraph)->degree;

    if ((*cgraph)->origIdx)
        delete[] (*cgraph)->origIdx;

    if ((*cgraph)->nodeNames) {
        delete[] (*cgraph)->nodeNames[0];
        delete (*cgraph)->nodeNames;
        (*cgraph)->nodeNames = nullptr;
    }

    if ((*cgraph)->w)
        delete[] (*cgraph)->w;

    clq_set_free(&((*cgraph)->clqSet));

    delete (*cgraph);

    (*cgraph) = nullptr;
}

void cgraph_load_dimensions(const char *fileName, int *nodes, int *edges) {
    FILE *f = fopen(fileName, "r");
    if (!f) {
        fprintf(stderr, "Could not open file %s\n", &(fileName[0]));
        exit(EXIT_FAILURE);
    }

    *nodes = -1;
    *edges = -1;

    char line[LINE_SIZE];
    while (fgets(line, LINE_SIZE, f)) {
        if ((line[0] == 'p') || (line[0] == 'P')) {
            const char *p;
            if (((p = strstr(line, "edge")) != nullptr)) {
                int shift = 7;
                if (p[5] == 's')
                    ++shift;
                if (sscanf(line + shift, "%d %d", nodes, edges) != 2) {
                    fprintf(stderr, "Error reading dimensions.\n");
                    exit(EXIT_FAILURE);
                }
            }

        }
    }

    fclose(f);
}

CGraph *cgraph_create_induced_subgraph( const CGraph *cgraph, const int *idxs, const int n ) {
    const int nOrig = cgraph_size(cgraph);
    char recomputeDegree = 0;
    int *ppIdx = new int[nOrig];
    std::fill(ppIdx, ppIdx + nOrig, -1);

    CGraph *result = cgraph_create(n);
    result->origIdx = new int[n];
    int *neighs = new int[n];

    int last = 0;
    for(int i = 0; i < n; i++) {
        const int origIdx = idxs[i];
        ppIdx[origIdx] = last;
        result->origIdx[last++] = origIdx;
    }

    if (cgraph->w)
        result->w = new int[n];

    if (cgraph->nodeNames) {
        result->nodeNames = new char*[n];
        result->nodeNames[0] = new char[MAX_NAME_SIZE * n];
        for (int i = 1; i < n; i++)
            result->nodeNames[i] = result->nodeNames[i - 1] + MAX_NAME_SIZE;
    }

    /* filling new conflicts */
    for (int i = 0; i < n; i++) {
        IntSet *nodeConflicts = cgraph->nodeConflicts[idxs[i]];
        const std::vector<int> &el = vint_set_get_elements(nodeConflicts);
        const int nEl = (int)el.size();
        int nNeighs = 0;
#ifdef DEBUG
        assert(idxs[i] >= 0 && idxs[i] < nOrig);
        assert(ppIdx[idxs[i]] >= 0 && ppIdx[idxs[i]] < n);
#endif
        for (int j = 0; j < nEl; j++) {
            const int origIdx = el[j];
            if (ppIdx[origIdx] != -1) {
                neighs[nNeighs++] = ppIdx[origIdx];
#ifdef DEBUG
                assert(ppIdx[origIdx] >= 0 && ppIdx[origIdx] < n);
#endif
            }
        }

        cgraph_add_node_conflicts_no_sim(result, ppIdx[idxs[i]], neighs, nNeighs);

        if (cgraph->w)
            result->w[ppIdx[idxs[i]]] = cgraph->w[idxs[i]];
        if (cgraph->nodeNames)
            strncpy(result->nodeNames[ppIdx[idxs[i]]], cgraph->nodeNames[idxs[i]], MAX_NAME_SIZE);
    }

    /* filling new cliques */
    const int nCliques = clq_set_number_of_cliques(cgraph->clqSet);
    for (int i = 0; i < nCliques; i++) {
        const int nEl = clq_set_clique_size( cgraph->clqSet, i );
        const int *el = clq_set_clique_elements( cgraph->clqSet, i );
        int newClqSize = 0;

        for (int j = 0; j < nEl; j++){
            if (ppIdx[el[j]] != -1) {
                neighs[newClqSize++] = ppIdx[el[j]];
#ifdef DEBUG
                assert(ppIdx[el[j]] >= 0 && ppIdx[el[j]] < n);
#endif
            }
        }

        if (newClqSize > 1) {
            if (newClqSize < MIN_CLIQUE_ROW)
                /* not adding as clique anymore */
                cgraph_add_clique_as_normal_conflicts(result, neighs, newClqSize);
            else
                cgraph_add_clique(result, neighs, newClqSize);
            recomputeDegree = 1;
        }
    }

    delete[] neighs;
    delete[] ppIdx;

    if (recomputeDegree)
        cgraph_recompute_degree(result);
    else
        cgraph_update_min_max_degree(result);

    return result;
}

int cgraph_get_node_weight(const CGraph *cgraph, int node) {
    if (cgraph->w)
        return cgraph->w[node];

    return 0;
}

void cgraph_set_node_weight(CGraph *cgraph, int node, int weight) {
    if (!cgraph->w)
        cgraph->w = new int[cgraph_size(cgraph)];

    cgraph->w[node] = weight;
}

const char *cgraph_get_node_name(const CGraph *cgraph, int node) {
    if (cgraph->nodeNames)
        return cgraph->nodeNames[node];

    return nullptr;
}

void cgraph_set_node_name(CGraph *cgraph, int node, const char *name) {
    if (!cgraph->nodeNames) {
        cgraph->nodeNames = new char*[cgraph_size(cgraph)];
        cgraph->nodeNames[0] = new char[MAX_NAME_SIZE * cgraph_size(cgraph)];
        for (int i = 1; (i < cgraph_size(cgraph)); ++i)
            cgraph->nodeNames[i] = cgraph->nodeNames[i - 1] + MAX_NAME_SIZE;
    }

    strncpy(cgraph->nodeNames[node], name, MAX_NAME_SIZE);
}

const int *cgraph_get_node_weights(const CGraph *cgraph) {
    return cgraph->w;
}

int cgraph_weight(const double w) {
    return (int) (w * 1000.0);
}

int cgraph_get_original_node_index(const CGraph *cgraph, const int node) {
    if (cgraph->origIdx)
        return cgraph->origIdx[node];

    return -1;
}


#ifdef DEBUG

void cgraph_check_node_cliques(const CGraph *cgraph) {
    /* all nodes */
    for (int i = 0; i < cgraph_size(cgraph); i++) {
        for (const int idx : cgraph->nodeCliques[i]) {
            if (!(clq_set_clique_has_element(cgraph->clqSet, idx, i))) {
                printf("\nnode %d should appear on clique %d but does not\n", i, idx);
                fflush(stdout);
            }
            assert(clq_set_clique_has_element(cgraph->clqSet, idx, i));
        }
    }

    /* printf("Information in node cliques indicate cliques which really have the node.\n"); */

    const int nc = clq_set_number_of_cliques(cgraph->clqSet);
    for (int i = 0; i < nc; i++) {
        const int *el = clq_set_clique_elements(cgraph->clqSet, i);
        for (int j = 0; j < clq_set_clique_size(cgraph->clqSet, i); j++) {
            const int currNode = el[j];
            /* this must appear in node clique */
            int l;
            const int nn = (int)cgraph->nodeCliques[currNode].size();

            for (l = 0; l < nn; l++)
                if (cgraph->nodeCliques[currNode][l] == i)
                    break;

            if (l == nn) {
                fprintf(stderr,
                        "ERROR: in clique %d node %d appears but this clique does not appears in its node list.\n", i,
                        currNode);
                fflush(stdout);
                fflush(stderr);
            }

            assert(l < nn);
        }
    }
}

void cgraph_check_neighs(const CGraph *cgraph) {
    const int n = cgraph_size(cgraph);
    int nNeighs;
    int *neighs = new int[n];

    for (int i = 0; i < n; i++) {
        nNeighs = cgraph_get_all_conflicting(cgraph, i, neighs, n);

        /* computed number of neighs */
        int cn = 0;
        for (int j = 0; j < n; j++)
            if (cgraph_conflicting_nodes(cgraph, i, j))
                cn++;

        assert(cn == nNeighs);

        for (int j = 0; j < n; j++)
            if (cgraph_conflicting_nodes(cgraph, i, j))
                assert(std::binary_search(neighs, neighs + nNeighs, j));
            else
                assert(!std::binary_search(neighs, neighs + nNeighs, j));
    }

    delete[] neighs;
}
#endif

const int *cgraph_get_original_node_indexes(const CGraph *cgraph) {
    return cgraph->origIdx;
}

void cgraph_print_summary(CGraph *cgraph, const char *name) {
    cgraph_update_min_max_degree(cgraph);
    int nodes = cgraph_size(cgraph);
    int nSmallDegree = 0;
    int edges = 0;

    for (int i = 0; i < nodes; i++) {
        edges += cgraph_degree(cgraph, i);
        if (cgraph_degree(cgraph, i) <= cgraph->lowDegree)
            nSmallDegree++;
    }

    printf("%s : n %d e %d - minD %d maxD %d smallD %d\n", name, cgraph_size(cgraph), edges, cgraph_min_degree(cgraph),
           cgraph_max_degree(cgraph), nSmallDegree);
}

void cgraph_set_low_degree(CGraph *cgraph, const int lowDegree) {
    cgraph->lowDegree = lowDegree;
}

void cgraph_add_clique_as_normal_conflicts(CGraph *cgraph, const int *idxs, const int size) {
    const int nconflicts = size - 1;
    int *confs = new int[size];
    const int last = idxs[size-1];

    std::copy(idxs, idxs + size, confs);

    for (int i = 0; i < size; i++) {
        const int node = confs[i];
        confs[i] = last;
        confs[nconflicts] = node;
        cgraph_add_node_conflicts_no_sim(cgraph, node, confs, nconflicts);
        confs[i] = node;
        confs[nconflicts] = last;
    }

    delete[] confs;
}

void cgraph_recompute_degree(CGraph *cgraph) {
    const int size = cgraph->nodeSize;
    char *iv = new char[size];

    cgraph->minDegree = std::numeric_limits<int>::max();
    cgraph->maxDegree = 0;
    memset(cgraph->degree, 0, sizeof(int) * size);

    for (int i = 0; i < size; i++) {
        memset(iv, 0, sizeof(char) * size);

        IntSet *isnc = cgraph->nodeConflicts[i];
        const std::vector<int> &el = vint_set_get_elements(isnc);
        const int nEl = (int)el.size();
        //individual conflicts
        cgraph->degree[i] += nEl;
        for (int j = 0; j < nEl; j++) {
            iv[el[j]] = 1;
        }

        //conflicts stored as cliques
        for (const int idxClique : cgraph->nodeCliques[i]) {
            const int* clqEl = clq_set_clique_elements(cgraph->clqSet, idxClique);
            for(int k = 0; k < clq_set_clique_size(cgraph->clqSet, idxClique); k++)
                if (!iv[clqEl[k]] && clqEl[k] != i) {
                    iv[clqEl[k]] = 1;
                    cgraph->degree[i]++;
                }
        }

        //updating min and max degree
        cgraph->minDegree = std::min(cgraph->minDegree, cgraph->degree[i]);
        cgraph->maxDegree = std::max(cgraph->maxDegree, cgraph->degree[i]);
    }

    delete[] iv;
}

typedef struct {
    int node;
    double cost;
} NodeCost;

bool compare_node_costs(const NodeCost &n1, const NodeCost &n2) {
    if(fabs(n1.cost - n2.cost) > 0.00001)
        return n1.cost < n2.cost;
    return n1.node < n2.node;
}

int cgraph_get_best_n_neighbors(const CGraph *cgraph, int node, const double *costs, int *neighs, int maxSize) {
    std::vector<NodeCost> candidates;
    candidates.reserve(vint_set_size(cgraph->nodeConflicts[node]) + (cgraph->nodeCliques[node].size() * 2));

#ifdef DEBUG
    assert( node >= 0 );
    assert( node < cgraph_size(cgraph) );
#endif

    NodeCost tmp;

    /* normal conflicts */
    IntSet *nodeConflicts = cgraph->nodeConflicts[node];
    const std::vector<int> &el = vint_set_get_elements(nodeConflicts);
    const int nEl = (int)el.size();
    for (int i = 0; i < nEl; i++) {
#ifdef DEBUG
        assert( (el[i]>=0) );
        assert( (el[i]<cgraph_size(cgraph)) );
#endif
        tmp.node = el[i];
        tmp.cost = costs[el[i]];
        candidates.push_back(tmp);
    }

    /* conflicts stored in cliques */
    const CliqueSet *clqSet = cgraph->clqSet;
    for (const int clique : cgraph->nodeCliques[node]) {
        const int *el = clq_set_clique_elements(clqSet, clique);
        for (int j = 0; j < clq_set_clique_size(clqSet, clique); j++) {
#ifdef DEBUG
            assert( (el[j]>=0) );
            assert( (el[j]<cgraph_size(cgraph)) );
#endif
            tmp.node = el[j];
            tmp.cost = costs[el[j]];
            candidates.push_back(tmp);
        }
    }

    int nNeighs = 0;

    if(candidates.size() < maxSize) {
        nNeighs = (int)candidates.size();
        std::sort(candidates.begin(), candidates.end(), compare_node_costs);
    }
    else {
        nNeighs = maxSize;
        std::partial_sort(candidates.begin(), candidates.begin() + maxSize, candidates.end(), compare_node_costs);
    }

    for(int i = 0; i < nNeighs; i++)
        neighs[i] = candidates[i].node;

    return nNeighs;
}
