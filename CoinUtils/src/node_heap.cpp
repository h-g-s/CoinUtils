#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <utility>
#include "node_heap.h"

/* linux specific feature for debugging */
#ifdef linux
#include <execinfo.h>
void nh_print_trace(FILE *out, const char *file, int line);
#endif

// position of root node in vector
#define rootPos(node) (((node + 1) / 2) - 1)
// position of the first child node in vector
#define childPos(node) ((node * 2) + 1)

typedef struct {
  int node;
  int cost;
} NodeCost;

struct _NodeHeap {
  NodeCost *pq; // priority queue itself
  int *pos; // indicating each node where it is in pq
  int nodes;
  int INFTY;
};

/* private functions */
void nh_swap(NodeHeap *npq, const int i1, const int i2)
{
  NodeCost t = npq->pq[i1];
  npq->pq[i1] = npq->pq[i2];
  npq->pq[i2] = t;

  npq->pos[npq->pq[i1].node] = i1;
  npq->pos[npq->pq[i2].node] = i2;
}

void nh_down(NodeHeap *npq, const int index);

void nh_up(NodeHeap *npq, const int index);

NodeHeap *nh_create(const int nodes, const int infinity)
{
  NodeHeap *result = new NodeHeap;

  result->nodes = nodes;
  result->INFTY = infinity;
  result->pq = new NodeCost[nodes];
  result->pos = new int[nodes];
  nh_reset(result);

  return result;
}

void nh_update(NodeHeap *npq, const int node, const int cost)
{
  const int pos = npq->pos[node];

  if (cost > npq->pq[pos].cost) {
    fprintf(stderr, "\nERROR:\n");
#ifdef linux
    nh_print_trace(stderr, __FILE__, __LINE__);
#else
    fprintf(stderr, "\t%s:%d\n", __FILE__, __LINE__);
#endif
    fprintf(stderr, "\tmonotone heap only accepts decreasing values.\n");
    fprintf(stderr, "\tnode %d old cost: %d new cost: %d.\n", node, npq->pq[pos].cost, cost);
    fprintf(stderr, "\texiting.\n\n");
    exit(EXIT_FAILURE);
  }

  npq->pq[pos].cost = cost;
  nh_up(npq, pos);
}

void nh_down(NodeHeap *npq, const int index)
{
  int root = index;
  int child;

  while ((child = childPos(root)) < npq->nodes) {
    // child with the smallest cost
    if ((child + 1 < npq->nodes) && (npq->pq[child].cost > npq->pq[child + 1].cost)) {
      child++;
    }

    if (npq->pq[root].cost > npq->pq[child].cost) {
      nh_swap(npq, root, child);
      root = child;
    } else {
      break;
    }
  }
}

void nh_up(NodeHeap *npq, const int index)
{
  int root, child = index;

  while ((root = rootPos(child)) >= 0) {
    if (npq->pq[root].cost > npq->pq[child].cost) {
      nh_swap(npq, child, root);
      child = root;
    } else {
      return;
    }
  }
}

const int nh_remove_first(NodeHeap *npq, int *node)
{
  const int posLastNode = npq->nodes - 1;
  int cost = npq->pq[0].cost;

  (*node) = npq->pq[0].node;
  npq->pq[0] = npq->pq[posLastNode];
  npq->pq[posLastNode].cost = npq->INFTY;
  npq->pq[posLastNode].node = (*node);
  npq->pos[npq->pq[0].node] = 0;
  npq->pos[(*node)] = posLastNode;
  nh_down(npq, 0);

  return cost;
}

void nh_reset(NodeHeap *npq)
{
  for (int i = 0; i < npq->nodes; i++) {
    npq->pq[i].node = i;
    npq->pq[i].cost = npq->INFTY;
    npq->pos[i] = i;
  }
}

int nh_get_dist(NodeHeap *npq, const int node)
{
  return npq->pq[npq->pos[node]].cost;
}

void nh_free(NodeHeap **nh)
{
  delete[](*nh)->pos;
  delete[](*nh)->pq;
  delete (*nh);
  (*nh) = nullptr;
}

#ifdef linux
void nh_print_trace(FILE *out, const char *file, int line)
{
  const int max_depth = 100;
  int stack_depth;
  void *stack_addrs[max_depth];
  char **stack_strings;

  stack_depth = backtrace(stack_addrs, max_depth);
  stack_strings = backtrace_symbols(stack_addrs, stack_depth);

  fprintf(out, "Call stack from %s:%d:\n", file, line);

  for (size_t i = 1; i < stack_depth; i++) {
    fprintf(out, "    %s\n", stack_strings[i]);
  }
  free(stack_strings); // malloc()ed by backtrace_symbols
  fflush(out);
}
#endif
