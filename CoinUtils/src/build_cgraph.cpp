#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include <cassert>
#include <cstdio>
#include "build_cgraph.h"

#define EPS 1e-8
#define MIN_CLIQUE_ROW 256 /* mininum size for a row to be considered a clique row */
#define DBL_EQUAL(v1, v2) (fabs(v1 - v2) < EPS)

const double LARGE_CONST = (std::numeric_limits< double >::max() / 10.0);

bool sort_sec_pair(const std::pair< int, int > &left, const std::pair< int, int > &right)
{
  if (left.second != right.second)
    return left.second < right.second;
  return left.first < right.first;
}

bool sort_columns(const std::pair< int, double > &left, const std::pair< int, double > &right)
{
  if (fabs(left.second - right.second) > EPS)
    return (left.second < right.second);
  return left.first < right.first;
}

void processClique(CGraph *cgraph, const int *idxs, const int size, char *recomputeDegree);

/* Returns the first position of columns which the lower bound for LHS (considering activation of variables) is greater than rhs */
/* colStart=initial position for search in columns, colEnd=last position for search in columns */
/* partialLHS = LHS calculated with only one variable */
int binary_search(const std::vector< std::pair< int, double > > &columns, double partialLHS, double rhs, int colStart,
  int colEnd);

/* Returns the first position of columns which the lower bound for LHS (considering deactivation of variables) is greater than rhs */
/* colStart=initial position for search in columns, colEnd=last position for search in columns */
/* partialLHS = LHS calculated with only one variable */
int binary_search_complement(const std::vector< std::pair< int, double > > &columns, double partialLHS, double rhs,
  int colStart,
  int colEnd);

/* Searches for cliques involving the activation of variables in this constraint. */
void cliqueDetection(CGraph *cgraph, const std::vector< std::pair< int, double > > &columns, double sumNegCoefs, double rhs,
  char *recomputeDegree);

/* Searches for cliques involving the complement of variables in this constraint. */
void cliqueComplementDetection(CGraph *cgraph, const std::vector< std::pair< int, double > > &columns, double sumNegCoefs,
  double rhs, char *recomputeDegree);

/* Searches for cliques involving variables and complements of variables in this constraint. */
void mixedCliqueDetection(CGraph *cgraph, const std::vector< std::pair< int, double > > &columns, double sumNegCoefs,
  double rhs,
  std::vector< std::vector< int > > &cvec);

void pairwiseAnalysis(CGraph *cgraph, const std::vector< std::pair< int, double > > &columns, const double sumNegCoefs,
  const double rhs, std::vector< std::vector< int > > &cvec)
{
  const int nElements = (int)columns.size();
  const int cgSize = cgraph_size(cgraph);
#ifdef DEBUG
  assert(cgSize % 2 == 0);
#endif
  int nCols = cgSize / 2;

  for (int j1 = 0; j1 < nElements; j1++) {
    const int cidx1 = columns[j1].first;
    const double coef1 = columns[j1].second;

    for (int j2 = j1 + 1; j2 < nElements; j2++) {
      const int cidx2 = columns[j2].first;
      const double coef2 = columns[j2].second;
      const double negDiscount = sumNegCoefs - std::min(0.0, coef1) - std::min(0.0, coef2);

      if (coef1 + coef2 + negDiscount > rhs + EPS) {
        cvec[cidx1].push_back(cidx2);
        cvec[cidx2].push_back(cidx1);
      }

      if (coef1 + negDiscount > rhs + EPS) { /* cidx1 = 1 and cidx2 = 0 */
        cvec[cidx1].push_back(cidx2 + nCols);
        cvec[cidx2 + nCols].push_back(cidx1);
      }

      if (coef2 + negDiscount > rhs + EPS) { /* cidx1 = 0 and cidx2 = 1 */
        cvec[cidx1 + nCols].push_back(cidx2);
        cvec[cidx2].push_back(cidx1 + nCols);
      }

      if (negDiscount > rhs + EPS) { /* cidx1 = 0 and cidx2 = 0 */
        cvec[cidx1 + nCols].push_back(cidx2 + nCols);
        cvec[cidx2 + nCols].push_back(cidx1 + nCols);
      }
    }
  }
}

void processClique(CGraph *cgraph, const int *idxs, const int size, char *recomputeDegree)
{
  if (size >= MIN_CLIQUE_ROW) {
    cgraph_add_clique(cgraph, idxs, size);
    *recomputeDegree = 1;
  } else {
    cgraph_add_clique_as_normal_conflicts(cgraph, idxs, size);
  }
}

int binary_search(const std::vector< std::pair< int, double > > &columns, double partialLHS, double rhs, int colStart,
  int colEnd)
{
  int mid, left = colStart, right = colEnd;

  while (left <= right) {
    mid = (left + right) / 2;
    double LHS = partialLHS - std::min(0.0, columns[mid].second) + columns[mid].second;

    if (rhs + EPS >= LHS)
      left = mid + 1;
    else
      right = mid - 1;
  }

  return right + 1;
}

int binary_search_complement(const std::vector< std::pair< int, double > > &columns, double partialLHS, double rhs,
  int colStart,
  int colEnd)
{
  int mid, left = colStart, right = colEnd;

  while (left <= right) {
    mid = (left + right) / 2;
    double LHS = partialLHS - std::min(0.0, columns[mid].second);

    if (rhs + EPS >= LHS)
      right = mid - 1;
    else
      left = mid + 1;
  }

  return left - 1;
}

void cliqueDetection(CGraph *cgraph, const std::vector< std::pair< int, double > > &columns, double sumNegCoefs, double rhs,
  char *recomputeDegree)
{
  int nElements = (int)columns.size(), cliqueStart = -1;
  double maxLHS; //maxLHS = lower bound for LHS when the two variables with highest coefficients are activated.
  int *idxs = new int[nElements];
  int cliqueSize = 0;

  maxLHS = sumNegCoefs - std::min(0.0, columns[nElements - 2].second) - std::min(0.0, columns[nElements - 1].second)
    + columns[nElements - 2].second + columns[nElements - 1].second;

  if (maxLHS <= rhs + EPS) {
    delete[] idxs;
    return; //there is no clique involving activation of variables in this constraint.
  }

  for (int i = 0; i < nElements - 1; i++) {
    double D = sumNegCoefs - std::min(0.0, columns[i].second) - std::min(0.0, columns[i + 1].second);
    double LHS = D + columns[i].second + columns[i + 1].second;

    if (LHS > rhs + EPS) {
      cliqueStart = i;
      break;
    }
  }

#ifdef DEBUG
  assert(cliqueStart >= 0 && cliqueStart < nElements - 1);
#endif
  for (int i = cliqueStart; i < nElements; i++) {
    idxs[cliqueSize++] = columns[i].first;
  }

  //process the first clique found
  processClique(cgraph, idxs, cliqueSize, recomputeDegree);

  //now we have to check the variables that are outside of the clique found.
  for (int i = cliqueStart - 1; i >= 0; i--) {
    const int idx = columns[i].first;
    const double coef = columns[i].second;
    const double partialLHS = sumNegCoefs - std::min(0.0, coef) + coef;

    maxLHS = sumNegCoefs - std::min(0.0, columns[i].second) - std::min(0.0, columns[nElements - 1].second)
      + columns[i].second + columns[nElements - 1].second;

    if (maxLHS <= rhs + EPS) {
      delete[] idxs;
      return;
    }

    int position = binary_search(columns, partialLHS, rhs, cliqueStart, nElements - 1);
#ifdef DEBUG
    assert(position >= 0 && position < nElements);
#endif
    cliqueSize = 0;
    idxs[cliqueSize++] = idx;
    for (int j = position; j < nElements; j++) {
      idxs[cliqueSize++] = columns[j].first;
    }
    processClique(cgraph, idxs, cliqueSize, recomputeDegree);
  }
  delete[] idxs;
}

void cliqueComplementDetection(CGraph *cgraph, const std::vector< std::pair< int, double > > &columns, double sumNegCoefs,
  double rhs, char *recomputeDegree)
{
  int nElements = (int)columns.size(), cliqueCompStart = -1;
  int nCols = cgraph_size(cgraph) / 2;
  double maxLHS; //maxLHS = lower bound for LHS when the two variaveis with smallest coefficients are deactivated.
  int *idxs = new int[nElements];
  int cliqueCompSize = 0;

  maxLHS = sumNegCoefs - std::min(0.0, columns[1].second) - std::min(0.0, columns[0].second);

  if (maxLHS <= rhs + EPS) {
    delete[] idxs;
    return; //there is no clique involving the complement of variables in this constraint.
  }

  for (int i = nElements - 1; i > 0; i--) {
    double D = sumNegCoefs - std::min(0.0, columns[i].second) - std::min(0.0, columns[i - 1].second);

    if (D > rhs + EPS) {
      cliqueCompStart = i;
      break;
    }
  }
#ifdef DEBUG
  assert(cliqueCompStart > 0 && cliqueCompStart < nElements);
#endif
  for (int i = 0; i < cliqueCompStart + 1; i++) {
    idxs[cliqueCompSize++] = columns[i].first + nCols; //binary complement
  }

  //process the first clique found
  processClique(cgraph, idxs, cliqueCompSize, recomputeDegree);

  //now we have to check the variables that are outside of the clique found.
  for (int i = cliqueCompStart + 1; i < nElements; i++) {
    int idx = columns[i].first;
    double coef = columns[i].second;
    double partialLHS = sumNegCoefs - std::min(0.0, coef);

    maxLHS = sumNegCoefs - std::min(0.0, columns[i].second) - std::min(0.0, columns[0].second);

    if (maxLHS <= rhs + EPS) {
      delete[] idxs;
      return;
    }

    int position = binary_search_complement(columns, partialLHS, rhs, 0, cliqueCompStart);

#ifdef DEBUG
    assert(position >= 0 && position < nElements);
#endif
    cliqueCompSize = 0;
    idxs[cliqueCompSize++] = idx + nCols;
    for (int j = 0; j <= position; j++) {
      idxs[cliqueCompSize++] = columns[j].first + nCols;
    }
    processClique(cgraph, idxs, cliqueCompSize, recomputeDegree);
  }

  delete[] idxs;
}

void mixedCliqueDetection(CGraph *cgraph, const std::vector< std::pair< int, double > > &columns, double sumNegCoefs,
  double rhs,
  std::vector< std::vector< int > > &cvec)
{
  int nElements = (int)columns.size();
  int nCols = cgraph_size(cgraph) / 2;

  //Looking for conflicts like (x = 0, y = 1)
  for (int i = 0; i < nElements - 1; i++) {
    int idx = columns[i].first;
    double coef = columns[i].second;
    double maxLHS = sumNegCoefs - std::min(0.0, coef) - std::min(0.0, columns[nElements - 1].second) + columns[nElements - 1].second;

    if (maxLHS <= rhs + EPS)
      break; //there are no conflicts here

    double partialLHS = sumNegCoefs - std::min(0.0, coef);
    int position = binary_search(columns, partialLHS, rhs, i + 1, nElements - 1);

#ifdef DEBUG
    assert(position >= i + 1 && position < nElements);
#endif

    //Conflicts between the complement of variable with index idx and
    //all the variables with indexes in the range [position, nElements-1]
    for (int j = position; j < nElements; j++) {
      cvec[idx + nCols].push_back(columns[j].first);
      cvec[columns[j].first].push_back(idx + nCols);
    }
  }

  //Looking for conflicts like (x = 1, y = 0)
  for (int i = nElements - 2; i >= 0; i--) {
    int idx = columns[i].first;
    double coef = columns[i].second;

    double maxLHS = sumNegCoefs - std::min(0.0, coef) - std::min(0.0, columns[i + 1].second) + coef;

    if (maxLHS <= rhs + EPS)
      break; //there are no conflicts here

    double partialLHS = sumNegCoefs - std::min(0.0, coef) + coef;
    int position = binary_search_complement(columns, partialLHS, rhs, i + 1, nElements - 1);

#ifdef DEBUG
    assert(position >= i + 1 && position < nElements);
#endif

    //Conflicts between the variable with index idx and
    //all the complements of variables with indexes in the range [i+1, position]
    for (int j = i + 1; j <= position; j++) {
      cvec[idx].push_back(columns[j].first + nCols);
      cvec[columns[j].first + nCols].push_back(idx);
    }
  }
}

#ifdef CGRAPH_LP
CGraph *build_cgraph(const LinearProgram *mip)
{
  if (lp_num_binary_cols(mip) < 2)
    return nullptr;

  const int nCols = lp_cols(mip);
  const int nRows = lp_rows(mip);
  const int cgraphSize = nCols * 2;
  CGraph *cgraph = cgraph_create(cgraphSize);
  int *idxs = new int[nCols];
  double *coefs = new double[nCols];
  char recomputeDegree = 0;
  std::vector< std::vector< int > > cvec(cgraphSize);

  for (int i = 0; i < nCols; i++) {
    cvec[i].reserve(128);
    /* inserting trivial conflicts: variable-complement */
    if (lp_is_binary(mip, i)) { //consider only binary variables
      cvec[i].push_back(i + nCols);
      cvec[i + nCols].push_back(i);
    }
  }

  for (int idxRow = 0; idxRow < nRows; idxRow++) {
    const int nElements = lp_row(mip, idxRow, idxs, coefs);
    const double rhs = lp_rhs(mip, idxRow);
    const char sense = lp_sense(mip, idxRow);
    std::vector< std::pair< int, double > > columns(nElements);
    int nBools = 0; // number of binary variables
    int nPos = 0; //number of positive coefficients
    double sumNegCoefs = 0.0; //sum of all negative coefficients
    double minCoef = std::numeric_limits< double >::max();
    double maxCoef = std::numeric_limits< double >::min();

    if ((nElements < 2) || (fabs(rhs) >= LARGE_CONST))
      continue;

    if (sense == 'R') { // lets not consider ranged constraints by now
      char name[256];
      printf("TODO: CHECK FOR RANGED CONSTRAINT (%s) rhs is %g\n", lp_row_name(mip, idxRow, name), rhs);
      continue;
    }

    double mult = (sense == 'G') ? -1.0 : 1.0;
    for (int i = 0; i < nElements; i++) {
      columns[i].first = idxs[i];
      columns[i].second = coefs[i] * mult;

      if (lp_is_binary(mip, columns[i].first))
        nBools++;

      if (columns[i].second <= -EPS)
        sumNegCoefs += columns[i].second;
      else
        nPos++;

      minCoef = std::min(minCoef, columns[i].second);
      maxCoef = std::max(maxCoef, columns[i].second);
    }

    if (nBools < nElements)
      continue;

    /* special case: GUB constraints */
    if (DBL_EQUAL(minCoef, maxCoef) && DBL_EQUAL(maxCoef, rhs * mult) && DBL_EQUAL(minCoef, 1.0) && ((sense == 'E') || (sense == 'L'))
      && (nElements > 3)) {
      processClique(cgraph, idxs, nElements, &recomputeDegree);
    } else {
      std::sort(columns.begin(), columns.end(), sort_columns);
      cliqueDetection(cgraph, columns, sumNegCoefs, rhs * mult, &recomputeDegree);
      cliqueComplementDetection(cgraph, columns, sumNegCoefs, rhs * mult, &recomputeDegree);
      mixedCliqueDetection(cgraph, columns, sumNegCoefs, rhs * mult, cvec);

      /*equality constraints are converted into two inequality constraints (<=).
            the first one is analyzed above and the second (multiplying the constraint by -1) is analyzed below.
            Example: x + y + z = 2 ==>  (x + y + z <= 2) and (- x - y - z <= -2)*/
      if (sense == 'E') {
        std::vector< std::pair< int, double > > newColumns(nElements);
        sumNegCoefs = 0.0;
        for (int i = 0; i < nElements; i++) {
          newColumns[i].first = columns[nElements - i - 1].first;
          newColumns[i].second = -1.0 * columns[nElements - i - 1].second;
          if (newColumns[i].second <= -EPS)
            sumNegCoefs += newColumns[i].second;
        }

        cliqueDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs, &recomputeDegree);
        cliqueComplementDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs, &recomputeDegree);
        mixedCliqueDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs, cvec);
      }
    }
  }

  for (int i = 0; i < (int)cvec.size(); i++) {
    if (!cvec[i].empty()) {
      cgraph_add_node_conflicts_no_sim(cgraph, i, &cvec[i][0], (int)cvec[i].size());
    }
  }

  if (recomputeDegree)
    cgraph_recompute_degree(cgraph);
  else
    cgraph_update_min_max_degree(cgraph);

  delete[] idxs;
  delete[] coefs;

  return cgraph;
}
#else
CGraph *build_cgraph(const CoinPackedMatrix *matrixByRow, const int numCols, const char *colType,
  const double *rhs, const char *sense)
{
  int cgraphSize = numCols * 2;
  CGraph *cgraph = cgraph_create(cgraphSize);
  int idxRow;
  char recomputeDegree = 0;
  std::vector< std::vector< int > > cvec(cgraphSize);

  for (int i = 0; i < numCols; i++) {
    cvec[i].reserve(128);
    /* inserting trivial conflicts: variable-complement */
    if (colType[i] == 1) { //consider only binary variables
      cvec[i].push_back(i + numCols);
      cvec[i + numCols].push_back(i);
    }
  }

  const int *idxs = matrixByRow->getIndices();
  const double *coefs = matrixByRow->getElements();
  const int *start = matrixByRow->getVectorStarts();
  const int *length = matrixByRow->getVectorLengths();

  for (idxRow = 0; idxRow < matrixByRow->getNumRows(); idxRow++) {
    const int nElements = length[idxRow];
    const int rowStart = start[idxRow];
    std::vector< std::pair< int, double > > columns(nElements);
    bool onlyBinaryVars = true;
    int nPos = 0; //number of positive coefficients
    double sumNegCoefs = 0.0; //sum of all negative coefficients
    double minCoef = (std::numeric_limits< double >::max() / 10.0);
    double maxCoef = -(std::numeric_limits< double >::max() / 10.0);

    if ((nElements < 2) || (fabs(rhs[idxRow]) >= LARGE_CONST))
      continue;

    if (sense[idxRow] == 'R') {
      //TODO: CHECK CONFLICTS IN RANGED CONSTRAINTS
      continue;
    }

    double mult = (sense[idxRow] == 'G') ? -1.0 : 1.0;
    for (int i = 0; i < nElements; i++) {
      columns[i].first = idxs[i + rowStart];
      columns[i].second = coefs[i + rowStart] * mult;

      if (colType[columns[i].first] != 1) {
        onlyBinaryVars = false;
        break;
      }

      if (columns[i].second <= -EPS) {
        sumNegCoefs += columns[i].second;
      } else {
        nPos++;
      }

      minCoef = std::min(minCoef, columns[i].second);
      maxCoef = std::max(maxCoef, columns[i].second);
    }

    if (!onlyBinaryVars) {
      continue;
    }

    /* special case: GUB constraints */
    if (DBL_EQUAL(minCoef, maxCoef) && DBL_EQUAL(maxCoef, rhs[idxRow] * mult) && DBL_EQUAL(minCoef, 1.0) && ((sense[idxRow] == 'E') || (sense[idxRow] == 'L'))
      && (nElements > 3)) {
      processClique(cgraph, idxs + rowStart, nElements, &recomputeDegree);
    } else {
      std::sort(columns.begin(), columns.end(), sort_columns);
      cliqueDetection(cgraph, columns, sumNegCoefs, rhs[idxRow] * mult, &recomputeDegree);
      cliqueComplementDetection(cgraph, columns, sumNegCoefs, rhs[idxRow] * mult, &recomputeDegree);
      mixedCliqueDetection(cgraph, columns, sumNegCoefs, rhs[idxRow] * mult, cvec);

      /*equality constraints are converted into two inequality constraints (<=).
            the first one is analyzed above and the second (multiplying the constraint by -1) is analyzed below.
            Example: x + y + z = 2 ==>  (x + y + z <= 2) and (- x - y - z <= -2)*/
      if (sense[idxRow] == 'E') {
        std::vector< std::pair< int, double > > newColumns(nElements);
        sumNegCoefs = 0.0;
        for (int i = 0; i < nElements; i++) {
          newColumns[i].first = columns[nElements - i - 1].first;
          newColumns[i].second = -1.0 * columns[nElements - i - 1].second;
          if (newColumns[i].second <= -EPS)
            sumNegCoefs += newColumns[i].second;
        }

        cliqueDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow], &recomputeDegree);
        cliqueComplementDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow], &recomputeDegree);
        mixedCliqueDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow], cvec);
      }
    }
  }

  for (int i = 0; i < (int)cvec.size(); i++) {
    if (!cvec[i].empty()) {
      cgraph_add_node_conflicts_no_sim(cgraph, i, &cvec[i][0], (int)cvec[i].size());
    }
  }

  if (recomputeDegree)
    cgraph_recompute_degree(cgraph);
  else
    cgraph_update_min_max_degree(cgraph);

  return cgraph;
}
#endif
