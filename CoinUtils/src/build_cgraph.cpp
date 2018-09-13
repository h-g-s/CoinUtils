#include <cstring>
#include <cmath>
#include <limits>
#include <cfloat>
#include <algorithm>
#include <climits>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include "build_cgraph.hpp"

extern "C"
{
    #include "memory.h"
    #include "macros.h"
}

using namespace std;

#define EPS 1e-8
#define MIN_CLIQUE_ROW 256 /* mininum size for a row to be considered a clique row */

const double LARGE_CONST = std::min( DBL_MAX/10.0, 1e20 );
char recomputeDegree = 0;

struct sort_columns
{
    bool operator()(const std::pair<int, double> &left, const std::pair<int,double> &right)
    {
        if ( fabs(left.second - right.second) > EPS )
            return ( left.second < right.second );

        return left.first < right.first;
    }
};

//standard probing technique
void pairwiseAnalysis(CGraph* cgraph, const int numCols, const vector<pair<int, double> >& rowColumns, const double sumNegCoefs, const double rhs);

void processClique(const int *idxs, const int rowStart, const int rowLength, CGraph *cgraph);

/* Returns the first position of columns which the lower bound for LHS (considering activation of variables) is greater than rhs */
/* colStart=initial position for search in columns, colEnd=last position for search in columns */
/* partialLHS = LHS calculated with only one variable */
int binary_search(const vector< pair<int, double> >& columns, double partialLHS, double rhs, int colStart, int colEnd);

/* Returns the first position of columns which the lower bound for LHS (considering deactivation of variables) is greater than rhs */
/* colStart=initial position for search in columns, colEnd=last position for search in columns */
/* partialLHS = LHS calculated with only one variable */
int binary_search_complement(const vector< pair<int, double> >& columns, double partialLHS, double rhs, int colStart, int colEnd);

/* Searches for cliques involving the activation of variables in this constraint. */
void cliqueDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, double sumNegCoefs, double rhs);

/* Searches for cliques involving the complement of variables in this constraint. */
void cliqueComplementDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, double sumNegCoefs, double rhs);

/* Searches for cliques involving variables and complements of variables in this constraint. */
void mixedCliqueDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, double sumNegCoefs, double rhs);

void pairwiseAnalysis(CGraph* cgraph, const int numCols, const vector<pair<int, double> >& rowColumns, const double sumNegCoefs, const double rhs)
{
    int nElements = (int)rowColumns.size();

    for(int j1 = 0; j1 < nElements; j1++) {
        const int cidx1 = rowColumns[j1].first;
        const double coef1 = rowColumns[j1].second;

        for(int j2 = j1+1; j2 < nElements; j2++) {
            const int cidx2 = rowColumns[j2].first;
            const double coef2 = rowColumns[j2].second;
            const double negDiscount = sumNegCoefs - min(0.0, coef1) - min(0.0, coef2);

            if(coef1 + coef2 + negDiscount > rhs + EPS) /* cidx1 = 1 and cidx2 = 1 */
                cgraph_add_node_conflict(cgraph, cidx1, cidx2);

            if(coef1 + negDiscount > rhs + EPS) /* cidx1 = 1 and cidx2 = 0 */
                cgraph_add_node_conflict(cgraph, cidx1, cidx2+numCols);

            if(coef2 + negDiscount > rhs + EPS) /* cidx1 = 0 and cidx2 = 1 */
                cgraph_add_node_conflict(cgraph, cidx1+numCols, cidx2);

            if(negDiscount > rhs + EPS) /* cidx1 = 0 and cidx2 = 0 */
                cgraph_add_node_conflict(cgraph, cidx1+numCols, cidx2+numCols);
        }
    }
}

void processClique(const int *idxs, const int rowStart, const int rowLength, CGraph *cgraph)
{
    if(rowLength >= MIN_CLIQUE_ROW) {
        cgraph_add_clique( cgraph, idxs + rowStart, rowLength );
        recomputeDegree = 1;
    } else {
        int i1, i2;
        for(i1 = 0; i1 < rowLength-1; i1++)
            for(i2 = i1 + 1; i2 < rowLength; i2++)
                cgraph_add_node_conflict(cgraph, idxs[i1+rowStart], idxs[i2+rowStart]);
    }
}

int binary_search(const vector< pair<int, double> >& columns, double partialLHS, double rhs, int colStart, int colEnd)
{
    int mid;
    while(colStart <= colEnd)
    {
        mid = (colStart + colEnd) / 2;
        double LHS = partialLHS - min(0.0, columns[mid].second) + columns[mid].second;

        if(rhs + EPS >= LHS)
            colStart = mid + 1;
        else
            colEnd = mid - 1;
    }

    return colEnd + 1;
}

int binary_search_complement(const vector< pair<int, double> >& columns, double partialLHS, double rhs, int colStart, int colEnd)
{
    int mid;
    while(colStart <= colEnd)
    {
        mid = (colStart + colEnd) / 2;
        double LHS = partialLHS - min(0.0, columns[mid].second);

        if(rhs + EPS >= LHS)
            colEnd = mid - 1;
        else
            colStart = mid + 1;
    }

    return colStart - 1;
}

void cliqueDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, double sumNegCoefs, double rhs)
{
    int nElements = (int)columns.size(), cliqueStart = -1;
    double maxLHS; //maxLHS = lower bound for LHS when the two variables with highest coefficients are activated.
    int cliqueSize = 0;

    maxLHS = sumNegCoefs - min(0.0, columns[nElements-2].second) - min(0.0, columns[nElements-1].second)
             + columns[nElements-2].second + columns[nElements-1].second;

    if(maxLHS <= rhs + EPS) return; //there is no clique involving activation of variables in this constraint.

    for(int i = 0; i < nElements - 1; i++)
    {
        double D = sumNegCoefs - min(0.0, columns[i].second) - min(0.0, columns[i+1].second);
        double LHS = D + columns[i].second + columns[i+1].second;

        if(LHS > rhs + EPS)
        {
            cliqueStart = i;
            break;
        }
    }

    assert(cliqueStart >= 0 && cliqueStart < nElements - 1);
    int n = nElements - cliqueStart, idxs[n];
    for(int i = cliqueStart, j = 0; i < nElements; i++)
    {
        idxs[j++] = columns[i].first;
        cliqueSize++;
    }
    //process the first clique found
    processClique(idxs, 0, cliqueSize, cgraph);

    //now we have to check the variables that are outside of the clique found.
    for(int i = cliqueStart - 1; i >= 0; i--)
    {
        int idx = columns[i].first;
        double coef = columns[i].second;
        double partialLHS = sumNegCoefs - min(0.0, coef) + coef;

        maxLHS = sumNegCoefs - min(0.0, columns[i].second) - min(0.0, columns[nElements-1].second)
                 + columns[i].second + columns[nElements-1].second;

        if(maxLHS <= rhs + EPS) return;

        int position = binary_search(columns, partialLHS, rhs, cliqueStart, nElements - 1);

        assert(position >= 0 && position < nElements);

        int n = nElements - position + 1, idxs[n];
        cliqueSize = 1;
        idxs[0] = idx;
        for(int i = position, j = 1; i < nElements; i++)
        {
            idxs[j++] = columns[i].first;
            cliqueSize++;
        }
        processClique(idxs, 0, cliqueSize, cgraph);
    }
}

void cliqueComplementDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, double sumNegCoefs, double rhs)
{
    int nElements = (int)columns.size(), cliqueCompStart = -1;
    int nCols = cgraph_size(cgraph) / 2;
    double maxLHS; //maxLHS = lower bound for LHS when the two variaveis with smallest coefficients are deactivated.
    int cliqueCompSize = 0;

    maxLHS = sumNegCoefs - min(0.0, columns[1].second) - min(0.0, columns[0].second);

    if(maxLHS <= rhs + EPS) return; //there is no clique involving the complement of variables in this constraint.

    for(int i = nElements - 1; i > 0; i--)
    {
        double D = sumNegCoefs - min(0.0, columns[i].second) - min(0.0, columns[i-1].second);

        if(D > rhs + EPS)
        {
            cliqueCompStart = i;
            break;
        }
    }

    assert(cliqueCompStart > 0 && cliqueCompStart < nElements);
    int n = cliqueCompStart + 1, idxs[n];
    for(int i = 0; i < n; i++)
    {
        idxs[i] = columns[i].first + nCols; //binary complement
        cliqueCompSize++;
    }

    //process the first clique found
    processClique(idxs, 0, cliqueCompSize, cgraph);

    //now we have to check the variables that are outside of the clique found.
    for(int i = cliqueCompStart + 1; i < nElements; i++)
    {
        int idx = columns[i].first;
        double coef = columns[i].second;
        double partialLHS = sumNegCoefs - min(0.0, coef);

        maxLHS = sumNegCoefs - min(0.0, columns[i].second) - min(0.0, columns[0].second);

        if(maxLHS <= rhs + EPS) return;

        int position = binary_search_complement(columns, partialLHS, rhs, 0, cliqueCompStart);

        assert(position >=0 && position < nElements);

        int n = position + 2, idxs[n];
        cliqueCompSize = 1;
        idxs[0] = idx + nCols;
        for(int i = 0, j = 1; i <= position; i++)
        {
            idxs[j++] = columns[i].first + nCols;
            cliqueCompSize++;
        }
        processClique(idxs, 0, cliqueCompSize, cgraph);
    }
}

void mixedCliqueDetection(CGraph* cgraph, const vector<pair<int, double> >& columns, double sumNegCoefs, double rhs)
{
    int nElements = (int)columns.size();
    int nCols = cgraph_size(cgraph) / 2;

    //Looking for conflicts like (x = 0, y = 1)
    for(int i = 0; i < nElements - 1; i++)
    {
        int idx = columns[i].first;
        double coef = columns[i].second;
        double maxLHS = sumNegCoefs - min(0.0, coef) - min(0.0, columns[nElements-1].second) + columns[nElements-1].second;

        if(maxLHS <= rhs + EPS) break; //there are no conflicts here

        double partialLHS = sumNegCoefs - min(0.0, coef);
        int position = binary_search(columns, partialLHS, rhs, i+1, nElements - 1);

        assert(position >= i+1 && position < nElements);

        //Conflicts between the complement of variable with index idx and
        //all the variables with indexes in the range [position, nElements-1]
        for(int j = position; j < nElements; j++)
            cgraph_add_node_conflict(cgraph, idx + nCols, columns[j].first);
    }

    //Looking for conflicts like (x = 1, y = 0)
    for(int i = nElements - 2; i >= 0; i--)
    {
        int idx = columns[i].first;
        double coef = columns[i].second;

        double maxLHS = sumNegCoefs - min(0.0, coef) - min(0.0, columns[i+1].second) + coef;

        if(maxLHS <= rhs + EPS) break; //there are no conflicts here

        double partialLHS = sumNegCoefs - min(0.0, coef) + coef;
        int position = binary_search_complement(columns, partialLHS, rhs, i+1, nElements - 1);

        assert(position >= i+1 && position < nElements);

        //Conflicts between the variable with index idx and
        //all the complements of variables with indexes in the range [i+1, position]
        for(int j = i+1; j <= position; j++)
            cgraph_add_node_conflict(cgraph, idx, columns[j].first + nCols);
    }
}

CGraph *build_cgraph(const CoinPackedMatrix *matrixByRow, const int numCols, const char *colType,
                     const double *rhs, const char *sense) {
    recomputeDegree = 0;

    int cgraphSize = numCols * 2;
    CGraph *cgraph = cgraph_create( cgraphSize );
    int idxRow;

    const int* idxs = matrixByRow->getIndices();
    const double *coefs = matrixByRow->getElements();
    const int *start = matrixByRow->getVectorStarts();
    const int *length = matrixByRow->getVectorLengths();

    for(int i = 0; i < numCols; i++)
        /* inserting trivial conflicts: variable-complement */
        if(colType[i] == 1) //consider only binary variables
            cgraph_add_node_conflict(cgraph, i, i + numCols);

    for(idxRow = 0; idxRow < matrixByRow->getNumRows(); idxRow++)
    {
        const int nElements = length[idxRow];
        const int rowStart = start[idxRow];
        vector< pair<int, double> > columns(nElements);
        bool onlyBinaryVars = true;
        int nPos = 0; //number of positive coefficients
        double sumNegCoefs = 0.0; //sum of all negative coefficients
        double minCoef = numeric_limits<double>::max();
        double maxCoef = numeric_limits<double>::min();

        if ( (nElements<2) || (fabs(rhs[idxRow])>=LARGE_CONST) )
            continue;

        if ( sense[idxRow] == 'R' )  // lets not consider ranged constraints by now
        {
            //TODO: CHECK CONFLICTS IN RANGED CONSTRAINTS
            continue;
        }

        double mult = (sense[idxRow] == 'G') ? -1.0 : 1.0;
        for(int i = 0; i < nElements; i++)
        {
            columns[i].first = idxs[i+rowStart];
            columns[i].second = coefs[i+rowStart] * mult;

            if(colType[columns[i].first] != 1) {
                onlyBinaryVars = false;
                break;
            }

            if(columns[i].second <= -EPS)
                sumNegCoefs += columns[i].second;
            else nPos++;

            minCoef = min(minCoef, columns[i].second);
            maxCoef = max(maxCoef, columns[i].second);
        }

        if(!onlyBinaryVars)
            continue;

        /* special case: GUB constraints */
        if ( DBL_EQUAL( minCoef, maxCoef ) &&  DBL_EQUAL( maxCoef, rhs[idxRow] * mult ) &&
             DBL_EQUAL(minCoef, 1.0) && ((sense[idxRow]=='E') || (sense[idxRow]=='L'))
             && (nElements > 3) )
        {
            processClique(idxs, rowStart, nElements, cgraph);
        }

        else
        {
            sort(columns.begin(), columns.end(), sort_columns());
            cliqueDetection(cgraph, columns, sumNegCoefs, rhs[idxRow] * mult);
            cliqueComplementDetection(cgraph, columns, sumNegCoefs, rhs[idxRow] * mult);
            mixedCliqueDetection(cgraph, columns, sumNegCoefs, rhs[idxRow] * mult);

            /*equality constraints are converted into two inequality constraints (<=).
            the first one is analyzed above and the second (multiplying the constraint by -1) is analyzed below.
            Example: x + y + z = 2 ==>  (x + y + z <= 2) and (- x - y - z <= -2)*/
            if(sense[idxRow] == 'E')
            {
                vector<pair<int, double> > newColumns(nElements);
                sumNegCoefs = 0.0;
                for(int i = 0; i < nElements; i++)
                {
                    newColumns[i].first = columns[nElements-i-1].first;
                    newColumns[i].second = -1.0 * columns[nElements-i-1].second;
                    if(newColumns[i].second <= -EPS)
                        sumNegCoefs += newColumns[i].second;
                }

                cliqueDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow]);
                cliqueComplementDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow]);
                mixedCliqueDetection(cgraph, newColumns, sumNegCoefs, -1.0 * rhs[idxRow]);
            }
        }
    }

    if(recomputeDegree)
        cgraph_recompute_degree( cgraph );
    else
        cgraph_update_min_max_degree( cgraph );

    return cgraph;
}