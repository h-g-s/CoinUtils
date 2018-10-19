#ifndef BUILDCGRAPH_H_INCLUDED
#define BUILDCGRAPH_H_INCLUDED

#include "cgraph.h"

#ifdef COIN_CGRAPH
#include "CoinPackedMatrix.hpp"
    CGraph *build_cgraph(const CoinPackedMatrix *matrixByRow, const int numCols, const char *colType,
                         const double *rhs, const char *sense);
#else
#include "lp.h"
    CGraph *build_cgraph(const LinearProgram *mip);
#endif

#endif
