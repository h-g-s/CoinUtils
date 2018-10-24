#ifndef BUILDCGRAPH_H_INCLUDED
#define BUILDCGRAPH_H_INCLUDED

#include "cgraph.h"

#ifdef CGRAPH_LP
	#include "lp.h"
	CGraph *build_cgraph(const LinearProgram *mip);
#else
	#include "CoinPackedMatrix.hpp"
    CGraph *build_cgraph(const CoinPackedMatrix *matrixByRow, const int numCols, const char *colType,
                         const double *rhs, const char *sense);
#endif

#endif
