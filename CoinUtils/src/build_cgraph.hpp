#ifndef BUILDCGRAPH_H_INCLUDED
#define BUILDCGRAPH_H_INCLUDED

extern "C" {
	#include "cgraph.h"
}
#include "CoinPackedMatrix.hpp"

CGraph *build_cgraph(const CoinPackedMatrix *matrixByRow, const int numCols, const char *colType,
                     const double *rhs, const char *sense);

#endif
