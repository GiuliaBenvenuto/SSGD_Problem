#ifndef EXTENDED_SOLVER
#define EXTENDED_SOLVER
#include "shortest_path.h"

geodesic_solver extended_solver(const DrawableTrimesh<> &m,
                                const dual_geodesic_solver &solver, int k);

#endif