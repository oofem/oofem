/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2025   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "remeshingcrit.h"
#include "errorestimator.h"
#include "domain.h"
#include "mathfem.h"

#ifdef __MPI_PARALLEL_MODE
 #include "problemcomm.h"
#endif

namespace oofem {
RemeshingCriteria :: RemeshingCriteria(int n, ErrorEstimator *e) : FEMComponent( n, e->giveDomain() )
{
    ee = e;
#ifdef __MPI_PARALLEL_MODE
    communicator = NULL;
    commBuff = NULL;
    initCommMap = true;
#endif
}

const char *__RemeshingStrategyToString(RemeshingStrategy s)
{
    switch ( s ) {
    case NoRemeshing_RS: return "no remeshing";
    case RemeshingFromCurrentState_RS: return "remeshing from current state";
    case RemeshingFromPreviousState_RS: return "remeshing from previous state";
    }
    return "unknown";
}


void
RemeshingCriteria :: printOutputAt(FILE *file, TimeStep *tStep)
{
    int nnodes = this->domain->giveNumberOfDofManagers(), nrefine = 0, ncoarsen = 0, ndet = 0;
    double minReq = 0., maxReq = 0., minCurr = 0., maxCurr = 0.;

    RemeshingStrategy strategy = this->giveRemeshingStrategy(tStep);
    for ( int i = 1; i <= nnodes; i++ ) {
        double req = this->giveRequiredDofManDensity(i, tStep);
        if ( req < 0. ) {
            continue; // undetermined
        }
        double curr = this->giveDofManDensity(i);
        if ( ndet == 0 ) {
            minReq = maxReq = req;
            minCurr = maxCurr = curr;
        } else {
            minReq = min(minReq, req);
            maxReq = max(maxReq, req);
            minCurr = min(minCurr, curr);
            maxCurr = max(maxCurr, curr);
        }
        ndet++;
        if ( req < curr ) {
            nrefine++;
        } else if ( req > curr ) {
            ncoarsen++;
        }
    }

    fprintf(file, "  Remeshing criteria          : %s\n", __RemeshingStrategyToString(strategy) );
    fprintf(file, "  Current element size        : %e - %e\n", minCurr, maxCurr);
    fprintf(file, "  Required element size       : %e - %e\n", minReq, maxReq);
    fprintf(file, "  Nodes requiring refinement  : %d of %d\n", nrefine, nnodes);
    fprintf(file, "  Nodes allowing coarsening   : %d of %d\n", ncoarsen, nnodes);
}


RemeshingCriteria :: ~RemeshingCriteria()
{
#ifdef __MPI_PARALLEL_MODE
    delete communicator;
    delete commBuff;
#endif
}
} // end namespace oofem
