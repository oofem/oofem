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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "dsssolver.h"

#ifdef __DSS_MODULE

 #include "dss.h"

namespace oofem {
DSSSolver :: DSSSolver(int i, Domain *d, EngngModel *m) :
    SparseLinearSystemNM(i, d, m) { }

DSSSolver ::  ~DSSSolver() { }

NM_Status
DSSSolver :: solve(SparseMtrx *A, FloatArray *b, FloatArray *x)
{
 #ifdef TIME_REPORT
    //clock_t tstart = clock();
    oofem_timeval tstart;
    :: getUtime(tstart);
 #endif


    DSSMatrix *_mtrx = dynamic_cast< DSSMatrix * >(A);
    if ( _mtrx ) {
        _mtrx->factorized();
        _mtrx->solve(b, x);
    } else {
        OOFEM_ERROR("DSSSolver::solve : incompatible sparse mtrx format");
    }

 #ifdef TIME_REPORT
    oofem_timeval ut;
    :: getRelativeUtime(ut, tstart);
    OOFEM_LOG_INFO( "DSSSolver info: user time consumed by solution: %.2fs\n", ( double ) ( ut.tv_sec + ut.tv_usec / ( double ) OOFEM_USEC_LIM ) );
 #endif

    return NM_Success;
}

IRResultType
DSSSolver :: initializeFrom(InputRecord *ir)
//
//
//
{
    return IRRT_OK;
}
} // end namespace oofem

#else // __DSS_MODULE

namespace oofem {
DSSSolver :: DSSSolver(int i, Domain *d, EngngModel *m) : SparseLinearSystemNM(i, d, m)
{
    _error("DSSSolver: can't create, DSS support not compiled");
}

DSSSolver :: ~DSSSolver() { }

IRResultType
DSSSolver :: initializeFrom(InputRecord *ir) { return IRRT_OK; }

NM_Status
DSSSolver :: solve(SparseMtrx *A, FloatArray *b, FloatArray *x) { return NM_NoSuccess; }
} // end namespace oofem
#endif
