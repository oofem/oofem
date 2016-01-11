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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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

#include "dsssolver.h"
#include "classfactory.h"

#include "dssmatrix.h"
#include "timer.h"

namespace oofem {

REGISTER_SparseLinSolver(DSSSolver, ST_DSS);

DSSSolver :: DSSSolver(Domain *d, EngngModel *m) :
    SparseLinearSystemNM(d, m) { }

DSSSolver :: ~DSSSolver() { }

NM_Status
DSSSolver :: solve(SparseMtrx &A, FloatArray &b, FloatArray &x)
{
 #ifdef TIME_REPORT
    Timer timer;
    timer.startTimer();
 #endif


    DSSMatrix *_mtrx = dynamic_cast< DSSMatrix * >(&A);
    if ( _mtrx ) {
        _mtrx->factorized();
        _mtrx->solve(b, x);
    } else {
        OOFEM_ERROR("incompatible sparse mtrx format");
    }

 #ifdef TIME_REPORT
    timer.stopTimer();
    OOFEM_LOG_INFO( "DSSSolver info: user time consumed by solution: %.2fs\n", timer.getUtime() );
 #endif

    return NM_Success;
}
} // end namespace oofem

