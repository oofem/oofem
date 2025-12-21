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

#ifndef dynamicrelaxationsolver_h
#define dynamicrelaxationsolver_h

#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "convergedreason.h"
#include "sparsemtrx.h"
#include "floatarray.h"
#include "linesearch.h"
#include "nrsolver.h"

#include <vector>


///@name Input fields for DynamicRelaxationSolver
//@{
#define _IFT_DynamicRelaxationSolver_Name "drsolver"
//@}

namespace oofem {
class Domain;
class EngngModel;


/**
 * Solves static equilibrium by means of explicit dynamic iterations.
 * @author Mikael Öhman
 */
class OOFEM_EXPORT DynamicRelaxationSolver : public NRSolver
{
public:
    DynamicRelaxationSolver(Domain * d, EngngModel * m);
    virtual ~DynamicRelaxationSolver() {}

    ConvergedReason solve(SparseMtrx &k, FloatArray &R, FloatArray *R0,
                    FloatArray &X, FloatArray &dX, FloatArray &F,
                    const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                    int &nite, TimeStep *) override;

    void initializeFrom(InputRecord &ir) override;

    const char *giveClassName() const override { return "DynamicRelaxationSolver"; }
    const char *giveInputRecordName() const override { return _IFT_DynamicRelaxationSolver_Name; }
};
} // end namespace oofem
#endif // dynamicrelaxationsolver_h
