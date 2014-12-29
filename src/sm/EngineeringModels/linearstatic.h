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

#ifndef linearstatic_h
#define linearstatic_h

#include "../sm/EngineeringModels/structengngmodel.h"
#include "sparselinsystemnm.h"
#include "sparsemtrxtype.h"

#define _IFT_LinearStatic_Name "linearstatic"

namespace oofem {
class SparseMtrx;

/**
 * This class implements linear static engineering problem.
 * Multiple loading works only if linear elastic material (such as isoLE)  is used.
 * (Other non-linear materials keep load history, so such multiple loading
 * will cause that next step will be assumed as new load increment,
 * not the total new load). Because they always compute real stresses according
 * to reached strain state, they are not able to respond to linear analysis.
 *
 * Solution of this problem is series of loading cases, maintained as sequence of
 * time-steps. This solution is in form of linear equation system Ax=b
 * Tasks:
 * - Creating Numerical method for solving @f$ K\cdot x=b @f$.
 * - Interfacing Numerical method to Elements.
 * - Managing time steps.
 */
class LinearStatic : public StructuralEngngModel
{
protected:
    std :: unique_ptr< SparseMtrx > stiffnessMatrix;
    FloatArray loadVector;
    FloatArray displacementVector;

    LinSystSolverType solverType;
    SparseMtrxType sparseMtrxType;
    /// Numerical method used to solve the problem.
    std :: unique_ptr< SparseLinearSystemNM > nMethod;

    int initFlag;

public:
    LinearStatic(int i, EngngModel * _master = NULL);
    virtual ~LinearStatic();

    virtual void solveYourself();
    virtual void solveYourselfAt(TimeStep *tStep);

    virtual double giveUnknownComponent(ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual void updateDomainLinks();

    virtual TimeStep *giveNextStep();
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);
    virtual void terminate(TimeStep *tStep);

    virtual IRResultType initializeFrom(InputRecord *ir);

    // identification
    virtual const char *giveInputRecordName() const { return _IFT_LinearStatic_Name; }
    virtual const char *giveClassName() const { return "LinearStatic"; }
    virtual fMode giveFormulation() { return TL; }

    virtual int estimateMaxPackSize(IntArray &commMap, DataStream &buff, int packUnpackType);
};
} // end namespace oofem
#endif // linearstatic_h
