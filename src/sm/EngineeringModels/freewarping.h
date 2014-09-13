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

#ifndef freewarping_h
#define freewarping_h

#include "../sm/EngineeringModels/structengngmodel.h"
#include "sparselinsystemnm.h"
#include "sparsemtrxtype.h"

#define _IFT_FreeWarping_Name "freewarping"

namespace oofem {
class SparseMtrx;

/**
 * This class implements the free warping engineering problem 
 * (evaluation of the warping function and torsional stiffness for a given cross section).
 * Only one time step is needed and the analysis is linear.
 * The material should be linear elastic, or any other material which has a well-defined
 * shear modulus of elasticity. 
 *
 * This problem leads to a linear equation system Ax=b
 *
 * Tasks:
 * - Creating Numerical method for solving @f$ K\cdot x=b @f$.
 * - Interfacing Numerical method to Elements.
 */
class FreeWarping : public StructuralEngngModel
{
protected:
    SparseMtrx *stiffnessMatrix;
    FloatArray loadVector;
    FloatArray displacementVector;

    LinSystSolverType solverType;
    SparseMtrxType sparseMtrxType;
    /// Numerical method used to solve the problem.
    SparseLinearSystemNM *nMethod;

    int initFlag; // needed?

public:
    FreeWarping(int i, EngngModel * _master = NULL);
    virtual ~FreeWarping();

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

    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep);

    // identification
    virtual const char *giveInputRecordName() const { return _IFT_FreeWarping_Name; }
    virtual const char *giveClassName() const { return "FreeWarping"; }
    virtual fMode giveFormulation() { return TL; }

#ifdef __PARALLEL_MODE
    virtual int estimateMaxPackSize(IntArray &commMap, CommunicationBuffer &buff, int packUnpackType);
#endif
};
} // end namespace oofem
#endif // freewarping_h
