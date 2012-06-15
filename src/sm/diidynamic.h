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

#ifndef diidynamic_h
#define diidynamic_h

#include "structengngmodel.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"

namespace oofem {
/**
 * This class implements Direct Implicit Integration of Dynamic problem
 *
 * Solution of this problem is series of loading cases, maintained as sequence of
 * time-steps. This solution is in form of linear equation system Ax=b
 * for Psi = 1 .... Newmark Method
 *     Psi >=1.37   Wilson  Method
 * dumping Matrix is assumed to be modeled as Rayleigh damping ( C = alpha*M + beta*K)
 *
 * We start to assemble governing equations at time step 0 (by boundary and initial conditions
 * we prescribe step -1 ) so solution is obtained first for step 0.
 * see deidynamic.h for difference.
 * When You specify initial conditions, you specify them in time step -1
 *
 * Tasks:
 * - Creating Numerical method for solving Ax=b.
 * - Interfacing Numerical method to Elements.
 * - Managing time steps.
 */
class DIIDynamic : public StructuralEngngModel
{
protected:
    bool initFlag;
    SparseMtrx *stiffnessMatrix;
    FloatArray loadVector, previousLoadVector, rhs;
    FloatArray displacementVector, velocityVector, accelerationVector;
    FloatArray previousDisplacementVector, previousVelocityVector, previousAccelerationVector;
    FloatArray help;
    double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
    double alpha, beta, deltaT;
    double Psi;

    LinSystSolverType solverType;
    SparseMtrxType sparseMtrxType;

    TimeDiscretizationType initialTimeDiscretization;

    /// Numerical method used to solve the problem
    SparseLinearSystemNM *nMethod;

public:
    DIIDynamic(int i, EngngModel *_master = NULL);
    virtual ~DIIDynamic();

    virtual void solveYourself();
    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep);
    virtual double giveUnknownComponent(EquationID eid, ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual TimeStep *giveNextStep();
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);

    // identification
    virtual const char *giveClassName() const { return "DIIDynamic"; }
    virtual classType giveClassID() const { return DIIDynamicClass; }
    virtual fMode giveFormulation() { return TL; }
    virtual int giveNumberOfTimeStepWhenIcApply() { return -1; }

    virtual void giveElementCharacteristicMatrix(FloatMatrix &answer, int num,
                                                 CharType type, TimeStep *tStep, Domain *domain);

    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);

    void timesMassMtrx(FloatArray &answer, FloatArray &vec, Domain *domain, TimeStep *tStep);
    void assembleLoadVector(FloatArray &_loadVector, Domain *domain, ValueModeType mode, TimeStep *tStep);
    void determineConstants(TimeStep *tStep);
    virtual int checkConsistency();
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};
} // end namespace oofem
#endif // diidynamic_h
