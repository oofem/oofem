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

#ifndef linearstability_h
#define linearstability_h

#include "structengngmodel.h"
#include "geneigvalsolvertype.h"
#include "sparsegeneigenvalsystemnm.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "nummet.h"

///@name Input fields for LinearStability
//@{
#define _IFT_LinearStability_Name "linearstability"
#define _IFT_LinearStability_nroot "nroot"
#define _IFT_LinearStability_rtolv "rtolv"
#define _IFT_LinearStability_stype "stype"
//@}

namespace oofem {
/**
 * This class implements way for examining critical load of structure.
 *
 * Solution of this problem is base on equation in the form of: @f$ K\cdot y=w (K_\sigma)y @f$.
 * Currently eigenvalue problem is solved using subspace iteration.
 * The linear static solution, determining normal forces is done in time = 0.
 *
 * Tasks:
 * - Assembling the governing equation in the form.
 * - Creating Numerical method for @f$ K\cdot y=w(K_\sigma)y @f$.
 * - Interfacing Numerical method to Elements.
 */
class LinearStability : public StructuralEngngModel
{
private:
    SparseMtrx *stiffnessMatrix;
    SparseMtrx *initialStressMatrix;
    FloatArray loadVector;
    FloatArray displacementVector;
    FloatMatrix eigVec;
    FloatArray eigVal;
    int numberOfRequiredEigenValues;
    double rtolv;
    /// Numerical method used to solve the problem.
    GenEigvalSolverType solverType;
    SparseGeneralEigenValueSystemNM *nMethod;
    /// Numerical method used to solve the static problem.
    SparseLinearSystemNM *nMethodLS;

public:
    LinearStability(int i, EngngModel *_master = NULL) : StructuralEngngModel(i, _master),
        loadVector(), displacementVector(), eigVec(), eigVal()
    {
        stiffnessMatrix = NULL;
        initialStressMatrix = NULL;
        numberOfSteps = 1;
        nMethodLS = NULL;
        nMethod = NULL;
        ndomains = 1;
    }
    virtual ~LinearStability() {
        delete  stiffnessMatrix;
        delete initialStressMatrix;
        if ( nMethodLS ) { delete nMethodLS; }
        if ( nMethod ) { delete nMethod; }
    }

    virtual void solveYourself();
    virtual void solveYourselfAt(TimeStep *tStep);

    virtual void terminate(TimeStep *tStep);
    void terminateLinStatic(TimeStep *tStep);
    int requiresNewLsh() { return 0; }
    virtual void updateYourself(TimeStep *tStep);

    // the intrinsic time of time step defines active eigen value and corresponding vector,
    // for which values can be requested using
    // giveUnknownComponent method.
    // When DisplacementVector is requested, then if time==0 linear elastic solution displacement are returned,
    // otherwise corresponding eigen vector is considered as displacement vector
    virtual double giveUnknownComponent(ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual TimeStep *giveNextStep();

    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);
    SparseLinearSystemNM *giveNumericalMethodForLinStaticProblem(TimeStep *tStep);

    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);

    // identification
    virtual const char *giveClassName() const { return "LinearStability"; }
    virtual classType giveClassID() const { return LinearStabilityClass; }
    virtual fMode giveFormulation() { return TL; }
};
} // end namespace oofem
#endif // linearstability_h
