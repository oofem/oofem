/* $Header: /home/cvs/bp/oofem/sm/src/linearstability.h,v 1.7 2003/04/06 14:08:30 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//
// Class LinearStability
//

#ifndef linearstability_h
#define linearstability_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "structengngmodel.h"
#include "sparsegeneigenvalsystemnm.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "nummet.h"
class LinearStability : public StructuralEngngModel
{
    /*
     * This class implements way for examining critical load of structure.
     * DESCRIPTION:
     * Solution of this problem is base on equation in the form of: Ky=w(K_sigma)y
     * Currently eigen value problem is solved using subspace iteration.
     * The linear static solution, determining normal forces is done in time = 0.
     * TASK:
     * Assembling the governing equation in the form  Ky=w(K_sigma)y
     * Creating Numerical method for Ky=w(K_sigma)y
     * Interfacing Numerical method to Elements
     */


private:
    SparseMtrx *stiffnessMatrix;
    SparseMtrx *initialStressMatrix;
    FloatArray loadVector;
    FloatArray displacementVector;
    FloatMatrix eigVec;
    FloatArray eigVal;
    int numberOfRequiredEigenValues;
    double rtolv;           // precision
    /// Numerical method used to solve the problem
    GenEigvalSolverType solverType;
    SparseGeneralEigenValueSystemNM *nMethod;
    /// Numerical method used to solve the static problem
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
    ~LinearStability() {
        delete  stiffnessMatrix;
        delete initialStressMatrix;
        if ( nMethodLS ) { delete nMethodLS; }

        if ( nMethod ) { delete nMethod; } }
    // solving
    void solveYourself();
    void solveYourselfAt(TimeStep *);

    void terminate(TimeStep *);
    void terminateLinStatic(TimeStep *);
    int requiresNewLsh() { return 0; }
    virtual void               updateYourself(TimeStep *);

    // the intrinsic time of time step defines active eigen value and corresponding vector,
    // for which values can be requested using
    // giveUnknownComponent method.
    // When DisplacementVector is requested, then if time==0 linear elastic solution displacement are returned,
    // otherwise corresponding eigen vector is considered as displacement vector
    double giveUnknownComponent(EquationID, ValueModeType, TimeStep *, Domain *, Dof *);
    double giveUnknownComponent(UnknownType, ValueModeType, TimeStep *, Domain *, Dof *);
    IRResultType initializeFrom(InputRecord *ir);
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    TimeStep *giveNextStep();

    NumericalMethod *giveNumericalMethod(TimeStep *);
    SparseLinearSystemNM *giveNumericalMethodForLinStaticProblem(TimeStep *);

    /** DOF printing routine. Called by DofManagers to print Dof specific part.
     * Dof class provides component printing routines, but emodel is responsible
     * for what will be printed at DOF level.
     * @param stream output stream
     * @param iDof dof to be processed
     * @param atTime solution step
     */
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);


    // identification
    const char *giveClassName() const { return "LinearStability"; }
    classType giveClassID()      const { return LinearStabilityClass; }
    fMode giveFormulation() { return TL; }
};

#endif // linearstability_h
