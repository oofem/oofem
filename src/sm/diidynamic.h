/* $Header: /home/cvs/bp/oofem/sm/src/diidynamic.h,v 1.7 2003/04/06 14:08:30 bp Exp $ */
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
// Class DIIDynamic
//

#ifndef diidynamic_h
#define diidynamic_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "structengngmodel.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"

class DIIDynamic : public StructuralEngngModel
{
    /*
     * This class implements Direct Implicit Integration of Dynamic problem
     * DESCRIPTION:
     * Solution of this problem is series of loading cases, maintained as sequence of
     * time-steps. This solution is in form of linear equation system Ax=b
     * for Psi = 1 .... Newmark Method
     *     PSI >=1.37   Wilson  Method
     * dumping Matrix is assumed to be modelled as Raileigh damping ( C = alpha*M + beta*K)
     *
     * we start to assemble governing equations at time step 0 (by boundary and initial conditions
     * we prescribe step -1 ) so solution is obtained first for step 0.
     * see deidynamic.h for difference.
     * When You specify initial conditions, you specify them in time step -1
     *
     * TASK:
     * Creating Numerical method for solving Ax=b
     * Interfacing Numerical method to Elements
     * Managing time  steps
     */

protected:
    SparseMtrx *stiffnessMatrix, *massMatrix;
    FloatArray loadVector, previousLoadVector, rhs;
    FloatArray displacementVector, velocityVector, accelerationVector;
    double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
    double alpha, beta, deltaT;
    double Psi;
    /// Numerical method used to solve the problem
    SparseLinearSystemNM *nMethod;


public:
    DIIDynamic(int i, EngngModel *_master = NULL) : StructuralEngngModel(i, _master),  loadVector(), previousLoadVector(),
        rhs(), displacementVector(), velocityVector(), accelerationVector()
    { stiffnessMatrix = NULL;
      massMatrix = NULL;
      ndomains = 1;
      nMethod = NULL; }
    ~DIIDynamic() { delete  stiffnessMatrix;
                    delete massMatrix;
                    if ( nMethod ) { delete nMethod; } }
    // solving
    //void solveYourself ();
    void solveYourselfAt(TimeStep *);
    //int requiresNewLhs () {return 0;}
    virtual void               updateYourself(TimeStep *);
    double giveUnknownComponent(EquationID, ValueModeType, TimeStep *, Domain *, Dof *);
    IRResultType initializeFrom(InputRecord *ir);
    TimeStep *giveNextStep();
    NumericalMethod *giveNumericalMethod(TimeStep *);

    // identification
    const char *giveClassName() const { return "DIIDynamic"; }
    classType giveClassID() const { return DIIDynamicClass; }
    fMode giveFormulation() { return TL; }
    virtual int        giveNumberOfFirstStep() { return 0; }
    virtual int        giveNumberOfTimeStepWhenIcApply() { return -1; }

    virtual void giveElementCharacteristicMatrix(FloatMatrix &answer, int num,
                                                 CharType type, TimeStep *tStep, Domain *domain);

    /** DOF printing routine. Called by DofManagers to print Dof specific part.
     * Dof class provides component printing routines, but emodel is responsible
     * for what will be printed at DOF level.
     * @param stream output stream
     * @param iDof dof to be processed
     * @param atTime solution step
     */
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);

protected:
};

#endif // diidynamic_h
