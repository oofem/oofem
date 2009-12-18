/* $Header: /home/cvs/bp/oofem/sm/src/deidynamic.h,v 1.5 2003/04/06 14:08:30 bp Exp $ */
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
// Class DEIDynamic - DirectExplicitIntegrationDynamic
//

#ifndef deidynamic_h
#define deidynamic_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "structengngmodel.h"
#include "skyline.h"

namespace oofem {

class DEIDynamic : public StructuralEngngModel
{
    /*
     * This class implements Linear (- may be changed) solution of dynamic
     * problems using Direct Explicit Integration scheme - Central Difference
     * Method. For efficiency reasons it uses diagonal mass matrix
     *
     * DESCRIPTION:
     * Solution of this problem is series of loading cases, maintained as sequence of
     * time-steps. For obtaining diagonal mass matrix from possibly non-diagonal one
     * returned from Element::giveMassMatrix() a FloatMatrix::Lumped() is called
     * to obtain diagonal form.
     *
     * we start assemble governing equations at time step 0 ( 0 given by boundary and initial cond.)
     * they result in response at time step 1.
     * for time step 0 we need special start code.
     * so we obtain solution for time step 1 and next.
     * because this method is explicit, when solving equations for step t, we obtain
     * solution in step t+dt. But printing is performed for step t.
     * see diidynamic.h for difference.
     * So, when You specify initial conditions, you specify them in time step 0.
     *
     * WARNING - FloatMatrix::Lumped() works only for elements with Linear displacement filed !
     *
     * TASK:
     * Creating Numerical method for solving Ax=b
     * Interfacing Numerical method to Elements
     * Managing time  steps
     */

protected:
    FloatArray massMatrix;
    FloatArray loadVector;
    FloatArray nextDisplacementVector;
    FloatArray displacementVector, velocityVector, accelerationVector;
    double dumpingCoef, deltaT;


public:
    DEIDynamic(int i, EngngModel *_master = NULL) : StructuralEngngModel(i, _master), massMatrix(), loadVector(),
        nextDisplacementVector(), displacementVector(), velocityVector(), accelerationVector() { ndomains = 1; }
    ~DEIDynamic();
    // solving
    //void solveYourself ();
    void solveYourselfAt(TimeStep *);
    //int requiresNewLhs () {return 0;}
    /**
     * Updates nodal values
     * (calls also this->updateDofUnknownsDictionary for updating dofs unknowns dictionaries
     * if model supports changes of static system). The element internal state update is also forced using
     * updateInternalState service.
     */
    virtual void               updateYourself(TimeStep *);
    double giveUnknownComponent(EquationID, ValueModeType, TimeStep *, Domain *, Dof *);
    IRResultType initializeFrom(InputRecord *ir);
    TimeStep *giveNextStep();
    NumericalMethod *giveNumericalMethod(TimeStep *);
    // void    terminate(TimeStep*);

    /** DOF printing routine. Called by DofManagers to print Dof specific part.
     * Dof class provides component printing routines, but emodel is responsible
     * for what will be printed at DOF level.
     * @param stream output stream
     * @param iDof dof to be processed
     * @param atTime solution step
     */
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);


    // identification
    const char *giveClassName() const { return "DEIDynamic"; }
    classType giveClassID() const { return DEIDynamicClass; }
    fMode giveFormulation() { return TL; }
    virtual int        giveNumberOfFirstStep() { return 0; }
    virtual int        giveNumberOfTimeStepWhenIcApply() { return 0; }
};

} // end namespace oofem
#endif // deidynamic_h
