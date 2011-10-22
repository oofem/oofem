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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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
// Class NlDEIDynamic - DirectExplicitIntegrationDynamic
//

#ifndef nldeidynamic_h
#define nldeidynamic_h

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif
#include "structengngmodel.h"
#include "skyline.h"

#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
#define LOCAL_ZERO_MASS_REPLACEMENT 1

/**
 * This class implements NonLinear (- may be changed) solution of dynamic
 * problems using Direct Explicit Integration scheme - Central Difference
 * Method. For efficiency reasons it uses diagonal mass matrix. It is formulated
 * in increments of displacements rather than in total variables.
 *
 * Description:
 * Solution of this problem is series of loading cases, maintained as sequence of
 * time-steps. For obtaining diagonal mass matrix from possibly non-diagonal one
 * returned from Element::giveMassMatrix() a FloatMatrix::Lumped() is called
 * to obtain diagonal form.
 *
 * We start assemble governing equations at time step 0 ( 0 given by boundary and initial cond.)
 * they result in response at time step 1.
 * for time step 0 we need special start code.
 * so we obtain solution for time step 1 and next.
 * because this method is explicit, when solving equations for step t, we obtain
 * solution in step t+dt. But printing is performed for step t.
 * see diidynamic.h for difference.
 * So, when You specify initial conditions, you specify them in time step 0.
 *
 * WARNING - FloatMatrix::Lumped() works only for elements with Linear displacement filed !
 */
class NlDEIDynamic : public StructuralEngngModel
{
protected:
    FloatArray massMatrix;
    FloatArray loadVector;
    FloatArray previousIncrementOfDisplacementVector;
    FloatArray incrementOfDisplacementVector, displacementVector,
               velocityVector, accelerationVector;
    double dumpingCoef, deltaT;
    /// Flag indicating the need for initialization
    int initFlag;

    // dynamic relaxation specific vars
    /// Flag indicating whether dynamic relaxation takes place
    int drFlag;
    /// Reference load vector
    FloatArray loadRefVector;
    /// Parameter determining rate of the loading process.
    double c;
    /// End of time interval.
    double Tau;
    /// Estimate of loadRefVector^T*displacementVector(Tau).
    double pyEstimate;
    /// Product of p^tM^(-1)p; where p is reference load vector.
    double pMp;

public:
    NlDEIDynamic(int i, EngngModel *_master = NULL) : StructuralEngngModel(i, _master), massMatrix(), loadVector(),
        previousIncrementOfDisplacementVector(), incrementOfDisplacementVector(),
        displacementVector(), velocityVector(), accelerationVector() {
        ndomains = 1;
        initFlag = 1;
    }
    ~NlDEIDynamic();

    void solveYourselfAt(TimeStep *tStep);

    virtual void updateYourself(TimeStep *tStep);
    double giveUnknownComponent(EquationID eid, ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);
    IRResultType initializeFrom(InputRecord *ir);
    TimeStep *giveNextStep();
    NumericalMethod *giveNumericalMethod(TimeStep *tStep);
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    void terminate(TimeStep *tStep);
    void giveInternalForces(FloatArray &answer, TimeStep *tStep);

    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);


    // identification
    const char *giveClassName() const { return "NlDEIDynamic"; }
    classType giveClassID() const { return NlDEIDynamicClass; }
    fMode giveFormulation() { return nonLinFormulation; }

    virtual int giveNumberOfFirstStep() { return 0; }
    virtual int giveNumberOfTimeStepWhenIcApply() { return 0; }
};
} // end namespace oofem
#endif // nldeidynamic_h
