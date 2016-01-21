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

#ifndef deidynamic_h
#define deidynamic_h

#include "../sm/EngineeringModels/structengngmodel.h"

///@name Input fields for DEIDynamic
//@{
#define _IFT_DEIDynamic_Name "deidynamic"
#define _IFT_DEIDynamic_dumpcoef "dumpcoef"
#define _IFT_DEIDynamic_deltat "deltat"
//@}

namespace oofem {
/**
 * This class implements Linear (- may be changed) solution of dynamic
 * problems using Direct Explicit Integration scheme - Central Difference
 * Method. For efficiency reasons it uses diagonal mass matrix
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
 * Tasks:
 * - Creating Numerical method for solving Ax=b
 * - Interfacing Numerical method to Elements
 * - Managing time steps
 */
class DEIDynamic : public StructuralEngngModel
{
protected:
    FloatArray massMatrix;
    FloatArray loadVector;
    FloatArray nextDisplacementVector;
    FloatArray displacementVector, velocityVector, accelerationVector;
    double dumpingCoef, deltaT;

public:
    DEIDynamic(int i, EngngModel * _master = NULL) : StructuralEngngModel(i, _master), massMatrix(), loadVector(),
        nextDisplacementVector(), displacementVector(), velocityVector(), accelerationVector() {
        ndomains = 1;
    }
    virtual ~DEIDynamic();

    virtual void solveYourselfAt(TimeStep *tStep);

    virtual double giveUnknownComponent(ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual TimeStep *giveNextStep();
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);

    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep);

    // identification
    virtual const char *giveClassName() const { return "DEIDynamic"; }
    virtual const char *giveInputRecordName() const { return _IFT_DEIDynamic_Name; }
    virtual fMode giveFormulation() { return TL; }
    virtual int giveNumberOfFirstStep(bool force = false) { return 0; }
};
} // end namespace oofem
#endif // deidynamic_h
