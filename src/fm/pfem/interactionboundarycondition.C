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

#include "interactionboundarycondition.h"
#include "timestep.h"
#include "function.h"
#include "verbose.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
//#include "dofmanager.h"
#include "interactionpfemparticle.h"
#include "dof.h"

namespace oofem {
REGISTER_BoundaryCondition(InteractionBoundaryCondition);

double InteractionBoundaryCondition :: give(Dof *dof, ValueModeType mode, TimeStep *stepN)
// Returns the value at stepN of the prescribed value of the kinematic
// unknown 'u'. Returns 0 if 'u' has no prescribed value.
{
    double value = 0.0;

    InteractionPFEMParticle *interactionParticle = dynamic_cast< InteractionPFEMParticle * >( dof->giveDofManager() );
    if ( interactionParticle ) {
        FloatArray velocities;
        interactionParticle->giveCoupledVelocities(velocities, stepN);
        if ( dof->giveDofID() == V_u ) {
            value = velocities.at(1);
        } else if ( dof->giveDofID() == V_v ) {
            value = velocities.at(2);
        }
    }

    return value;
}
} // end namespace oofem
