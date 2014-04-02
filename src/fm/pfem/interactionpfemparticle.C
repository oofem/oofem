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


#include "interactionpfemparticle.h"
#include "timestep.h"
#include "classfactory.h"

#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
#endif



namespace oofem {
REGISTER_DofManager(InteractionPFEMParticle);
/**
 * Constructor. Creates a particle with number n, belonging to aDomain.
 */
InteractionPFEMParticle :: InteractionPFEMParticle(int n, Domain *aDomain) : PFEMParticle(n, aDomain)
{ }
//from hanging node

/**
 * Gets from the source line from the data file all the data of the receiver.
 */
IRResultType
InteractionPFEMParticle :: initializeFrom(InputRecord *ir)
{
    return PFEMParticle :: initializeFrom(ir);
}

/**
 * Checks internal data consistency in node.
 */
int
InteractionPFEMParticle :: checkConsistency()
{
	return PFEMParticle :: checkConsistency();
}

void
InteractionPFEMParticle :: updateYourself(TimeStep *tStep)
{
    PFEMParticle :: updateYourself(tStep);
}

void
InteractionPFEMParticle :: printOutputAt(FILE *stream, TimeStep *stepN)
{
    PFEMParticle :: printOutputAt(stream, stepN);
}

#ifdef __OOFEG
void InteractionPFEMParticle :: drawScalar(oofegGraphicContext &gc)
{
	PFEMParticle :: drawScalar(gc);
}
#endif
} // end namespace oofem
