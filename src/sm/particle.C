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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "particle.h"
#include "classfactory.h"
#include "error.h"
#include "floatmatrix.h"
#include "domain.h"
#include "parametermanager.h"
#include "paramkey.h"

namespace oofem {
REGISTER_DofManager(Particle);

ParamKey Particle::IPK_Particle_rad("rad");

Particle :: Particle(int n, Domain *aDomain) : Node(n, aDomain)
{ }


void
Particle :: initializeFrom(InputRecord &ir, int priority)
{
    ParameterManager &ppm =  domain->dofmanPPM;

    Node :: initializeFrom(ir, priority);
    PM_UPDATE_PARAMETER(radius, ppm, ir, this->number, IPK_Particle_rad, priority) ;
}

void Particle::postInitialize () {
    ParameterManager &ppm =  this->giveDomain()->dofmanPPM;

    Node::postInitialize();
    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_Particle_rad) ;
    if ( radius < 0.0 ) {
        throw ComponentInputException(IPK_Particle_rad.getName(), ComponentInputException::ComponentType::ctDofManager, this->number, "must be positive");
    }
 
}

} // namespace oofem
