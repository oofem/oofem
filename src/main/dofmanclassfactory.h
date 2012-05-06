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

#include "elementside.h"
#include "rigidarmnode.h"
#include "hangingnode.h"
#include "slavenode.h"
#include "node.h"

#ifdef __SM_MODULE
 #include "particle.h"
#endif

REGISTER_CLASS(Node, "node", NodeClass)
REGISTER_CLASS(ElementSide, "elementside", ElementSideClass)
REGISTER_CLASS(RigidArmNode, "rigidarmnode", RigidArmNodeClass)
REGISTER_CLASS(HangingNode, "hangingnode", HangingNodeClass)
REGISTER_CLASS(SlaveNode, "slavenode", SlaveNodeClass)
#ifdef __SM_MODULE
REGISTER_CLASS(Particle, "particle", ParticleClass)
#endif

