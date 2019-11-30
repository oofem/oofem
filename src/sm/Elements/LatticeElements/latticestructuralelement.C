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

#include "sm/Elements/LatticeElements/latticestructuralelement.h"

namespace oofem {
LatticeStructuralElement :: LatticeStructuralElement(int n, Domain *aDomain) : StructuralElement(n, aDomain)
{ }

void
LatticeStructuralElement :: initializeFrom(InputRecord &ir)
{
    StructuralElement :: initializeFrom(ir);
}

void
LatticeStructuralElement :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralElement :: printOutputAt(file, tStep);

    /// FIXME: This output should just be moved to the elements themselves. But, they don't exist yet? / Mikael
    FloatArray forces;
    if ( this->giveClassName() == std::string("LatticeBeam3d") ) {
        this->giveInternalForcesVector(forces, tStep, 0);
        fprintf(file, "LatticeBeam forces = %e %e %e %e %e %e.\n", forces.at(7), forces.at(8), forces.at(9), forces.at(10), forces.at(11), forces.at(12) );
    } else if ( this->giveClassName() == std::string("LatticeBeam3dBoundary") ) {
        this->giveInternalForcesVector(forces, tStep, 0);
        fprintf(file, "LatticeBeam3dBoundary forces = %e %e %e %e %e %e.\n", forces.at(7), forces.at(8), forces.at(9), forces.at(10), forces.at(11), forces.at(12) );
    }
}

} // end namespace oofem
