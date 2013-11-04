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

#include "structuralinterfacecrosssection.h"
#include "structuralinterfacematerialstatus.h"
#include "gausspoint.h"
#include "element.h"
#include "structuralinterfacematerial.h"
#include "floatarray.h"

namespace oofem {

StructuralInterfaceMaterial 
*StructuralInterfaceCrossSection :: giveInterfaceMaterial() 
{
    return static_cast< StructuralInterfaceMaterial *>  ( this->giveDomain()->giveMaterial( this->materialNum ) );
}

const FloatArray 
&StructuralInterfaceCrossSection :: giveTraction(IntegrationPoint *ip) 
{ 
    // Returns the traction vector stored in the material status
    return static_cast< StructuralInterfaceMaterialStatus *> ( this->giveInterfaceMaterial()->giveStatus(ip) )->giveTraction(); 
}

int
StructuralInterfaceCrossSection :: checkConsistency()
{
    // Checks if the given cross section material is a 'StructuralInterfaceMaterial'
    int result = 1;
    Material *mat = this->giveDomain()->giveMaterial( this->materialNum );
    if ( !dynamic_cast< StructuralInterfaceMaterial * >( mat ) ) {
        OOFEM_ERROR2("StructuralInterfaceCrossSection :: checkConsistency : material %s is not a structural interface material", mat->giveClassName());
        result = 0;
    }

    return result;
}

} // end namespace oofem
