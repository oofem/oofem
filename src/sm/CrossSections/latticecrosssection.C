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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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

#include "sm/CrossSections/latticestructuralcrosssection.h"
#include "sm/Materials/LatticeMaterials/latticematstatus.h"
#include "sm/Materials/LatticeMaterials/latticestructuralmaterial.h"
#include "gausspoint.h"
#include "element.h"
#include "floatarray.h"

namespace oofem {
REGISTER_CrossSection(LatticeCrossSection);

StructuralInterfaceMaterial *
LatticeCrossSection :: giveLatticeMaterial() const
{
    return static_cast< LatticeStructuralMaterial * >( this->giveDomain()->giveMaterial(this->materialNum) );
}

int
LatticeCrossSection :: checkConsistency()
{
    // Checks if the given cross section material is a 'StructuralInterfaceMaterial'
    Material *mat = this->giveDomain()->giveMaterial(this->materialNum);
    if ( !dynamic_cast< StructuralInterfaceMaterial * >(mat) ) {
        OOFEM_WARNING("material %s is not a structural interface material", mat->giveClassName() );
        return 0;
    }

    return 1;
}

IRResultType
LatticeCrossSection :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    CrossSection :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, this->materialNum, _IFT_LatticeCrossSection_Material);

    double thickness = 0.0;
    if ( ir->hasField(_IFT_LatticeCrossSection_thickness) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, thickness, _IFT_LatticeCrossSection_thickness);
        propertyDictionary.add(CS_Thickness, thickness);
    }

    return IRRT_OK;
}


void
LatticeStructuralCrossSection :: giveCharMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
  LatticeMaterial *mat = dynamic_cast< LatticeMaterial * >( this->giveMaterial(gp) );
  if ( mode == _3dLattice ) {
    mat->give3dLatticeStiffnessMatrix(answer, rMode, gp, tStep);
  } else if ( mode == _2dLattice ) {
    mat->give2dLatticeStiffMtrx(answer, rMode, gp, tStep);
  } else if ( mode == _1dLattice ) {
    mat->give1dLatticeStiffMtrx(answer, rMode, gp, tStep);
  }  else {
    mat->giveStiffnessMatrix(answer, rMode, gp, tStep);
  }
}

int
LatticeCrossSection :: giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_CrossSectionNumber ) {
        answer.resize(1);
        answer.at(1) = this->giveNumber();
        return 1;
    }
    return this->giveInterfaceMaterial()->giveIPValue(answer, ip, type, tStep);
}


Material *LatticeCrossSection :: giveMaterial(IntegrationPoint *ip)
{
    if ( this->giveMaterialNumber() ) {
        return this->giveDomain()->giveMaterial( this->giveMaterialNumber() );
    } else {
        return ip->giveElement()->giveMaterial();
    }
}


int
LatticeCrossSection :: packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->giveInterfaceMaterial()->packUnknowns(buff, tStep, gp);
}

int
LatticeCrossSection :: unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->giveInterfaceMaterial()->unpackAndUpdateUnknowns(buff, tStep, gp);
}

int
LatticeCrossSection :: estimatePackSize(DataStream &buff, GaussPoint *gp)
{
    return this->giveInterfaceMaterial()->estimatePackSize(buff, gp);
}

} // end namespace oofem
