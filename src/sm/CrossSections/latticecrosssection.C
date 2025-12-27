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

#include "latticecrosssection.h"
#include "sm/Materials/LatticeMaterials/latticematstatus.h"
#include "sm/Materials/LatticeMaterials/latticestructuralmaterial.h"
#include "gausspoint.h"
#include "element.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"

namespace oofem {
REGISTER_CrossSection(LatticeCrossSection);

LatticeStructuralMaterial *
LatticeCrossSection :: giveLatticeMaterial() const
{
    return static_cast< LatticeStructuralMaterial * >( this->giveDomain()->giveMaterial(this->materialNum) );
}

int
LatticeCrossSection :: checkConsistency()
{
    // Checks if the given cross section material is a 'LatticeStructuralMaterial'
    Material *mat = this->giveDomain()->giveMaterial(this->materialNum);
    if ( !dynamic_cast< LatticeStructuralMaterial * >( mat ) ) {
        OOFEM_WARNING("material %s is not a structural interface material", mat->giveClassName() );
        return 0;
    }

    return 1;
}

void
LatticeCrossSection :: initializeFrom(InputRecord &ir)
{
    CrossSection :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, this->materialNum, _IFT_LatticeCrossSection_Material);

    double thickness = 0.0;
    if ( ir.hasField(_IFT_LatticeCrossSection_thickness) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, thickness, _IFT_LatticeCrossSection_thickness);
        propertyDictionary.add(CS_Thickness, thickness);
    }
}

double LatticeCrossSection :: giveLatticeStress1d(double strain, GaussPoint *gp, TimeStep *tStep) const
{
    return this->giveLatticeMaterial()->giveLatticeStress1d(strain, gp, tStep);
}

FloatArrayF< 3 >LatticeCrossSection :: giveLatticeStress2d(const FloatArrayF< 3 > &strain, GaussPoint *gp, TimeStep *tStep) const
{
    return this->giveLatticeMaterial()->giveLatticeStress2d(strain, gp, tStep);
}

FloatArrayF< 6 >LatticeCrossSection :: giveLatticeStress3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const
{
    return this->giveLatticeMaterial()->giveLatticeStress3d(strain, gp, tStep);
}



FloatMatrixF< 1, 1 >
LatticeCrossSection :: give1dStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    LatticeStructuralMaterial *mat = this->giveLatticeMaterial();
    if ( mat->hasAnalyticalTangentStiffness() ) {
        return mat->give1dLatticeStiffnessMatrix(rMode, gp, tStep);
    } else {
        OOFEM_ERROR("not implemented");
    }
}

FloatMatrixF< 3, 3 >
LatticeCrossSection :: give2dStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    LatticeStructuralMaterial *mat = this->giveLatticeMaterial();
    if ( mat->hasAnalyticalTangentStiffness() ) {
        return mat->give2dLatticeStiffnessMatrix(rMode, gp, tStep);
    } else {
        OOFEM_ERROR("not implemented");
    }
}

FloatMatrixF< 6, 6 >
LatticeCrossSection :: give3dStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    LatticeStructuralMaterial *mat = this->giveLatticeMaterial();
    if ( mat->hasAnalyticalTangentStiffness() ) {
        return mat->give3dLatticeStiffnessMatrix(rMode, gp, tStep);
    } else {
        OOFEM_ERROR("not implemented");
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
    return this->giveLatticeMaterial()->giveIPValue(answer, ip, type, tStep);
}

double
LatticeCrossSection :: give(int aProperty, GaussPoint *gp) const
{
    return this->giveMaterial(gp)->give(aProperty, gp);
}

Material *LatticeCrossSection :: giveMaterial(IntegrationPoint *ip) const
{
    if ( this->giveMaterialNumber() ) {
        return this->giveDomain()->giveMaterial(this->giveMaterialNumber() );
    } else {
        return ip->giveElement()->giveMaterial();
    }
}


int
LatticeCrossSection :: packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->giveLatticeMaterial()->packUnknowns(buff, tStep, gp);
}

int
LatticeCrossSection :: unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->giveLatticeMaterial()->unpackAndUpdateUnknowns(buff, tStep, gp);
}

int
LatticeCrossSection :: estimatePackSize(DataStream &buff, GaussPoint *gp)
{
    return this->giveLatticeMaterial()->estimatePackSize(buff, gp);
}
} // end namespace oofem
