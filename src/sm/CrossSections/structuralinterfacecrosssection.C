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

#include "../sm/CrossSections/structuralinterfacecrosssection.h"
#include "../sm/Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"
#include "../sm/Materials/InterfaceMaterials/structuralinterfacematerial.h"
#include "gausspoint.h"
#include "element.h"
#include "floatarray.h"

namespace oofem {
REGISTER_CrossSection(StructuralInterfaceCrossSection);

StructuralInterfaceMaterial *
StructuralInterfaceCrossSection :: giveInterfaceMaterial()
{
    return static_cast< StructuralInterfaceMaterial * >( this->giveDomain()->giveMaterial(this->materialNum) );
}

const FloatArray &
StructuralInterfaceCrossSection :: giveTraction(IntegrationPoint *ip)
{
    // Returns the traction vector stored in the material status
    //return static_cast< StructuralInterfaceMaterialStatus *> ( this->giveInterfaceMaterial()->giveStatus(ip) )->giveTraction();
    return static_cast< StructuralInterfaceMaterialStatus * >( this->giveInterfaceMaterial()->giveStatus(ip) )->giveFirstPKTraction();
}


int
StructuralInterfaceCrossSection :: checkConsistency()
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
StructuralInterfaceCrossSection :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    CrossSection :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, this->materialNum, _IFT_StructuralInterfaceCrossSection_Material);

    double thickness = 0.0;
    if ( ir->hasField(_IFT_StructuralInterfaceCrossSection_thickness) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, thickness, _IFT_StructuralInterfaceCrossSection_thickness);
        propertyDictionary.add(CS_Thickness, thickness);
    }

    return IRRT_OK;
}


void
StructuralInterfaceCrossSection :: give1dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralInterfaceMaterial *mat = this->giveInterfaceMaterial();
    if ( mat->useNumericalTangent ) {
        mat->give1dStiffnessMatrix_Eng_Num( answer, gp, tStep );
    } else if( mat->hasAnalyticalTangentStiffness() ) {
        mat->give1dStiffnessMatrix_Eng( answer, rMode, gp, tStep );
    } else {
        OOFEM_ERROR("Not implemented - use numerical tangent instead (keyword: 'use_num_tangent') ");
    }
}

void
StructuralInterfaceCrossSection :: give2dStiffnessMatrix_Eng(FloatMatrix &answer,  MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralInterfaceMaterial *mat = this->giveInterfaceMaterial();
    if ( mat->useNumericalTangent ) {
        mat->give2dStiffnessMatrix_Eng_Num( answer, gp, tStep );
    } else if ( mat->hasAnalyticalTangentStiffness() ) {
        mat->give2dStiffnessMatrix_Eng(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR("not implemented - use numerical tangent instead (keyword: 'use_num_tangent') ");
    }
}

void
StructuralInterfaceCrossSection :: give3dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralInterfaceMaterial *mat = this->giveInterfaceMaterial( );
    if ( mat->useNumericalTangent ) {
        mat->give3dStiffnessMatrix_Eng_Num(answer, gp, tStep);
    } else if ( mat->hasAnalyticalTangentStiffness() ) {
        mat->give3dStiffnessMatrix_Eng(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR("Not implemented - use numerical tangent instead (keyword: 'use_num_tangent') ");
    }
}




void
StructuralInterfaceCrossSection :: give1dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralInterfaceMaterial *mat = this->giveInterfaceMaterial();
    if ( mat->useNumericalTangent ) {
        mat->give1dStiffnessMatrix_dTdj_Num( answer, gp, tStep);
    } else if ( mat->hasAnalyticalTangentStiffness() ) {
        mat->give1dStiffnessMatrix_dTdj(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR("not implemented - use numerical tangent instead (keyword: 'use_num_tangent') ");
    }
}

void
StructuralInterfaceCrossSection :: give2dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralInterfaceMaterial *mat = this->giveInterfaceMaterial();
    if ( mat->useNumericalTangent ) {
        mat->give2dStiffnessMatrix_dTdj_Num( answer, gp, tStep );
    } else if ( mat->hasAnalyticalTangentStiffness() ) {
        mat->give2dStiffnessMatrix_dTdj(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR("not implemented - use numerical tangent instead (keyword: 'use_num_tangent') ");
    }
}


void
StructuralInterfaceCrossSection :: give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralInterfaceMaterial *mat = this->giveInterfaceMaterial();
    if ( mat->useNumericalTangent ) {
        mat->give3dStiffnessMatrix_dTdj_Num(answer, gp, tStep);
    } else if ( mat->hasAnalyticalTangentStiffness() ) {
        mat->give3dStiffnessMatrix_dTdj(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR("not implemented - use numerical tangent instead (keyword: 'use_num_tangent') ");
    }
}


int
StructuralInterfaceCrossSection :: giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_CrossSectionNumber ) {
        answer.resize(1);
        answer.at(1) = this->giveNumber();
        return 1;
    }
    return this->giveInterfaceMaterial()->giveIPValue(answer, ip, type, tStep);
}

int
StructuralInterfaceCrossSection :: packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->giveInterfaceMaterial()->packUnknowns(buff, tStep, gp);
}

int
StructuralInterfaceCrossSection :: unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->giveInterfaceMaterial()->unpackAndUpdateUnknowns(buff, tStep, gp);
}

int
StructuralInterfaceCrossSection :: estimatePackSize(DataStream &buff, GaussPoint *gp)
{
    return this->giveInterfaceMaterial()->estimatePackSize(buff, gp);
}

} // end namespace oofem
