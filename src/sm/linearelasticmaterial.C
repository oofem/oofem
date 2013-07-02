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

#include "linearelasticmaterial.h"
#include "gausspoint.h"
#include "simplecrosssection.h"
#include "structuralms.h"

namespace oofem {
int
LinearElasticMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    if ( ( mode == _3dMat ) || ( mode == _PlaneStress ) ||
        ( mode == _PlaneStrain ) || ( mode == _1dMat ) ||
        ( mode == _2dPlateLayer ) || ( mode == _2dBeamLayer ) ||
        ( mode == _3dShellLayer ) || ( mode == _2dPlate ) ||
        ( mode == _3dShell ) || ( mode == _PlaneStressRot ) ) {
        return 1;
    }

    return 0;
}


void
LinearElasticMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                  MatResponseMode rMode,
                                                  GaussPoint *gp,
                                                  TimeStep *atTime)
{
    //
    // Returns characteristic material stiffness matrix of the receiver
    //
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _2dPlate:
        this->give2dPlateStiffMtrx(answer, rMode, gp, atTime);
        break;
    case _PlaneStressRot:
        this->give2dPlaneStressRotStiffMtrx(answer, rMode, gp, atTime);
        break;
    case _3dShell:
        this->give3dShellStiffMtrx(answer, rMode, gp, atTime);
        break;
    default:
        StructuralMaterial :: giveStiffnessMatrix(answer, rMode, gp, atTime);
    }
}


void
LinearElasticMaterial :: give2dPlateStiffMtrx(FloatMatrix &answer,
                                              MatResponseMode rMode,
                                              GaussPoint *gp,
                                              TimeStep *tStep)
//
// return material stiffness matrix for derived types of stressStreinState
//
{
    MaterialMode mode = gp->giveMaterialMode();
    SimpleCrossSection *crossSection =  dynamic_cast< SimpleCrossSection * >( gp->giveCrossSection() );
    FloatMatrix mat3d;
    double thickness3, thickness;
    int i, j;

    if ( mode != _2dPlate ) {
        _error("Give2dPlateStiffMtrx : unsupported mode");
    }

    if ( crossSection == NULL ) {
        _error(" Give2dPlateStiffMtrx : no SimpleCrossSection");
    }

    this->givePlaneStressStiffMtrx(mat3d, rMode, gp, tStep);
    thickness = crossSection->give(CS_Thickness);
    thickness3 = thickness * thickness * thickness;

    answer.resize(5, 5);
    answer.zero();

    for ( i = 1; i <= 2; i++ ) {
        for ( j = 1; j <= 2; j++ ) {
            answer.at(i, j) = mat3d.at(i, j) * thickness3 / 12.;
        }
    }

    answer.at(3, 3) = mat3d.at(3, 3) * thickness3 / 12.;
    answer.at(4, 4) = mat3d.at(3, 3) * thickness * ( 5. / 6. );
    answer.at(5, 5) = answer.at(4, 4);
}


void
LinearElasticMaterial :: give2dPlaneStressRotStiffMtrx(FloatMatrix &answer,
                                                       MatResponseMode rMode,
                                                       GaussPoint *gp,
                                                       TimeStep *tStep)
//
// return material stiffness matrix for derived types of stressStreinState
//
{
    MaterialMode mode = gp->giveMaterialMode();
    FloatMatrix mat;

    if ( mode != _PlaneStressRot ) {
        _error("Give2dPlaneStressRotStiffMtrx : unsupported mode");
    }

    this->givePlaneStressStiffMtrx(mat, rMode, gp, tStep);

    answer.resize(4, 4);
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(i, j) = mat.at(i, j);
        }
    }

    answer.at(4, 4) = mat.at(3, 3);
}


void
LinearElasticMaterial :: give3dShellStiffMtrx(FloatMatrix &answer,
                                              MatResponseMode rMode,
                                              GaussPoint *gp,
                                              TimeStep *tStep)
//
// return material stiffness matrix for derived types of stressStreinState
//
{
    MaterialMode mode = gp->giveMaterialMode();
    SimpleCrossSection *crossSection = dynamic_cast< SimpleCrossSection * >( gp->giveCrossSection() );
    FloatMatrix mat3d;
    double thickness3, thickness;

    if ( mode != _3dShell ) {
        _error("Give3dShellMaterialStiffness : unsupported mode");
    }

    if ( crossSection == NULL ) {
        _error(" Give2dBeamStiffMtrx : no SimpleCrossSection");
    }

    this->givePlaneStressStiffMtrx(mat3d, rMode, gp, tStep);
    thickness = crossSection->give(CS_Thickness);
    thickness3 = thickness * thickness * thickness;

    answer.resize(8, 8);
    answer.zero();


    for ( int i = 1; i <= 2; i++ ) {
        for ( int j = 1; j <= 2; j++ ) {
            answer.at(i, j) = mat3d.at(i, j) * thickness;
            answer.at(i + 3, j + 3) = mat3d.at(i, j) * thickness3 / 12.0;
        }
    }

    answer.at(3, 1) = mat3d.at(3, 1) * thickness;
    answer.at(3, 2) = mat3d.at(3, 2) * thickness;
    answer.at(3, 3) = mat3d.at(3, 3) * thickness;
    answer.at(6, 4) = mat3d.at(3, 1) * thickness3 / 12.0;
    answer.at(6, 5) = mat3d.at(3, 2) * thickness3 / 12.0;
    answer.at(6, 6) = mat3d.at(3, 3) * thickness3 / 12.0;

    answer.at(7, 7) = mat3d.at(3, 3) * thickness * ( 5. / 6. );
    answer.at(8, 8) = answer.at(7, 7);
}


void
LinearElasticMaterial :: giveRealStressVector(FloatArray &answer,
                                              GaussPoint *gp,
                                              const FloatArray &reducedStrain,
                                              TimeStep *atTime)
{
    FloatArray stressIncrement, stressVector, strainVector;
    FloatMatrix d;
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(strainVector, gp,
                                                reducedStrain,
                                                atTime, VM_Total);

    if ( status->giveStressVector().giveSize() ) {
        stressVector = status->giveStressVector();
    } else {
        stressVector.resize( strainVector.giveSize() );
    }

    this->giveCharacteristicMatrix(d, TangentStiffness, gp, atTime);
    stressVector.beProductOf(d, strainVector);

    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(stressVector);

    answer = stressVector;
}

} // end namespace oofem
