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

#include "linearelasticmaterial.h"
#include "gausspoint.h"
#include "simplecrosssection.h"
#include "structuralms.h"

namespace oofem {
void
LinearElasticMaterial :: giveRealStressVector(FloatArray &answer,
                                              GaussPoint *gp,
                                              const FloatArray &reducedStrain,
                                              TimeStep *tStep)
{
    FloatArray strainVector;
    FloatMatrix d;
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(strainVector, gp,
                                                reducedStrain,
                                                tStep, VM_Total);

    this->giveStiffnessMatrix(d, TangentStiffness, gp, tStep);
    answer.beProductOf(d, strainVector);

    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
}


void
LinearElasticMaterial :: giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    FloatArray strainVector;
    FloatMatrix d;
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Total);

    this->give3dMaterialStiffnessMatrix(d, TangentStiffness, gp, tStep);
    answer.beProductOf(d, strainVector);

    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
}


void
LinearElasticMaterial :: giveRealStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    FloatArray strainVector;
    FloatMatrix d;
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Total);

    this->givePlaneStrainStiffMtrx(d, TangentStiffness, gp, tStep);
    answer.beProductOf(d, strainVector);

    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
}


void
LinearElasticMaterial :: giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    FloatArray strainVector;
    FloatMatrix d;
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Total);

    this->givePlaneStressStiffMtrx(d, TangentStiffness, gp, tStep);
    answer.beProductOf(d, strainVector);

    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
}


void
LinearElasticMaterial :: giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    FloatArray strainVector;
    FloatMatrix d;
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Total);

    this->give1dStressStiffMtrx(d, TangentStiffness, gp, tStep);
    answer.beProductOf(d, strainVector);

    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
}


void
LinearElasticMaterial :: giveRealStressVector_2dBeamLayer(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    FloatArray strainVector;
    FloatMatrix d;
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Total);

    this->give2dBeamLayerStiffMtrx(d, TangentStiffness, gp, tStep);
    answer.beProductOf(d, strainVector);

    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
}


void
LinearElasticMaterial :: giveRealStressVector_PlateLayer(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    FloatArray strainVector;
    FloatMatrix d;
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Total);

    this->givePlateLayerStiffMtrx(d, TangentStiffness, gp, tStep);
    answer.beProductOf(d, strainVector);

    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
}


void
LinearElasticMaterial :: giveRealStressVector_Fiber(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    FloatArray strainVector;
    FloatMatrix d;
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Total);

    this->giveFiberStiffMtrx(d, TangentStiffness, gp, tStep);
    answer.beProductOf(d, strainVector);

    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
}
} // end namespace oofem
