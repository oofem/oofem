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
#include "../sm/CrossSections/simplecrosssection.h"
#include "../sm/Materials/structuralms.h"

namespace oofem {

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
LinearElasticMaterial :: giveRealStressVector_3dDegeneratedShell(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    FloatArray strainVector;
    FloatMatrix d;
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Total);

    this->give3dMaterialStiffnessMatrix(d, TangentStiffness, gp, tStep);

    d.at(1, 1) -= d.at(1, 3) * d.at(3, 1) / d.at(3, 3);
    d.at(2, 1) -= d.at(2, 3) * d.at(3, 1) / d.at(3, 3);
    d.at(1, 2) -= d.at(1, 3) * d.at(3, 2) / d.at(3, 3);
    d.at(2, 2) -= d.at(2, 3) * d.at(3, 2) / d.at(3, 3);

    d.at(3, 1) = 0.0;
    d.at(3, 2) = 0.0;
    d.at(3, 3) = 0.0;
    d.at(2, 3) = 0.0;
    d.at(1, 3) = 0.0;


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
LinearElasticMaterial :: giveRealStressVector_Warping(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    // reducedStrain contains the stress components tau_zy and tau_zx, computed as    //        tau_zy = G * theta * ( x + dPsi/dy )
    //        tau_zx = G * theta * (-y + dPsi/dx )
    // where x and y are the global coordinates of the Gauss point (the origin must be at the centroid of the cross section)
    //       G is the shear modulus of elasticity and theta is the relative twist (dPhi_z/dz)
    FloatArray gcoords;
    Element *elem = gp->giveElement();
    elem->computeGlobalCoordinates( gcoords, gp->giveNaturalCoordinates() );
    answer.resize(2);
    answer.at(1) = reducedStrain.at(1);
    answer.at(2) = reducedStrain.at(2);
    double G = this->giveShearModulus();
    answer.times(G);

    // update gp
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
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

void
LinearElasticMaterial :: giveEshelbyStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedF, TimeStep *tStep)
{
    FloatArray fullFv;
    StructuralMaterial :: giveFullVectorFormF(fullFv, reducedF, _PlaneStrain);

    FloatMatrix H;
    H.beMatrixForm(fullFv);

    H(0, 0) -= 1.0;
    H(1, 1) -= 1.0;
    H(2, 2) -= 1.0;


    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    const FloatArray &stressV = status->giveTempStressVector();

    FloatArray fullStressV;
    StructuralMaterial :: giveFullSymVectorForm(fullStressV, stressV, _PlaneStrain);

    FloatMatrix stress;
    stress.beMatrixFormOfStress(fullStressV);


    double energyDens = giveEnergyDensity(gp, tStep);


    FloatMatrix eshelbyStress(3, 3);
    eshelbyStress.beTProductOf(H, stress);
    eshelbyStress.negated();

    eshelbyStress(0, 0) += energyDens;
    eshelbyStress(1, 1) += energyDens;
    eshelbyStress(2, 2) += energyDens;


    FloatArray eshelbyStressV;
    eshelbyStressV.beVectorForm(eshelbyStress);
    StructuralMaterial :: giveReducedVectorForm(answer, eshelbyStressV, _PlaneStrain);
}

double
LinearElasticMaterial :: giveEnergyDensity(GaussPoint *gp, TimeStep *tStep)
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    const FloatArray &strain = status->giveTempStrainVector();
    const FloatArray &stress = status->giveTempStressVector();

    return 0.5 * stress.dotProduct(strain);
}


MaterialStatus *
LinearElasticMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(1, this->giveDomain(), gp);
}
} // end namespace oofem
