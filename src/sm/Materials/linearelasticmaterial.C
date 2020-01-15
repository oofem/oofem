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
#include "sm/CrossSections/simplecrosssection.h"
#include "sm/Materials/structuralms.h"
#include "dynamicinputrecord.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"

namespace oofem {
void
LinearElasticMaterial :: initializeFrom(InputRecord &ir)
{
    StructuralMaterial :: initializeFrom(ir);

    preCastStiffnessReduction = 0.99999999;
    IR_GIVE_OPTIONAL_FIELD(ir, preCastStiffnessReduction, _IFT_LinearElasticMaterial_preCastStiffRed);
}


void
LinearElasticMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);
    input.setField(this->preCastStiffnessReduction, _IFT_LinearElasticMaterial_preCastStiffRed);
}

void
LinearElasticMaterial :: computesSubTangents()
{
    tangentPlaneStrain = {
        tangent(0, 0), tangent(1, 0), tangent(2, 0), tangent(3, 0),
        tangent(0, 1), tangent(1, 1), tangent(2, 1), tangent(3, 1),
        tangent(0, 2), tangent(1, 2), tangent(2, 2), tangent(3, 2),
        tangent(0, 3), tangent(1, 3), tangent(2, 3), tangent(3, 3),
    };

    auto c = inv(tangent);
    FloatMatrixF<3,3> reduced = {
        c(0,0), c(0,1), c(0,5),
        c(1,0), c(1,1), c(1,5),
        c(5,0), c(5,1), c(5,5),
    };
    tangentPlaneStress = inv(reduced);
}


FloatMatrixF<6,6>
LinearElasticMaterial :: give3dMaterialStiffnessMatrix(MatResponseMode mode,
                                                       GaussPoint *gp,
                                                       TimeStep *tStep) const
{
    if ( tStep->giveIntrinsicTime() < this->castingTime ) {
        return tangent *  (1. - this->preCastStiffnessReduction);
    } else {
        return tangent;
    }
}


FloatArrayF<6>
LinearElasticMaterial :: giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const
{
    return alpha;
}


FloatArrayF<6>
LinearElasticMaterial :: giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    auto d = this->give3dMaterialStiffnessMatrix(TangentStiffness, gp, tStep);

    FloatArrayF<6> stress;
    if ( this->castingTime < 0. ) { // no changes in material stiffness ->> total formulation
        auto thermalStrain = this->computeStressIndependentStrainVector_3d(gp, tStep, VM_Total);
        auto strainVector = strain - thermalStrain;
        stress = dot(d, strainVector);
    } else { // changes in material stiffness ->> incremental formulation
        auto thermalStrain = this->computeStressIndependentStrainVector_3d(gp, tStep, VM_Incremental);
        auto strainIncrement = strain - thermalStrain - FloatArrayF<6>(status->giveStrainVector());

        stress = dot(d, strainIncrement) + status->giveStressVector();
    }

    // update gp
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(stress);
    return stress;
}


FloatArrayF<6>
LinearElasticMaterial :: giveRealStressVector_3dDegeneratedShell(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    auto d = this->give3dMaterialStiffnessMatrix(TangentStiffness, gp, tStep);

    d.at(1, 1) -= d.at(1, 3) * d.at(3, 1) / d.at(3, 3);
    d.at(2, 1) -= d.at(2, 3) * d.at(3, 1) / d.at(3, 3);
    d.at(1, 2) -= d.at(1, 3) * d.at(3, 2) / d.at(3, 3);
    d.at(2, 2) -= d.at(2, 3) * d.at(3, 2) / d.at(3, 3);

    d.at(3, 1) = 0.0;
    d.at(3, 2) = 0.0;
    d.at(3, 3) = 0.0;
    d.at(2, 3) = 0.0;
    d.at(1, 3) = 0.0;

    FloatArrayF<6> stress, strainVector;
    if ( this->castingTime < 0. ) { // no changes in material stiffness ->> total formulation
        auto thermalStrain = this->computeStressIndependentStrainVector_3d(gp, tStep, VM_Total);
        strainVector = strain - thermalStrain;
        stress = dot(d, strainVector);
    } else { // changes in material stiffness ->> incremental formulation
        auto thermalStrain = this->computeStressIndependentStrainVector_3d(gp, tStep, VM_Incremental);
        auto strainIncrement = strain - thermalStrain - FloatArrayF<6>(status->giveStrainVector());

        stress = dot(d, strainIncrement) + status->giveStressVector();
    }

    stress = dot(d, strainVector);

    // update gp
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(stress);
    return stress;
}


FloatArrayF<3>
LinearElasticMaterial :: giveRealStressVector_PlaneStress(const FloatArrayF<3> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    auto d = this->givePlaneStressStiffMtrx(TangentStiffness, gp, tStep);

    FloatArray answer;
    FloatArray strainVector;
    if ( this->castingTime < 0. ) { // no changes in material stiffness ->> total formulation
        this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Total);
        answer.beProductOf(d, strainVector);
    } else { // changes in material stiffness ->> incremental formulation
        this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Incremental);
        auto strainIncrement = strainVector - status->giveStrainVector();

        answer.beProductOf(d, strainIncrement);
        answer.add( status->giveStressVector() );
    }

    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
    return answer;
}


FloatArrayF<1>
LinearElasticMaterial :: giveRealStressVector_1d(const FloatArrayF<1> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    auto d = this->give1dStressStiffMtrx(TangentStiffness, gp, tStep);

    FloatArray answer;
    FloatArray strainVector;
    if ( this->castingTime < 0. ) {  // no changes in material stiffness ->> total formulation
        this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Total);
        answer.beProductOf(d, strainVector);
    } else { // changes in material stiffness ->> incremental formulation
        this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Incremental);
        auto strainIncrement = strainVector - status->giveStrainVector();

        answer.beProductOf(d, strainIncrement);
        answer.add( status->giveStressVector() );
    }

    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
    return answer;
}


FloatArrayF<2>
LinearElasticMaterial :: giveRealStressVector_Warping(const FloatArrayF<2> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const
{
    // reducedStrain contains the stress components tau_zy and tau_zx, computed as    //        tau_zy = G * theta * ( x + dPsi/dy )
    //        tau_zx = G * theta * (-y + dPsi/dx )
    // where x and y are the global coordinates of the Gauss point (the origin must be at the centroid of the cross section)
    //       G is the shear modulus of elasticity and theta is the relative twist (dPhi_z/dz)
    
    ///@todo Why is this warping method implemented here? It seems to assume isotropic linear elastic material, so it should be implemented in the subclass instead of relying on "giveShearModulus" which is only implemented there anyway.
    
    auto status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    double G = this->giveShearModulus();
    FloatArray gcoords;
    Element *elem = gp->giveElement();
    elem->computeGlobalCoordinates( gcoords, gp->giveNaturalCoordinates() );
    FloatArrayF<2> answer;

    if ( this->castingTime < 0. ) { // no changes in material stiffness ->> total formulation
        answer.at(1) = reducedStrain.at(1);
        answer.at(2) = reducedStrain.at(2);

        answer *= G;
    } else { // changes in material stiffness ->> incremental formulation
        FloatArray strainIncrement;
        strainIncrement.beDifferenceOf( reducedStrain, status->giveStrainVector() );

        answer.at(1) = strainIncrement.at(1);
        answer.at(2) = strainIncrement.at(2);
        answer *= G;

        if ( ( tStep->giveIntrinsicTime() < this->castingTime ) ) {
            answer *= (1. - this->preCastStiffnessReduction);
        }

        answer += FloatArrayF<2>(status->giveStressVector());
    }

    // update gp

    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
    return answer;
}


FloatArrayF<2>
LinearElasticMaterial :: giveRealStressVector_2dBeamLayer(const FloatArrayF<2> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    auto d = this->give2dBeamLayerStiffMtrx(TangentStiffness, gp, tStep);

    FloatArray answer;
    FloatArray strainVector, strainIncrement;
    if ( this->castingTime < 0. ) { // no changes in material stiffness ->> total formulation
        this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Total);
        answer.beProductOf(d, strainVector);
    } else { // changes in material stiffness ->> incremental formulation
        this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Incremental);
        strainIncrement.beDifferenceOf( strainVector, status->giveStrainVector() );

        answer.beProductOf(d, strainIncrement);
        answer.add( status->giveStressVector() );
    }

    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
    return answer;
}


FloatArrayF<5>
LinearElasticMaterial :: giveRealStressVector_PlateLayer(const FloatArrayF<5> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    auto d = this->givePlateLayerStiffMtrx(TangentStiffness, gp, tStep);

    FloatArray answer;
    FloatArray strainVector, strainIncrement;
    if ( this->castingTime < 0. ) { // no changes in material stiffness ->> total formulation
        this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Total);
        answer.beProductOf(d, strainVector);
    } else { // changes in material stiffness ->> incremental formulation
        this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Incremental);
        strainIncrement.beDifferenceOf( strainVector, status->giveStrainVector() );

        answer.beProductOf(d, strainIncrement);
        answer.add( status->giveStressVector() );
    }

    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
    return answer;
}


FloatArrayF<3>
LinearElasticMaterial :: giveRealStressVector_Fiber(const FloatArrayF<3> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    auto d = this->giveFiberStiffMtrx(TangentStiffness, gp, tStep);

    FloatArray answer;
    FloatArray strainVector;
    if ( this->castingTime < 0. ) { // no changes in material stiffness ->> total formulation
        this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Total);
        answer.beProductOf(d, strainVector);
    } else { // changes in material stiffness ->> incremental formulation
        this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Incremental);
        auto strainIncrement = strainVector - status->giveStrainVector();

        answer.beProductOf(d, strainIncrement);
        answer.add( status->giveStressVector() );
    }

    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
    return answer;
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

int
LinearElasticMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_ElasticStrainTensor ) {
        this->giveStressDependentPartOfStrainVector(answer, gp, status->giveStrainVector(), tStep, VM_Total);
        return 1;
    } else if ( type == IST_ThermalStrainTensor ) {
        answer = this->computeStressIndependentStrainVector_3d(gp, tStep, VM_Total);
        return 1;
    } else if ( type == IST_CreepStrainTensor ) {
        answer.resize(6);
        answer.zero();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}


MaterialStatus *
LinearElasticMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(gp);
}
} // end namespace oofem
