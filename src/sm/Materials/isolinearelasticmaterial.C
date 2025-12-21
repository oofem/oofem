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

#include "linearelasticmaterial.h"
#include "isolinearelasticmaterial.h"
#include "sm/CrossSections/simplecrosssection.h"
#include "sm/Materials/structuralms.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "contextioerr.h"

namespace oofem {
REGISTER_Material(IsotropicLinearElasticMaterial);

IsotropicLinearElasticMaterial :: IsotropicLinearElasticMaterial(int n, Domain *d) :
    LinearElasticMaterial(n, d), MixedPressureMaterialExtensionInterface(d)
{ }

IsotropicLinearElasticMaterial :: IsotropicLinearElasticMaterial(int n, Domain *d,
                                                                 double _E, double _nu) :
    LinearElasticMaterial(n, d), MixedPressureMaterialExtensionInterface(d),
    E(_E),
    nu(_nu),
    G( E / ( 2.0 * ( 1. + nu ) ) ),
    a(0.)
{
    this->initTangents();
}


void
IsotropicLinearElasticMaterial :: initializeFrom(InputRecord &ir)
{
    LinearElasticMaterial :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, E, _IFT_IsotropicLinearElasticMaterial_e);
    IR_GIVE_FIELD(ir, nu, _IFT_IsotropicLinearElasticMaterial_n);
    IR_GIVE_FIELD(ir, a, _IFT_IsotropicLinearElasticMaterial_talpha);
    //    propertyDictionary.add(tAlpha, a);
    // compute  value of shear modulus
    G = E / ( 2.0 * ( 1. + nu ) );
    this->initTangents();
}


void
IsotropicLinearElasticMaterial :: initTangents()
{
    double K = E / ( 3.0 * ( 1. - 2. * nu ) );
    this->tangent = 2 * G * I_dev6 + K * I6_I6;
    this->computesSubTangents();
    this->alpha = {
        a, a, a, 0., 0., 0.
    };
}

void
IsotropicLinearElasticMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    this->LinearElasticMaterial :: giveInputRecord(input);
    StructuralMaterial :: giveInputRecord(input);

    input.setField(this->E, _IFT_IsotropicLinearElasticMaterial_e);
    input.setField(this->nu, _IFT_IsotropicLinearElasticMaterial_n);
    input.setField(this->a, _IFT_IsotropicLinearElasticMaterial_talpha);
}


void IsotropicLinearElasticMaterial :: saveContext(DataStream &stream, ContextMode mode)
{
    LinearElasticMaterial :: saveContext(stream, mode);

    if ( ( mode & CM_Definition ) ) {
        if ( !stream.write(E) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(nu) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(a) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
}


void IsotropicLinearElasticMaterial :: restoreContext(DataStream &stream, ContextMode mode)
{
    LinearElasticMaterial :: restoreContext(stream, mode);

    if ( mode & CM_Definition ) {
        if ( !stream.read(E) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(nu) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(a) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
    G = E / ( 2.0 * ( 1. + nu ) );
    this->initTangents();
}



double
IsotropicLinearElasticMaterial :: give(int aProperty, GaussPoint *gp) const
{
    if ( ( aProperty == NYxy ) || ( aProperty == NYxz ) || ( aProperty == NYyz ) ) {
        return nu;
    }

    if ( ( aProperty == 'G' ) || ( aProperty == Gyz ) || ( aProperty == Gxz ) ||
         ( aProperty == Gxy ) ) {
        return G;
    }

    if ( ( aProperty == 'E' ) || ( aProperty == Ex ) || ( aProperty == Ey ) ||
         ( aProperty == Ez ) ) {
        return E;
    }

    if ( ( aProperty == 'n' ) || ( aProperty == NYzx ) || ( aProperty == NYzy ) ||
         ( aProperty == NYyx ) ) {
        return nu;
    }

    return this->Material :: give(aProperty, gp);
}


FloatMatrixF<3,3>
IsotropicLinearElasticMaterial :: givePlaneStressStiffMtrx(MatResponseMode mode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep) const
{
    if ( ( tStep->giveIntrinsicTime() < this->castingTime ) ) {
        return tangentPlaneStress * (1. - this->preCastStiffnessReduction);
    } else {
        return tangentPlaneStress;
    }
}


FloatMatrixF<4,4>
IsotropicLinearElasticMaterial :: givePlaneStrainStiffMtrx(MatResponseMode mode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep) const
{
    if ( ( tStep->giveIntrinsicTime() < this->castingTime ) ) {
        return tangentPlaneStrain * (1. - this->preCastStiffnessReduction);
    } else {
        return tangentPlaneStrain;
    }
}


FloatMatrixF<1,1>
IsotropicLinearElasticMaterial :: give1dStressStiffMtrx(MatResponseMode mode,
                                                        GaussPoint *gp,
                                                        TimeStep *tStep) const
{
    double e = this->E;
    if ( tStep->giveIntrinsicTime() < this->castingTime ) {
        e *= 1. - this->preCastStiffnessReduction;
    }
    return {e};
}








void
IsotropicLinearElasticMaterial :: giveDeviatoric3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                                          MatResponseMode mode,
                                                                          GaussPoint *gp,
                                                                          TimeStep *tStep) const
//
// forceElasticResponse ignored - always elastic
//
{
    double mu = this->E / 2 / ( 1. + this->nu );


    answer.resize(6, 6);
    answer.zero();

    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = 4. / 3.;
    answer.at(1, 2) = answer.at(1, 3) = -2. / 3.;
    answer.at(2, 1) = answer.at(2, 3) = -2. / 3.;
    answer.at(3, 1) = answer.at(3, 2) = -2. / 3.;

    answer.at(4, 4) = answer.at(5, 5) = answer.at(6, 6) = 1.0;

    answer.times(mu);
}


void
IsotropicLinearElasticMaterial :: giveDeviatoricPlaneStrainStiffMtrx(FloatMatrix &answer,
                                                                     MatResponseMode mode,
                                                                     GaussPoint *gp,
                                                                     TimeStep *tStep) const
{
    answer.resize(4, 4);

    double mu = this->E / 2 / ( 1. + this->nu );

    answer.zero();
    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) =  4. / 3.;
    answer.at(1, 2) = answer.at(1, 3) = -2. / 3.;
    answer.at(2, 1) = answer.at(2, 3) = -2. / 3.;
    answer.at(3, 1) = answer.at(3, 2) = -2. / 3.;

    answer.at(4, 4) = 1.0;

    answer.times(mu);
}



void
IsotropicLinearElasticMaterial :: giveRealStressVectorUP_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, double pressure, TimeStep *tStep) const
{
    FloatArray strainVector;
    FloatMatrix d;
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Total);

    this->giveDeviatoric3dMaterialStiffnessMatrix(d, TangentStiffness, gp, tStep);
    answer.beProductOf(d, strainVector);
    answer.at(1) -= pressure;
    answer.at(2) -= pressure;
    answer.at(3) -= pressure;

    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
}



void
IsotropicLinearElasticMaterial :: giveRealStressVectorUP_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, double pressure, TimeStep *tStep) const
{
    FloatArray strainVector;
    FloatMatrix d;
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Total);

    this->giveDeviatoricPlaneStrainStiffMtrx(d, TangentStiffness, gp, tStep);
    answer.beProductOf(d, strainVector);
    answer.at(1) -= pressure;
    answer.at(2) -= pressure;
    answer.at(3) -= pressure;


    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
}

double 
IsotropicLinearElasticMaterial::giveCharacteristicValue(MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const {
    switch (type) {
        case ElasticBulkModulus:
            return giveBulkModulus(); 
        case ElasticBulkModulusInverse:
            return 1./giveBulkModulus();
        case MRM_ScalarOne:
            return 1.0;
        default:
            return this->Material::giveCharacteristicValue(type, gp, tStep);
    }
}

void
IsotropicLinearElasticMaterial::giveCharacteristicVector(FloatArray &answer, FloatArray& flux, MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const {
    if (type == Stress) {
        return LinearElasticMaterial::giveRealStressVector(answer, gp, flux, tStep);
    } else if (type == DeviatoricStress) {
        FloatMatrix d;
        this->giveDeviatoricConstitutiveMatrix(d, TangentStiffness, gp, tStep);
        answer.beProductOf(d, flux);
        return; 
    } else {
        OOFEM_ERROR("Not implemented");
    } 
}

void 
IsotropicLinearElasticMaterial::giveCharacteristicMatrix(FloatMatrix &answer, MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const {
    if (type == DeviatoricStiffness) {
        return this->giveDeviatoricConstitutiveMatrix(answer, type, gp, tStep);
    } else {
        return this->LinearElasticMaterial::giveCharacteristicMatrix(answer, type, gp, tStep);
    }
}

} // end namespace oofem
