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
#include "isolinearelasticmaterial.h"
#include "../sm/CrossSections/simplecrosssection.h"
#include "../sm/Materials/structuralms.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(IsotropicLinearElasticMaterial);

IsotropicLinearElasticMaterial :: IsotropicLinearElasticMaterial(int n, Domain *d,
                                                                 double _E, double _nu) :
    LinearElasticMaterial(n, d)
{
    E = _E;
    nu = _nu;

    // compute  value of shear modulus
    G = E / ( 2.0 * ( 1. + nu ) );
}


IRResultType
IsotropicLinearElasticMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    double value;
    // we use rather object's member data than to store data into slow
    // key-val dictionary with lot of memory allocations

    IR_GIVE_FIELD(ir, E, _IFT_IsotropicLinearElasticMaterial_e);
    IR_GIVE_FIELD(ir, nu, _IFT_IsotropicLinearElasticMaterial_n);
    IR_GIVE_FIELD(ir, value, _IFT_IsotropicLinearElasticMaterial_talpha);
    propertyDictionary.add(tAlpha, value);
    // compute  value of shear modulus
    G = E / ( 2.0 * ( 1. + nu ) );

    return LinearElasticMaterial :: initializeFrom(ir);;
}


void
IsotropicLinearElasticMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    this->LinearElasticMaterial :: giveInputRecord(input);
    StructuralMaterial :: giveInputRecord(input);

    input.setField(this->E, _IFT_IsotropicLinearElasticMaterial_e);
    input.setField(this->nu, _IFT_IsotropicLinearElasticMaterial_n);
    input.setField(this->propertyDictionary.at(tAlpha), _IFT_IsotropicLinearElasticMaterial_talpha);
}



double
IsotropicLinearElasticMaterial :: give(int aProperty, GaussPoint *gp)
//
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
//
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


void
IsotropicLinearElasticMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                                MatResponseMode mode,
                                                                GaussPoint *gp,
                                                                TimeStep *tStep)
//
// forceElasticResponse ignored - always elastic
//
{
    double e, nu, ee;

    e  = this->E;
    nu = this->nu;

    ee = e / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
    answer.resize(6, 6);
    answer.zero();

    answer.at(1, 1) =  1. - nu;
    answer.at(1, 2) =  nu;
    answer.at(1, 3) =  nu;
    answer.at(2, 1) =  nu;
    answer.at(2, 2) =  1. - nu;
    answer.at(2, 3) =  nu;
    answer.at(3, 1) =  nu;
    answer.at(3, 2) =  nu;
    answer.at(3, 3) =  1. - nu;

    answer.at(4, 4) =  ( 1. - 2. * nu ) * 0.5;
    answer.at(5, 5) =  ( 1. - 2. * nu ) * 0.5;
    answer.at(6, 6) =  ( 1. - 2. * nu ) * 0.5;

    answer.times(ee);
}


void
IsotropicLinearElasticMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep)
{
    this->giveStatus(gp);
    double e, nu, ee, shear;

    e     = this->E;
    nu    = this->nu;
    ee    = e / ( 1. - nu * nu );
    shear = this->G;

    answer.resize(3, 3);
    answer.zero();

    answer.at(1, 1) = ee;
    answer.at(1, 2) = nu * ee;
    answer.at(2, 1) = nu * ee;
    answer.at(2, 2) = ee;
    answer.at(3, 3) = shear;
}


void
IsotropicLinearElasticMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep)
{
    double e, nu, ee, shear;

    e     = this->E;
    nu    = this->nu;
    ee    = e / ( 1.0 + nu ) / ( 1. - 2.0 * nu );
    shear = this->G;

    answer.resize(4, 4);
    answer.zero();

    answer.at(1, 1) = ee * ( 1.0 - nu );
    answer.at(1, 2) = nu * ee;
    answer.at(1, 3) = nu * ee;
    answer.at(2, 1) = nu * ee;
    answer.at(2, 2) = ee * ( 1.0 - nu );
    answer.at(2, 3) = nu * ee;
    answer.at(3, 1) = ee * nu;
    answer.at(3, 2) = ee * nu;
    answer.at(3, 3) = ee * ( 1.0 - nu );
    answer.at(4, 4) = shear;
}


void
IsotropicLinearElasticMaterial :: give1dStressStiffMtrx(FloatMatrix &answer,
                                                        MatResponseMode mode,
                                                        GaussPoint *gp,
                                                        TimeStep *tStep)
{
    answer.resize(1, 1);
    answer.at(1, 1) = this->E;
}


void
IsotropicLinearElasticMaterial :: giveThermalDilatationVector(FloatArray &answer,
                                                              GaussPoint *gp,  TimeStep *tStep)
//
// returns a FloatArray(6) of initial strain vector
// eps_0 = {exx_0, eyy_0, ezz_0, gyz_0, gxz_0, gxy_0}^T
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    double alpha = this->give(tAlpha, gp);
    answer.resize(6);
    answer.zero();
    answer.at(1) = alpha;
    answer.at(2) = alpha;
    answer.at(3) = alpha;
}

} // end namespace oofem
