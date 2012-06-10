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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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
#include "isolinearelasticmaterial.h"
#include "simplecrosssection.h"
#include "material.h"
#include "structuralms.h"
#include "flotmtrx.h"
#include "gausspnt.h"

namespace oofem {
IsotropicLinearElasticMaterial :: IsotropicLinearElasticMaterial(int n, Domain *d,
                                                                 double _E, double _nu) :
    LinearElasticMaterial(n, d)
{
    E = _E;
    nu = _nu;

    // compute  value of shear modulus
    G = E / ( 2.0 * ( 1. + nu ) );
}

int
IsotropicLinearElasticMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    if ( ( mode == _3dMat ) || ( mode == _PlaneStress ) ||
        ( mode == _PlaneStrain ) || ( mode == _1dMat ) ||
        ( mode == _2dPlateLayer ) || ( mode == _2dBeamLayer ) ||
        ( mode == _3dShellLayer ) || ( mode == _2dPlate ) ||
        ( mode == _2dBeam ) || ( mode == _3dShell ) ||
        ( mode == _3dBeam ) || ( mode == _PlaneStressRot ) ||
        ( mode == _1dFiber ) ) {
        return 1;
    }

    return 0;
}


void
IsotropicLinearElasticMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                           MatResponseForm form,
                                                           MatResponseMode rMode,
                                                           GaussPoint *gp,
                                                           TimeStep *atTime)
{
    //
    // Returns characteristic material stiffness matrix of the receiver
    //
    MaterialMode mMode = gp->giveMaterialMode();

    switch ( mMode ) {
    case _2dBeam:
        this->give2dBeamStiffMtrx(answer, form, rMode, gp, atTime);
        break;
    case _3dBeam:
        this->give3dBeamStiffMtrx(answer, form, rMode, gp, atTime);
        break;
    default:
        LinearElasticMaterial :: giveCharacteristicMatrix(answer, form, rMode, gp, atTime);
    }

    if ( !isActivated(atTime) ) {
        answer.times(0.0);
    }
}



IRResultType
IsotropicLinearElasticMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    double value;

    this->LinearElasticMaterial :: initializeFrom(ir);
    // we use rather object's member data than to store data into slow
    // key-val dictionary with lot of memory allocations

    IR_GIVE_FIELD(ir, E, IFT_IsotropicLinearElasticMaterial_e, "e"); // Macro
    IR_GIVE_FIELD(ir, nu, IFT_IsotropicLinearElasticMaterial_n, "n"); // Macro
    IR_GIVE_FIELD(ir, value, IFT_IsotropicLinearElasticMaterial_talpha, "talpha"); // Macro
    propertyDictionary->add(tAlpha, value);
    // compute  value of shear modulus
    G = E / ( 2.0 * ( 1. + nu ) );

    return IRRT_OK;
}


int
IsotropicLinearElasticMaterial :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    LinearElasticMaterial :: giveInputRecordString(str, keyword);
    sprintf( buff, " e %e n %e talpha %e", this->E, this->nu, propertyDictionary->at(tAlpha) );
    str += buff;

    return 1;
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
                                                                MatResponseForm form,
                                                                MatResponseMode mode,
                                                                GaussPoint *gp,
                                                                TimeStep *atTime)
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
                                                           MatResponseForm form,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp,
                                                           TimeStep *atTime)
{
    double e, nu, ee, shear;

    e     = this->E;
    nu    = this->nu;
    ee    = e / ( 1. - nu * nu );
    shear = this->G;

    if ( form == FullForm ) {
        answer.resize(6, 6);
        answer.zero();

        answer.at(1, 1) = ee;
        answer.at(1, 2) = nu * ee;
        answer.at(2, 1) = nu * ee;
        answer.at(2, 2) = ee;
        answer.at(6, 6) = shear;
    } else {
        answer.resize(3, 3);
        answer.zero();

        answer.at(1, 1) = ee;
        answer.at(1, 2) = nu * ee;
        answer.at(2, 1) = nu * ee;
        answer.at(2, 2) = ee;
        answer.at(3, 3) = shear;
    }
}


void
IsotropicLinearElasticMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                                           MatResponseForm form,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp,
                                                           TimeStep *atTime)
{
    double e, nu, ee, shear;

    e     = this->E;
    nu    = this->nu;
    ee    = e / ( 1.0 + nu ) / ( 1. - 2.0 * nu );
    shear = this->G;

    if ( form == FullForm ) {
        answer.resize(6, 6);
        answer.zero();

        answer.at(1, 1) = ee * ( 1.0 - nu );
        answer.at(1, 2) = nu * ee;
        answer.at(1, 3) = nu * ee;
        answer.at(2, 1) = nu * ee;
        answer.at(2, 2) = ee * ( 1.0 - nu );
        answer.at(2, 3) = nu * ee;
        answer.at(3, 1) = nu * ee;
        answer.at(3, 2) = nu * ee;
        answer.at(3, 3) = ee * ( 1.0 - nu );
        answer.at(6, 6) = shear;
    } else {
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
}


void
IsotropicLinearElasticMaterial :: give1dStressStiffMtrx(FloatMatrix &answer,
                                                        MatResponseForm form,
                                                        MatResponseMode mode,
                                                        GaussPoint *gp,
                                                        TimeStep *atTime)
{
    double e;

    e     = this->E;

    if ( form == FullForm ) {
        answer.resize(6, 6);
        answer.zero();

        answer.at(1, 1) = e;
    } else {
        answer.resize(1, 1);
        answer.zero();

        answer.at(1, 1) = e;
    }
}


void
IsotropicLinearElasticMaterial :: give2dBeamStiffMtrx(FloatMatrix &answer,
                                                      MatResponseForm form,
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
    double area, Iy, shearAreaz;


    if ( mode != _2dBeam ) {
        _error("Give2dBeamStiffMtrx : unsupported mode");
    }

    if ( crossSection == NULL ) {
        _error(" Give2dBeamStiffMtrx : no SimpleCrossSection");
    }

    this->give1dStressStiffMtrx(mat3d, FullForm, rMode, gp, tStep);
    area = crossSection->give(CS_Area);
    Iy   = crossSection->give(CS_InertiaMomentY);
    shearAreaz = crossSection->give(CS_SHEAR_AREA_Z);

    if ( form == ReducedForm ) {
        answer.resize(3, 3);
        answer.zero();

        answer.at(1, 1) = mat3d.at(1, 1) * area;
        answer.at(2, 2) = mat3d.at(1, 1) * Iy;
        answer.at(3, 3) = shearAreaz * mat3d.at(1, 1) / ( 2. * ( 1 + nu ) );
    } else {
        answer.resize(8, 8);
        answer.zero();

        answer.at(1, 1) = mat3d.at(1, 1) * area;
        answer.at(5, 5) = mat3d.at(1, 1) * Iy;
        answer.at(7, 7) = shearAreaz * mat3d.at(1, 1) / ( 2. * ( 1 + nu ) );
    }
}


void
IsotropicLinearElasticMaterial :: give3dBeamStiffMtrx(FloatMatrix &answer,
                                                      MatResponseForm form,
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
    double area, E, Iy, Iz, Ik;
    double shearAreay, shearAreaz;

    if ( mode != _3dBeam ) {
        _error("give3dBeamStiffMtrx : unsupported mode");
    }

    if ( crossSection == NULL ) {
        _error("give3dBeamStiffMtrx : no SimpleCrossSection");
    }

    this->give1dStressStiffMtrx(mat3d, FullForm, rMode, gp, tStep);
    E    = mat3d.at(1, 1);
    area = crossSection->give(CS_Area);
    Iy   = crossSection->give(CS_InertiaMomentY);
    Iz   = crossSection->give(CS_InertiaMomentZ);
    Ik   = crossSection->give(CS_TorsionMomentX);

    //shearCoeff = crossSection->give(CS_BeamShearCoeff);
    shearAreay = crossSection->give(CS_SHEAR_AREA_Y);
    shearAreaz = crossSection->give(CS_SHEAR_AREA_Z);

    answer.resize(6, 6);
    answer.zero();

    answer.at(1, 1) = E * area;
    answer.at(2, 2) = shearAreay * this->give('G',gp);
    answer.at(3, 3) = shearAreaz * this->give('G',gp);
    //answer.at(2, 2) = shearCoeff * this->give('G', gp) * area;
    //answer.at(3, 3) = shearCoeff * this->give('G', gp) * area;
    answer.at(4, 4) = this->give('G', gp) * Ik;
    answer.at(5, 5) = E * Iy;
    answer.at(6, 6) = E * Iz;
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
    answer.resize(6);
    answer.zero();
    answer.at(1) = this->give(tAlpha, gp);
    answer.at(2) = this->give(tAlpha, gp);
    answer.at(3) = this->give(tAlpha, gp);
}


MaterialStatus *
IsotropicLinearElasticMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(1, this->giveDomain(), gp);
}
} // end namespace oofem
