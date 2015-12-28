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

#include "hemobaznajmat.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {

REGISTER_Material( HeMoBazNajMaterial );

int
HeMoBazNajMaterial :: hasMaterialModeCapability(MaterialMode mode)
{
    if ( ( mode == _2dHeMo ) || ( mode == _3dHeMo ) ) {
        return 1;
    }

    return 0;
}


IRResultType
HeMoBazNajMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, C1, _IFT_HeMoBazNajMaterial_c1);
    IR_GIVE_FIELD(ir, n, _IFT_HeMoBazNajMaterial_n);
    IR_GIVE_FIELD(ir, alpha0, _IFT_HeMoBazNajMaterial_alpha0);
    IR_GIVE_FIELD(ir, hC, _IFT_HeMoBazNajMaterial_hc);

    this->moistureCapacity = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, moistureCapacity, _IFT_HeMoBazNajMaterial_capa);

    IR_GIVE_FIELD(ir, heatConductivity, _IFT_HeMoBazNajMaterial_k);
    IR_GIVE_FIELD(ir, heatCapacity, _IFT_HeMoBazNajMaterial_c);

     return Material :: initializeFrom(ir);
}


double
HeMoBazNajMaterial :: give(int aProperty, GaussPoint *gp)
//
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
//
{
    return this->Material :: give(aProperty, gp);
}


void
HeMoBazNajMaterial :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep)
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );

    FloatArray s;
    s = ms->giveTempField();
    if ( s.isEmpty() ) {
        OOFEM_ERROR("undefined state vector");
    }
    double h = s.at(2);
    double t = s.at(1);

    FloatArray ans_w, ans_t;
    FloatArray grad_w, grad_t;
    int size = grad.giveSize() / 2;
    for (int i = 1; i <= size; ++i) {
        grad_w.at(i) = grad.at(i);
    }
    for (int i = size+1; i <= size*2; ++i) {
        grad_t.at(i) = grad.at(i);
    }

    ans_w.beScaled(perm_mm(h, t), grad_w);
    ans_w.beScaled(perm_mh(h, t), grad_t);
    ans_t.beScaled(perm_hm(h, t), grad_w);
    ans_t.beScaled(perm_hh(h, t), grad_t);
    
    answer.resize(size * 2);
    answer.zero();
    answer.addSubVector(ans_w, 1);
    answer.addSubVector(ans_t, size+1);

    ms->setTempField(field);
    ms->setTempGradient(grad);
    ms->setTempFlux(answer);
}


void
HeMoBazNajMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime)
{
    /*
     * returns constitutive matrix of receiver
     */
    if ( ( mode == Conductivity_ww ) || ( mode == Conductivity_hh ) || ( mode == Conductivity_hw ) || ( mode == Conductivity_wh ) ) {
        this->computeConductivityMtrx(answer, mode, gp, atTime);
    } else {
        OOFEM_ERROR( "giveCharacteristicMatrix : unknown mode (%s)", __MatResponseModeToString(mode) );
    }
}


double
HeMoBazNajMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *atTime)
{
    return this->computeCapacityCoeff(mode, gp, atTime);
}


void HeMoBazNajMaterial :: computeConductivityMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    MaterialMode mmode = gp->giveMaterialMode();
    switch ( mmode ) {
    case _2dHeMo:
        this->matcond2d(answer, gp, mode, atTime);
        return;

    case _3dHeMo:
        this->matcond3d(answer, gp, mode, atTime);
        return;

    default:
        OOFEM_ERROR("Unsupported MaterialMode");
    }
}


void
HeMoBazNajMaterial :: matcond1d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime)
{
    double k = 0.0, h = 0.0, t = 0.0;
    TransportMaterialStatus *status = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    
    FloatArray s;
    s = status->giveTempField();
    if ( s.isEmpty() ) {
        OOFEM_ERROR("undefined state vector");
    }
    h = s.at(2);
    t = s.at(1);


    if ( mode == Conductivity_ww ) {
        k = perm_mm(h, t);
    } else if ( mode == Conductivity_wh ) {
        k = perm_mh(h, t);
    } else if ( mode == Conductivity_hw ) {
        k = perm_hm(h, t);
    } else if ( mode == Conductivity_hh ) {
        k = perm_hh(h, t);
    } else {
        OOFEM_ERROR("Unknown MatResponseMode");
    }

    d.resize(1, 1);
    d.at(1, 1) = k;
}

void
HeMoBazNajMaterial :: matcond2d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime)
{
    double k = 0.0, h = 0.0, t = 0.0;
    TransportMaterialStatus *status = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );

    FloatArray s;
    s = status->giveTempField();
    if ( s.isEmpty() ) {
        OOFEM_ERROR("undefined state vector");
    }
    h = s.at(2);
    t = s.at(1);
    
    if ( mode == Conductivity_ww ) {
        k = perm_mm(h, t);
    } else if ( mode == Conductivity_wh ) {
        k = perm_mh(h, t);
    } else if ( mode == Conductivity_hw ) {
        k = perm_hm(h, t);
    } else if ( mode == Conductivity_hh ) {
        k = perm_hh(h, t);
    } else {
        OOFEM_ERROR("Unknown MatResponseMode");
    }

    d.resize(2, 2);
    d.at(1, 1) = k;
    d.at(1, 2) = 0.0;
    d.at(2, 1) = 0.0;
    d.at(2, 2) = k;
}

void
HeMoBazNajMaterial :: matcond3d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime)
{
    double k = 0.0, h = 0.0, t = 0.0;
    TransportMaterialStatus *status = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );

    FloatArray s;
    s = status->giveTempField();
    if ( s.isEmpty() ) {
        OOFEM_ERROR("undefined state vector");
    }
    h = s.at(2);
    t = s.at(1);

    if ( mode == Conductivity_ww ) {
        k = perm_mm(h, t);
    } else if ( mode == Conductivity_wh ) {
        k = perm_mh(h, t);
    } else if ( mode == Conductivity_hw ) {
        k = perm_hm(h, t);
    } else if ( mode == Conductivity_hh ) {
        k = perm_hh(h, t);
    } else {
        OOFEM_ERROR("Unknown MatResponseMode");
    }

    d.resize(3, 3);
    d.at(1, 1) = k;
    d.at(1, 2) = 0.0;
    d.at(1, 3) = 0.0;
    d.at(2, 1) = 0.0;
    d.at(2, 2) = k;
    d.at(2, 3) = 0.0;
    d.at(3, 1) = 0.0;
    d.at(3, 2) = 0.0;
    d.at(3, 3) = k;
}


double HeMoBazNajMaterial :: computeCapacityCoeff(MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    if ( mode == Capacity_ww ) {
      return this->moistureCapacity;

    } else if ( mode == Capacity_wh ) {
      return 0.0;

    } else if ( mode == Capacity_hw ) {
      return 0.0;

    } else if ( mode == Capacity_hh ) {
      return this->heatCapacity;

    } else {
        OOFEM_ERROR("Unknown MatResponseMode");
    }

    return 0.0; // to make compiler happy
}


double
HeMoBazNajMaterial :: perm_mm(double h, double T)
{
    return ( C1 * ( alpha0 + ( 1. - alpha0 ) / ( 1. + pow( ( 1. - h ) / ( 1. - hC ), n ) ) ) );
}

double
HeMoBazNajMaterial :: perm_mh(double h, double T)
{
    return (0.);
}

double
HeMoBazNajMaterial :: perm_hm(double h, double T)
{
    return (0.);
}

double
HeMoBazNajMaterial :: perm_hh(double h, double T)
{
  return this->heatConductivity;
}

bool
HeMoBazNajMaterial :: isCharacteristicMtrxSymmetric(MatResponseMode mode)
{
    if ( ( mode == Conductivity_ww ) || ( mode == Conductivity_hh ) || ( mode == Conductivity_hw ) || ( mode == Conductivity_wh ) ) {
        return true;
    } else {
        OOFEM_ERROR( "isCharacteristicMtrxSymmetric : unknown mode (%s)", __MatResponseModeToString(mode) );
    }

    return false; // to make compiler happy
}

int
HeMoBazNajMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime)
// IST_Humidity overriden to use inverse_sorption_isotherm
{
    if ( type == IST_Humidity ) {
        answer.resize(1);
        answer.at(1) = giveHumidity(gp, VM_Velocity); // VM_Previous = equilibrated value of humidity
        return 1;
    } else {
        return TransportMaterial :: giveIPValue(answer, gp, type, atTime);
    }
}
} // end namespace oofem
