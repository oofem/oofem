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

#include "nlisomoisturemat.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {

REGISTER_Material( NlIsoMoistureMaterial );

IRResultType
NlIsoMoistureMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    int type = 0;
    IR_GIVE_FIELD ( ir, type, _IFT_NlIsoMoistureMaterial_isothermtype );

    if ( type >= 6 ) {
        _error ( "initializeFrom: isothermType must be equal to 0, 1, 2 ... 5" );
    }

    this->Isotherm = ( isothermType ) type;

    if ( this->Isotherm == linear ) { // linear isotherm = type 0
        IR_GIVE_FIELD ( ir, moistureCapacity, _IFT_NlIsoMoistureMaterial_capa );
    } else if ( this->Isotherm == multilinear ) { // multilinear isotherm = type 1
        IR_GIVE_FIELD ( ir, iso_h, _IFT_NlIsoMoistureMaterial_iso_h );
        IR_GIVE_FIELD ( ir, iso_wh, _IFT_NlIsoMoistureMaterial_iso_wh );

        if ( ! ( iso_h.giveSize() == iso_wh.giveSize() ) )
            _error ( "initializeFrom: the size of 'iso_h' and 'iso_w(h)' must be the same" );

        for ( int i = 1; i < iso_h.giveSize(); i++ ) {
            if ( ( iso_h.at ( i ) < 0. ) || ( iso_h.at ( i ) > 1. ) )
                _error ( "initializeFrom: iso_h must be in the range <0; 1>" );
        }
    } else if ( this->Isotherm == Ricken ) { // reference mentioned in Kuenzel
        IR_GIVE_FIELD ( ir, dd, _IFT_NlIsoMoistureMaterial_dd );
    } else if ( this->Isotherm == Kuenzel ) {
        IR_GIVE_FIELD ( ir, wf, _IFT_NlIsoMoistureMaterial_wf );
        IR_GIVE_FIELD ( ir, b, _IFT_NlIsoMoistureMaterial_b );
    } else if ( this->Isotherm == Hansen ) {
        IR_GIVE_FIELD ( ir, uh, _IFT_NlIsoMoistureMaterial_uh );
        IR_GIVE_FIELD ( ir, A, _IFT_NlIsoMoistureMaterial_a );
        IR_GIVE_FIELD ( ir, nn, _IFT_NlIsoMoistureMaterial_nn );
        IR_GIVE_FIELD ( ir, rhodry, _IFT_NlIsoMoistureMaterial_rhodry );
    } else if ( this->Isotherm == BSB ) {
        IR_GIVE_FIELD ( ir, c, _IFT_NlIsoMoistureMaterial_c );
        IR_GIVE_FIELD ( ir, k, _IFT_NlIsoMoistureMaterial_k );
        IR_GIVE_FIELD ( ir, Vm, _IFT_NlIsoMoistureMaterial_vm );
        IR_GIVE_FIELD ( ir, rhodry, _IFT_NlIsoMoistureMaterial_rhodry );
    } else {
        _error ( "initializeFrom: isothermType must be equal to 0, 1, 2 ... 5" );
    }

    IR_GIVE_FIELD ( ir, type, _IFT_NlIsoMoistureMaterial_permeabilitytype );

    if ( type >= 3 ) {
        _error ( "initializeFrom: isothermType must be equal to 0, 1, 2" );
    }

    this->Permeability = ( permeabilityType ) type;

    if ( this->Permeability == multilin ) {

        IR_GIVE_FIELD ( ir, perm_h, _IFT_NlIsoMoistureMaterial_perm_h );
        IR_GIVE_FIELD ( ir, perm_ch, _IFT_NlIsoMoistureMaterial_perm_ch );

        if ( ! ( perm_h.giveSize() == perm_ch.giveSize() ) )
            _error ( "initializeFrom: the size of 'perm_h' and 'perm_c(h)' must be the same" );

        for ( int i = 1; i < perm_h.giveSize(); i++ ) {
            if ( ( perm_h.at ( i ) < 0. ) || ( perm_h.at ( i ) > 1. ) )
                _error ( "initializeFrom: perm_h must be in the range <0; 1>" );
        }

    } else if ( this->Permeability == Bazant ) {
        IR_GIVE_FIELD ( ir, C1, _IFT_NlIsoMoistureMaterial_c1 );
        IR_GIVE_FIELD ( ir, n, _IFT_NlIsoMoistureMaterial_n );
        IR_GIVE_FIELD ( ir, alpha0, _IFT_NlIsoMoistureMaterial_alpha0 );
        IR_GIVE_FIELD ( ir, hC, _IFT_NlIsoMoistureMaterial_hc );
    } else if ( this->Permeability == Xi ) {
        IR_GIVE_FIELD ( ir, alphah, _IFT_NlIsoMoistureMaterial_alphah );
        IR_GIVE_FIELD ( ir, betah, _IFT_NlIsoMoistureMaterial_betah );
        IR_GIVE_FIELD ( ir, gammah, _IFT_NlIsoMoistureMaterial_gammah );
    } else {
        _error ( "initializeFrom: permeabilityType must be equal to 0, 1 or 2" );
    }

    IsotropicMoistureTransferMaterial :: initializeFrom ( ir );
    return IRRT_OK;
}

double
NlIsoMoistureMaterial :: giveMoistureCapacity(GaussPoint *gp, TimeStep *atTime)
{
    double humidity = this->giveHumidity(gp, VM_Total);

    if ( this->Isotherm == linear ) {
        return moistureCapacity;

    } else if ( this->Isotherm == multilinear ) {
        double tol = 1.e-10;
        for (int i = 1; i <= iso_h.giveSize(); i++) {
        if ( ( humidity - iso_h.at(i) ) < tol ) {
            return ((iso_wh.at(i)-iso_wh.at(i-1))/(iso_h.at(i)-iso_h.at(i-1)));
        }
    }

    } else if ( this->Isotherm == Ricken ) {
        return 1./(dd*(1.-humidity));

    } else if ( this->Isotherm == Kuenzel ) {
        return wf*(b-1.)*b / ((b-humidity)*(b-humidity));

    } else if ( this->Isotherm == Hansen ) {
        return rhodry*uh/(A*nn*humidity*pow((1-log(humidity)/A),(1/nn + 1)));

    } else if ( this->Isotherm == BSB ) {
        double nominator, denominator;
        nominator = c*k*Vm*rhodry*(1. + k*k*humidity*humidity*c - k*k*humidity*humidity);
        denominator = (1.-k*humidity)*(1.-k*humidity) * (1.+(c-1.)*k*humidity)*(1.+(c-1.)*k*humidity);
        return nominator/denominator;

    } else {
        _error("initializeFrom: isothermType must be equal to 0, 1, 2 ... 5");
    }

    return 0.;
}

double
NlIsoMoistureMaterial :: givePermeability(GaussPoint *gp, TimeStep *atTime)
{
    double permeability = 0.;
    double humidity = this->giveHumidity(gp, VM_Total);

    if ( this->Permeability == multilin ) {
        double tol = 1.e-10;
        for (int i = 1; i <= perm_h.giveSize(); i++) {
            if ((humidity - perm_h.at(i))<tol) {
                permeability = perm_ch.at(i-1) + (perm_ch.at(i)-perm_ch.at(i-1))*(humidity-perm_h.at(i-1))/(perm_h.at(i)-perm_h.at(i-1));
                break;
            }
        }

    } else if ( this->Permeability == Bazant ){
        permeability = C1 * (alpha0 + (1.-alpha0)/(1. + pow( (1.-humidity)/(1.-hC), n) ) ) ;

    } else if ( this->Permeability == Xi ){
        double power = pow(10., gammah*(humidity-1.));
        permeability = alphah + betah * (1. - pow(2., -power));

    } else {
        _error("initializeFrom: permeabilityType must be equal to 0, 1 or 2");
    }

    return permeability;
}

double
NlIsoMoistureMaterial :: giveHumidity(GaussPoint *gp, ValueModeType mode)
{
    const FloatArray &tempState = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) )->giveTempStateVector();
    if ((tempState.at(1) > 1.0)||(tempState.at(1) < 0.0)) {
        OOFEM_ERROR2("NlIsoMoistureMaterial :: giveHumidity : Relative humidity %.3f is out of range", tempState.at(1));
        return 0.0;
    } else {
        return tempState.at(1);
    }
}


} // end namespace oofem
