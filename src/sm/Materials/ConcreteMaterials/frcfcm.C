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
 *               Copyright (C) 1993 - 2015   Borek Patzak
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

#include "frcfcm.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"
#include "contextioerr.h"
#include "datastream.h"

namespace oofem {
REGISTER_Material(FRCFCM);

FRCFCM :: FRCFCM(int n, Domain *d) : ConcreteFCM(n, d)
{
    fiberShearStrengthType = FSS_Unknown;
}



IRResultType
FRCFCM :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = ConcreteFCM :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }


    this->fibreActivationOpening = 0.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, fibreActivationOpening, _IFT_FRCFCM_fibreActivationOpening);
    if ( this->fibreActivationOpening < 0. ) {
        OOFEM_ERROR("FibreActivationOpening must be positive");
    }

    this->dw0 = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, dw0, _IFT_FRCFCM_dw0);
    if ( ( this->dw0 < 0. ) || ( this->dw0 > this->fibreActivationOpening ) ) {
        OOFEM_ERROR("dw0 must be positive and smaller than fibreActivationOpening");
    }

    this->dw1 = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, dw1, _IFT_FRCFCM_dw1);
    if ( ( this->dw1 < 0. ) ) {
        OOFEM_ERROR("dw1 must be positive");
    }

    this->smoothen = false;
    if ( ( ( this->dw1 == 0. ) && ( this->dw0 > 0. ) ) || ( ( this->dw0 == 0. ) && ( this->dw1 > 0. ) ) ) {
        OOFEM_ERROR("both dw0 and dw1 must be specified at the same time");
    } else {
        this->smoothen = true;
        this->dw0 *= -1.;
    }

    // type of SHEAR STRESS FOR FIBER PULLOUT
    int type = 0;

    IR_GIVE_OPTIONAL_FIELD(ir, type, _IFT_FRCFCM_fssType);
    switch ( type ) {
    case 0:
        fiberShearStrengthType = FSS_NONE;
        break;

    case 1:
        fiberShearStrengthType = FSS_Sajdlova;
        IR_GIVE_FIELD(ir, b0, _IFT_FRCFCM_b0);
        break;

    case 2:
        fiberShearStrengthType = FSS_Kabele;
        IR_GIVE_FIELD(ir, b1, _IFT_FRCFCM_b1);
        IR_GIVE_FIELD(ir, b2, _IFT_FRCFCM_b2);
        IR_GIVE_FIELD(ir, b3, _IFT_FRCFCM_b3);
        break;

    case 3:
        fiberShearStrengthType = FSS_Havlasek;
        IR_GIVE_FIELD(ir, b1, _IFT_FRCFCM_b1);
        IR_GIVE_FIELD(ir, b2, _IFT_FRCFCM_b2);
        IR_GIVE_FIELD(ir, b3, _IFT_FRCFCM_b3);
        break;

    default:
        fiberShearStrengthType = FSS_Unknown;
        OOFEM_WARNING("Softening type number %d is unknown", type);
        return IRRT_BAD_FORMAT;
    }



    // type of FIBER DAMAGE
    type = 0;

    IR_GIVE_OPTIONAL_FIELD(ir, type, _IFT_FRCFCM_fDamType);
    switch ( type ) {
    case 0:
        fiberDamageType = FDAM_NONE;
        break;

    case 1:
        fiberDamageType = FDAM_GammaCrackLin;
        IR_GIVE_FIELD(ir, gammaCrackFail, _IFT_FRCFCM_gammaCrack);
        break;

    case 2:
        fiberDamageType = FDAM_GammaCrackExp;
        IR_GIVE_FIELD(ir, gammaCrackFail, _IFT_FRCFCM_gammaCrack);
        break;

    default:
        fiberDamageType = FDAM_Unknown;
        OOFEM_WARNING("Fibre damage type number %d is unknown", type);
        return IRRT_BAD_FORMAT;
    }


    IR_GIVE_OPTIONAL_FIELD(ir, type, _IFT_FRCFCM_fiberType);
    switch ( type ) {
    case 0:
        fiberType = FT_CAF;
        break;

    case 1:
        fiberType = FT_SAF;
        break;

    case 2:
        fiberType = FT_SRF;
        break;

    default:
        fiberType = FT_Unknown;
        OOFEM_WARNING("Fibre type number %d is unknown", type);
        return IRRT_BAD_FORMAT;
    }

    if ( ( fiberType == FT_CAF ) || ( fiberType == FT_SAF ) ) {
        if  ( ir->hasField(_IFT_FRCFCM_orientationVector) ) {
            IR_GIVE_FIELD(ir, orientationVector, _IFT_FRCFCM_orientationVector);

            if ( !( ( this->orientationVector.giveSize() == 2 ) || ( this->orientationVector.giveSize() == 3 ) ) ) {
                OOFEM_ERROR("length of the fibre orientation vector must be 2 for 2D and 3 for 3D analysis");
            }

            // we normalize the user-defined orientation vector
            double length = 0.;
            for ( int i = 1; i <= this->orientationVector.giveSize(); i++ ) {
                length += pow(this->orientationVector.at(i), 2);
            }
            length = sqrt(length);

            this->orientationVector.times(1 / length);
        } else {
            this->orientationVector.resize(3);
            this->orientationVector.zero();
            this->orientationVector.at(1) = 1.;
        }
    }

    // general properties
    IR_GIVE_FIELD(ir, Vf, _IFT_FRCFCM_Vf);
    if ( Vf < 0. ) {
        OOFEM_ERROR("fibre volume content must not be negative");
    }

    IR_GIVE_FIELD(ir, Df, _IFT_FRCFCM_Df);
    if ( Df <= 0. ) {
        OOFEM_ERROR("fibre diameter must be positive");
    }

    IR_GIVE_FIELD(ir, Ef, _IFT_FRCFCM_Ef);
    if ( Ef <= 0. ) {
        OOFEM_ERROR("fibre stiffness must be positive");
    }

    // compute or read shear modulus of fibers
    double nuf = 0.;
    Gfib = 0.;
    if  ( ir->hasField(_IFT_FRCFCM_nuf) ) {
        IR_GIVE_FIELD(ir, nuf, _IFT_FRCFCM_nuf);
        Gfib = Ef / ( 2. * ( 1. + nuf ) );
    } else {
        IR_GIVE_FIELD(ir, Gfib, _IFT_FRCFCM_Gfib);
    }

    this->kfib = 0.9;
    IR_GIVE_OPTIONAL_FIELD(ir, kfib, _IFT_FRCFCM_kfib);

    // debonding
    IR_GIVE_FIELD(ir, tau_0, _IFT_FRCFCM_tau_0);
    if ( tau_0 <= 0. ) {
        OOFEM_ERROR("shear strength must be positive");
    }

    // snubbing coefficient
    IR_GIVE_FIELD(ir, f, _IFT_FRCFCM_f);
    if ( f < 0. ) {
        OOFEM_ERROR("snubbing coefficient must not be negative");
    }

    // for SRF only
    if ( fiberType == FT_SRF ) {
        this->g = 2. * ( 1. + exp(M_PI * f / 2.) ) / ( 4. + f * f );
    }

    double Em;
    IR_GIVE_FIELD(ir, Em, _IFT_IsotropicLinearElasticMaterial_e);

    this->eta = this->Ef * this->Vf / ( Em * ( 1. - this->Vf ) );

    this->M = 4; // exponent in unloading-reloading law
    IR_GIVE_OPTIONAL_FIELD(ir, M, _IFT_FRCFCM_M);


    if ( ( fiberType == FT_SAF ) || ( fiberType == FT_SRF ) ) {
        IR_GIVE_FIELD(ir, Lf, _IFT_FRCFCM_Lf);
        // transitional opening at which sigma_b = max (zero derivative) for SRF and SAF
        this->w_star = this->Lf * this->Lf * this->tau_0 / ( ( 1. + this->eta ) * this->Ef * this->Df );
    }

    if  ( ir->hasField(_IFT_FRCFCM_computeCrackSpacing) ) {
        this->crackSpacing = this->computeCrackSpacing();
    }

    return IRRT_OK;
}



double
FRCFCM :: computeCrackFibreAngle(GaussPoint *gp, int index)
{
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    double theta = 0.;

    for ( int i = 1; i <= min( status->giveCrackDirs().giveNumberOfRows(), this->orientationVector.giveSize() ); i++ ) {
        theta += status->giveCrackDirs().at(i, index) * this->orientationVector.at(i);
    }

    // can be exceeded due to truncation error
    if ( theta > 1 ) {
        return 0.;
    } else if ( theta < -1. ) {
        return 0.;
    }

    // transform into angle
    theta = acos(theta);


    if ( theta > M_PI / 2. ) {
        theta = M_PI - theta;
    }


    return theta;
}



double
FRCFCM :: giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp, int i)
//
// returns current cracking modulus according to crackStrain for i-th
// crackplane
//
{
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    double Cfc; //stiffness of cracked concrete and fibers
    double Cff = 0.;
    double theta = 0.;
    double ec, emax, Le, Ncr, w, w_max, factor, omega;
    double max_traction; // for secant stiffness
    double sig_dw1, Dsig_dw1, x, dw_smooth, C1, C2; // smoothen


    // adjust concrete modulus by fibre volume fraction and by the crack orientation wrt fibres
    // it turns out that the concrete area is in the end independent of the angle
    // larger angle causes less fiber to cross the crack but also larger area of the fiber section area (an ellipse)
    // these two effects cancel out and what only matters is the fiber volume Vf
    Cfc = ConcreteFCM :: giveCrackingModulus(rMode, gp, i);
    Cfc *= ( 1. - this->Vf );


    if ( ( this->fiberType == FT_CAF ) || ( this->fiberType == FT_SAF ) ) {
        theta = fabs( this->computeCrackFibreAngle(gp, i) );
    }

    // in case of continuous aligned fibres and short aligned fibers the TANGENT stiffness of the fibers is multiplied by cos(theta) to reflect decreased number of crossing fibers and by exp(theta * f) to consider the snubbing effect. The secant stiffness uses a function for computing the normal stress where "theta" is used too.

    // virgin material -> return stiffness of concrete
    if ( status->giveTempCrackStatus(i) == pscm_NONE ) {
        return Cfc;
    }

    // modification to take into account crack spacing - makes sense only in crack-opening-based materials

    Ncr = this->giveNumberOfCracksInDirection(gp, i);

    ec =  status->giveTempCrackStrain(i) / Ncr;
    emax = max(status->giveMaxCrackStrain(i) / Ncr, ec);

    Le = status->giveCharLength(i);

    w_max = emax * Le; //maximal crack opening
    w_max -= fibreActivationOpening;

    w = ec * Le; //current crack opening
    w -= fibreActivationOpening;


    if ( rMode == TangentStiffness ) { // assumption of constant tau_s(w) in the derivative
        // for zero or negative strain has been taken care of before

        if ( this->smoothen ) {
            if ( ( w_max < this->dw0 ) || ( w < this->dw0 ) ) {
                return Cfc;
            }
        } else {
            if ( ( w_max < 0. ) || ( w < 0. ) ) {
                return Cfc;
            }
        }


        if ( ( this->smoothen ) && ( w_max > this->dw0 ) && ( w_max < this->dw1 ) ) {
            if ( w == w_max ) {
                omega = this->computeTempDamage(gp);

                x = w_max - this->dw0;
                dw_smooth = this->dw1 - this->dw0;

                if ( fiberType == FT_CAF ) { // continuous aligned fibers
                    sig_dw1 =  2. *this->Vf *sqrt(this->dw1 *this->Ef * ( 1. + this->eta ) *this->tau_0 / this->Df); // stress in fibers per unit area of concrete
                    Dsig_dw1 = this->Vf * sqrt( this->Ef * ( 1. + this->eta ) * this->tau_0 / ( this->Df * this->dw1 ) ); // stress derivative

                    C1 = Dsig_dw1 / pow(dw_smooth, 2) - ( 2 * sig_dw1 ) / pow(dw_smooth, 3);
                    C2 = ( 3 * sig_dw1 ) / pow(dw_smooth, 2) - Dsig_dw1 / dw_smooth;

                    Cff = 3. *C1 *pow(x, 2) + 2. * C2 * x;
                    Cff *= Le;

                    // reflect fibre orientation wrt crack
                    theta = fabs( this->computeCrackFibreAngle(gp, i) );
                    Cff *= fabs( cos(theta) ) * exp(theta * this->f);

                    Cff *= ( 1. - omega );
                } else if ( fiberType == FT_SAF ) { // short aligned fibers
                    sig_dw1 =  this->Vf * ( sqrt(4. * this->Ef * ( 1. + this->eta ) * this->dw1 * this->tau_0 / this->Df) - this->Ef * ( 1. + this->eta ) * this->dw1 / this->Lf );
                    Dsig_dw1 = this->Vf * ( sqrt( this->Ef * ( 1. + this->eta ) * this->tau_0 / ( this->Df * this->dw1 ) ) - this->Ef * ( 1. + this->eta ) / this->Lf );

                    C1 = Dsig_dw1 / pow(dw_smooth, 2) - ( 2 * sig_dw1 ) / pow(dw_smooth, 3);
                    C2 = ( 3 * sig_dw1 ) / pow(dw_smooth, 2) - Dsig_dw1 / dw_smooth;

                    Cff = 3. *C1 *pow(x, 2) + 2. * C2 * x;
                    Cff *= Le;

                    // reflect fibre orientation wrt crack
                    theta = fabs( this->computeCrackFibreAngle(gp, i) );
                    Cff *= fabs( cos(theta) ) * exp(theta * this->f);

                    Cff *= ( 1. - omega );
                } else if ( fiberType == FT_SRF ) { // short random fibers
                    factor = ( 1. - omega ) * this->g * this->Vf * this->Lf / ( 2. * this->Df );

                    sig_dw1 = factor * this->tau_0 * ( 2. * sqrt(this->dw1 / this->w_star) - this->dw1 / this->w_star );
                    Dsig_dw1 = factor * this->tau_0 * ( 1. / ( this->w_star * sqrt(this->dw1 / this->w_star) ) - 1. / this->w_star );

                    C1 = Dsig_dw1 / pow(dw_smooth, 2) - ( 2 * sig_dw1 ) / pow(dw_smooth, 3);
                    C2 = ( 3 * sig_dw1 ) / pow(dw_smooth, 2) - Dsig_dw1 / dw_smooth;

                    Cff = 3. *C1 *pow(x, 2) + 2. * C2 * x;
                    Cff *= Le;
                }
            } else { // unlo-relo with smoothing
                // compute traction for ec == emax
                double max_traction = this->computeStressInFibersInCracked(gp, emax * Ncr, i);
                Cff = ( this->M * max_traction * Le / ( w_max - this->dw0 ) ) * pow( ( w - this->dw0 ) / ( w_max - this->dw0 ), ( this->M - 1 ) );
            }
        } else if ( w_max == 0. ) { //zero strain
            w_max = 1.e-10;

            if ( fiberType == FT_CAF ) { // continuous aligned fibers
                Cff = this->Vf * Le * sqrt( this->Ef * ( 1. + this->eta ) * this->tau_0 / ( this->Df * w_max ) );

                // reflect fibre orientation wrt crack
                Cff *= fabs( cos(theta) ) * exp(theta * this->f);
            } else if ( fiberType == FT_SAF ) { // short aligned fibers
                Cff = this->Vf * Le * ( sqrt( this->Ef * ( 1. + this->eta ) * this->tau_0 / ( this->Df * w_max ) ) - this->Ef * ( 1. + this->eta ) / this->Lf );

                // reflect fibre orientation wrt crack
                Cff *= fabs( cos(theta) ) * exp(theta * this->f);
            } else if ( fiberType == FT_SRF ) { // short random fibers
                factor = this->g * this->Vf * this->Lf / ( 2. * this->Df );
                Cff = factor * this->tau_0 * ( Le / ( this->w_star * sqrt(w_max / this->w_star) ) - Le / this->w_star );
            } else {
                OOFEM_ERROR("Unknown fiber type");
            }
        } else if ( w == w_max ) { // softening
            if ( fiberType == FT_CAF ) { // continuous aligned fibers
                Cff = this->Vf * Le * sqrt( this->Ef * ( 1. + this->eta ) * this->tau_0 / ( this->Df * w_max ) );
                // reflect fibre orientation wrt crack
                Cff *= fabs( cos(theta) ) * exp(theta * this->f);

                omega = this->computeTempDamage(gp);
                Cff *= ( 1. - omega );
            } else if ( fiberType == FT_SAF ) { // short aligned fibers
                omega = this->computeTempDamage(gp);


                if ( w_max < this->w_star ) { // debonding + pullout
                    Cff = ( 1. - omega ) * this->Vf * Le * ( sqrt( this->Ef * ( 1. + this->eta ) * this->tau_0 / ( this->Df * w_max ) ) - this->Ef * ( 1. + this->eta ) / this->Lf );

                    // reflect fibre orientation wrt crack
                    Cff *= fabs( cos(theta) ) * exp(theta * this->f);
                } else if ( w_max <= this->Lf / 2. ) { // pullout
                    factor = ( 1. - omega ) * this->Vf * this->Lf / this->Df;

                    Cff =  factor * this->computeFiberBond(w_max) * ( 8. * w_max * Le / ( this->Lf * this->Lf ) - 4. * Le / this->Lf );

                    // reflect fibre orientation wrt crack
                    Cff *= fabs( cos(theta) ) * exp(theta * this->f);
                } else { // fully pulled out
                    Cff = 0.;
                }
            } else if ( fiberType == FT_SRF ) { // short random fibers
                omega = this->computeTempDamage(gp);
                factor = ( 1. - omega ) * this->g * this->Vf * this->Lf / ( 2. * this->Df );

                if ( w_max < this->w_star ) { // debonding + pullout
                    Cff = factor * this->tau_0 * ( Le / ( this->w_star * sqrt(w_max / this->w_star) ) - Le / this->w_star );
                } else if ( w_max <= this->Lf / 2. ) { // pullout
                    Cff = factor * this->computeFiberBond(w_max) * ( 8. * w_max * Le / ( this->Lf * this->Lf ) - 4. * Le / this->Lf );
                } else { // fully pulled out
                    Cff = 0.;
                }
            } else {
                OOFEM_ERROR("Unknown fiber type");
            }
        } else { // unlo-relo
            // compute traction for ec == emax
            max_traction = this->computeStressInFibersInCracked(gp, emax * Ncr, i);

            Cff = ( this->M * max_traction * Le / w_max ) * pow( ( w / w_max ), ( this->M - 1 ) );
        }
    } else if ( rMode == SecantStiffness ) {
        if ( this->smoothen ) {
            if  ( ( w_max < this->dw0 ) || ( w < this->dw0 ) ) { // fibres are not activated
                Cff = 0.;
            } else if ( ec == emax ) { // softening
                Cff = computeStressInFibersInCracked(gp, emax * Ncr, i); // gets stress
                Cff /= ( emax ); // converted into stiffness
            } else { // unlo-relo
                // compute traction for ec == emax
                max_traction = this->computeStressInFibersInCracked(gp, emax * Ncr, i);
                Cff =  max_traction * pow( ( ( w - this->dw0 ) / ( w_max - this->dw0 ) ), this->M ) / ec;
            }
        } else {
            if ( ( w_max < 0. ) || ( w < 0. ) ) { // fibres are not activated
                Cff = 0.;
            } else if ( w_max == 0. ) { //zero strain
                w_max = 1.e-10 / Ncr;
                emax = w_max / Le;

                Cff = computeStressInFibersInCracked(gp, emax * Ncr, i); // gets stress
                Cff /= ( emax ); // converted into stiffness
            } else if ( ec == emax ) { // softening
                Cff = computeStressInFibersInCracked(gp, emax * Ncr, i); // gets stress
                Cff /= ( emax ); // converted into stiffness
            } else { // unlo-relo
                // compute traction for ec == emax
                max_traction = this->computeStressInFibersInCracked(gp, emax * Ncr, i);

                Cff =  max_traction * pow( ( w / w_max ), this->M ) / ec;
            }
        }
    } else {
        OOFEM_ERROR("Unknown material response mode");
    }

    // to take into account crack-spacing
    Cff /= Ncr;

    return Cfc + Cff;
}



double
FRCFCM :: computeFiberBond(double w)
{
    double tau_s = 0.;
    double dw, tau_tilde = 0.;

    if ( ( w <= 0. ) || ( fiberShearStrengthType == FSS_NONE ) || ( fiberType == FT_CAF ) ) {
        return this->tau_0;
    }

    if  ( fiberShearStrengthType == FSS_Sajdlova ) { // Sajdlova
        tau_s = this->tau_0 * ( 1. + sgn(this->b0) * ( 1. - exp(-fabs(this->b0) * w / this->Df) ) );
    } else if ( fiberShearStrengthType == FSS_Kabele ) { // Kabele
        tau_s = this->tau_0 * ( 1 + this->b1 * ( w / this->Df ) + this->b2 * ( w / this->Df ) * ( w / this->Df ) + this->b3 * ( w / this->Df ) * ( w / this->Df ) * ( w / this->Df ) );
    } else if ( fiberShearStrengthType == FSS_Havlasek ) { // Havlasek
        dw = w - this->w_star;

        if ( fiberType == FT_SRF ) {
            tau_tilde = this->tau_0 / ( ( 1. - this->w_star / this->Lf ) * ( 1. - this->w_star / this->Lf ) );
        } else { // SAF
            tau_tilde = this->tau_0 * this->Ef * ( 1. + this->eta ) * this->Df / ( this->Ef * ( 1. + this->eta ) * this->Df - 2. * this->Lf * this->tau_0 );
        }

        tau_s = tau_tilde + this->tau_0 * ( this->b1 * ( dw / this->Df ) +
                                            this->b2 * ( dw / this->Df ) * ( dw / this->Df ) +
                                            this->b3 * ( dw / this->Df ) * ( dw / this->Df ) * ( dw / this->Df ) );
    } else {
        OOFEM_ERROR("Unknown FiberShearStrengthType");
    }

    return tau_s;
}




double
FRCFCM :: giveNormalCrackingStress(GaussPoint *gp, double ec, int i)
//
// returns receivers Normal Stress in crack i  for given cracking strain
//
{
    double traction_fc, traction_ff; //normal stress in cracked concrete and fibers

    traction_fc =  ConcreteFCM :: giveNormalCrackingStress(gp, ec, i);
    traction_fc *= ( 1. - this->Vf );

    traction_ff = this->computeStressInFibersInCracked(gp, ec, i);


    return traction_fc + traction_ff;
}


double
FRCFCM :: computeStressInFibersInCracked(GaussPoint *gp, double ec, int i)
//
// returns receivers Normal Stress in crack i  for given cracking strain
//
{
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    double traction_ff = 0.; //normal stress in cracked concrete in its FIBERS (continuum point of view)

    double emax, Le, w_max, w, Ncr, omega, factor;
    double theta = 0.;

    double sig_dw1, Dsig_dw1, x, dw_smooth, C1, C2; // smoothen


    if ( status->giveTempCrackStatus(i) == pscm_NONE ) {
        return 0.;
    }

    emax = max(status->giveMaxCrackStrain(i), ec);

    if ( emax == 0. ) {
        return 0.;
    }

    Ncr = this->giveNumberOfCracksInDirection(gp, i);
    emax /= Ncr;
    ec /= Ncr;

    Le = status->giveCharLength(i);

    w_max = emax * Le; // maximal crack opening
    // w_max already corresponds to one parallel
    w_max -= fibreActivationOpening; // offset the initial crack opening when fibers do not extend but matrix has already cracked

    w = ec * Le; // crack opening
    w -= fibreActivationOpening; // offset the initial crack opening when fibers do not extend but matrix has already cracked


    if ( this->smoothen ) {
        if ( ( w_max < this->dw0 ) || ( w < this->dw0 ) ) {
            return 0.;
        }
    } else {
        if ( ( w_max <= 0. ) || ( w <= 0. ) ) {
            return 0.;
        }
    }


    if ( ( this->fiberType == FT_CAF ) || ( this->fiberType == FT_SAF ) ) {
        theta = fabs( this->computeCrackFibreAngle(gp, i) );
    }

    if ( fiberType == FT_CAF ) { // continuous aligned fibers
        // smooth
        if ( ( this->smoothen ) && ( w_max > this->dw0 ) && ( w_max < this->dw1 ) ) {
            sig_dw1 =  2. *this->Vf *sqrt(this->dw1 *this->Ef * ( 1. + this->eta ) *this->tau_0 / this->Df); // stress in fibers per unit area of concrete
            Dsig_dw1 = this->Vf * sqrt( this->Ef * ( 1. + this->eta ) * this->tau_0 / ( this->Df * this->dw1 ) ); // stress derivative

            x = w_max - this->dw0;
            dw_smooth = this->dw1 - this->dw0;

            C1 = Dsig_dw1 / pow(dw_smooth, 2) - ( 2 * sig_dw1 ) / pow(dw_smooth, 3);
            C2 = ( 3 * sig_dw1 ) / pow(dw_smooth, 2) - Dsig_dw1 / dw_smooth;

            traction_ff = C1 * pow(x, 3) + C2 *pow(x, 2);

            // sharp
        } else {
            // expressed with respect to crack opening
            //
            traction_ff =  2. *this->Vf *sqrt(w_max *this->Ef * ( 1. + this->eta ) *this->tau_0 / this->Df); // stress in fibers per unit area of concrete
        }

        // reflect fibre orientation wrt crack
        traction_ff *= fabs( cos(theta) ) * exp(theta * this->f);

        // reflect damage
        omega = this->computeTempDamage(gp);
        traction_ff *= ( 1. - omega );
    } else if ( fiberType == FT_SAF ) { // short aligned fibers
        omega = this->computeTempDamage(gp);

        // smooth
        if ( ( this->smoothen ) && ( w_max > this->dw0 ) && ( w_max < this->dw1 ) ) {
            sig_dw1 =  this->Vf * ( sqrt(4. * this->Ef * ( 1. + this->eta ) * this->dw1 * this->tau_0 / this->Df) - this->Ef * ( 1. + this->eta ) * this->dw1 / this->Lf );
            Dsig_dw1 = this->Vf * ( sqrt( this->Ef * ( 1. + this->eta ) * this->tau_0 / ( this->Df * this->dw1 ) ) - this->Ef * ( 1. + this->eta ) / this->Lf );

            x = w_max - this->dw0;
            dw_smooth = this->dw1 - this->dw0;

            C1 = Dsig_dw1 / pow(dw_smooth, 2) - ( 2 * sig_dw1 ) / pow(dw_smooth, 3);
            C2 = ( 3 * sig_dw1 ) / pow(dw_smooth, 2) - Dsig_dw1 / dw_smooth;

            traction_ff = C1 * pow(x, 3) + C2 *pow(x, 2);

            // reflect fibre orientation wrt crack
            traction_ff *= fabs( cos(theta) ) * exp(theta * this->f);

            traction_ff *= ( 1. - omega );
        } else if ( w_max < this->w_star ) { // debonding + pullout
            traction_ff = this->Vf * ( sqrt(4. * this->Ef * ( 1. + this->eta ) * w_max * this->tau_0 / this->Df) - this->Ef * ( 1. + this->eta ) * w_max / this->Lf );

            // reflect fibre orientation wrt crack
            traction_ff *= fabs( cos(theta) ) * exp(theta * this->f);

            traction_ff *= ( 1. - omega );
        } else if ( w_max <= this->Lf / 2. ) { // pullout
            factor = this->Vf * this->Lf / this->Df;

            traction_ff = factor * this->computeFiberBond(w_max) * ( 1. - 2. * w_max / this->Lf ) * ( 1. - 2. * w_max / this->Lf );

            // reflect fibre orientation wrt crack
            traction_ff *= fabs( cos(theta) ) * exp(theta * this->f);

            traction_ff *= ( 1. - omega );
        } else { // fully pulled out
            traction_ff = 0.;
        }
    } else if ( fiberType == FT_SRF ) { // short random fibers
        omega = this->computeTempDamage(gp);
        factor = ( 1. - omega ) * this->g * this->Vf * this->Lf / ( 2. * this->Df );

        // smooth
        if ( ( this->smoothen ) && ( w_max > this->dw0 ) && ( w_max < this->dw1 ) ) {
            sig_dw1 = factor * this->tau_0 * ( 2. * sqrt(this->dw1 / this->w_star) - this->dw1 / this->w_star );
            Dsig_dw1 = factor * this->tau_0 * ( 1. / ( this->w_star * sqrt(this->dw1 / this->w_star) ) - 1. / this->w_star );

            x = w_max - this->dw0;
            dw_smooth = this->dw1 - this->dw0;

            C1 = Dsig_dw1 / pow(dw_smooth, 2) - ( 2 * sig_dw1 ) / pow(dw_smooth, 3);
            C2 = ( 3 * sig_dw1 ) / pow(dw_smooth, 2) - Dsig_dw1 / dw_smooth;

            traction_ff = C1 * pow(x, 3) + C2 *pow(x, 2);
        } else if ( w_max < this->w_star ) { // debonding + pullout
            traction_ff = factor * this->tau_0 * ( 2. * sqrt(w_max / this->w_star) - w_max / this->w_star );
        } else if ( w_max <= this->Lf / 2. ) { // pullout
            traction_ff = factor * this->computeFiberBond(w_max) * ( 1. - 2. * w_max / this->Lf ) * ( 1. - 2. * w_max / this->Lf );
        } else { // fully pulled out
            traction_ff = 0.;
        }
    } else {
        OOFEM_ERROR("Unknown fiber type");
    }


    if ( ec < emax ) { // unloading-reloading
        if ( smoothen ) {
            traction_ff *=  pow( ( w - this->dw0 ) / ( w_max - this->dw0 ), this->M );
        } else {
            traction_ff *=  pow(w / w_max, this->M);
        }
    }

    return traction_ff;
}



double
FRCFCM :: computeEffectiveShearModulus(GaussPoint *gp, int shearDirection)
{
    double G, Geff;
    double beta_mf;
    double D2_1, D2_2, D2;
    int crackA, crackB;

    G = this->computeOverallElasticShearModulus();

    if ( this->isIntactForShear(gp, shearDirection) ) {
        Geff = G;
    } else {
        if ( this->shearType == SHR_NONE ) {
            Geff = G;
        } else {
            if ( shearDirection == 4 ) {
                crackA = 2;
                crackB = 3;
            } else if ( shearDirection == 5 ) {
                crackA = 1;
                crackB = 3;
            }  else if ( shearDirection == 6 ) {
                crackA = 1;
                crackB = 2;
            } else {
                OOFEM_ERROR("Unexpected value of index i (4, 5, 6 permitted only)");
                crackA = crackB = 0; // happy compiler
            }

            if ( ( this->isIntact(gp, crackA) ) || ( this->isIntact(gp, crackB) ) ) {
                if ( this->isIntact(gp, crackA) ) {
                    D2 = this->computeD2ModulusForCrack(gp, crackB);
                } else {
                    D2 = this->computeD2ModulusForCrack(gp, crackA);
                }

                beta_mf = D2 / ( D2 + G );
            } else {
                D2_1 = this->computeD2ModulusForCrack(gp, crackA);
                D2_2 = this->computeD2ModulusForCrack(gp, crackB);

                if ( multipleCrackShear ) {
                    beta_mf = 1. / ( 1. + G * ( 1 / D2_1 + 1 / D2_2 ) );
                } else {
                    D2 = min(D2_1, D2_2);
                    beta_mf = D2 / ( D2 + G );
                }
            }

            Geff = G * beta_mf;
        }
    } // not intact for shear

    return Geff;
}


double
FRCFCM :: computeD2ModulusForCrack(GaussPoint *gp, int icrack)
{
    double cos_theta = 0.5; // to reduce Vf for SRF to one half
    double E = this->giveLinearElasticMaterial()->giveYoungsModulus();
    FRCFCMStatus *status = static_cast< FRCFCMStatus * >( this->giveStatus(gp) );
    double crackStrain;
    double D2m, D2f;
    double omega;

    crackStrain = status->giveTempMaxCrackStrain(icrack);

    if ( ( this->isIntact(gp, icrack) ) || ( crackStrain <= 0. ) ) {
        return E * fcm_BIGNUMBER;
    } else {
        if ( ( this->fiberType == FT_CAF ) || ( this->fiberType == FT_SAF ) ) {
            cos_theta = fabs( cos( this->computeCrackFibreAngle(gp, icrack) ) );
        }

        D2m = ConcreteFCM :: computeD2ModulusForCrack(gp, icrack);
        D2m *= ( 1. - this->Vf );

        omega = this->computeTempDamage(gp);

        // fiber shear stiffness is not influenced by the number of parallel cracks
        D2f = ( 1. - omega ) * this->Vf * cos_theta * this->kfib * this->Gfib / crackStrain;

        D2f = min(E * fcm_BIGNUMBER, D2f);

        return D2m + D2f;
    }
}

// the same function as "computeD2ModulusForCrack", only without current value of fiber damage.
double
FRCFCM :: estimateD2ModulusForCrack(GaussPoint *gp, int icrack)
{
    double cos_theta = 0.5; // to reduce Vf for SRF to one half
    double E = this->giveLinearElasticMaterial()->giveYoungsModulus();
    FRCFCMStatus *status = static_cast< FRCFCMStatus * >( this->giveStatus(gp) );
    double crackStrain;
    double D2m, D2f;
    double omega;

    crackStrain = status->giveTempMaxCrackStrain(icrack);

    if ( ( this->isIntact(gp, icrack) ) || ( crackStrain <= 0. ) ) {
        return E * fcm_BIGNUMBER;
    } else {
        if ( ( this->fiberType == FT_CAF ) || ( this->fiberType == FT_SAF ) ) {
            cos_theta = fabs( cos( this->computeCrackFibreAngle(gp, icrack) ) );
        }

        D2m = ConcreteFCM :: computeD2ModulusForCrack(gp, icrack);
        D2m *= ( 1. - this->Vf );

        omega = status->giveDamage();

        // fiber shear stiffness is not influenced by the number of parallel cracks
        D2f = ( 1. - omega ) * this->Vf * cos_theta * this->kfib * this->Gfib / crackStrain;

        D2f = min(E * fcm_BIGNUMBER, D2f);

        return D2m + D2f;
    }
}


double
FRCFCM :: computeTempDamage(GaussPoint *gp) {
    // we assume that fibre damage is the same for all crack planes
    FRCFCMStatus *status = static_cast< FRCFCMStatus * >( this->giveStatus(gp) );

    double omega = 0.;
    double gammaCrack = 0.;

    double slip, opening;

    int numberOfActiveCracks = status->giveNumberOfTempCracks();


    if ( fiberDamageType != FDAM_NONE ) { // fiber damage is allowed
        for ( int i = 1; i <= numberOfActiveCracks; i++ ) {
            if ( !this->isIntact(gp, i) ) {
                opening =  this->computeMaxNormalCrackOpening(gp, i);

                //	if (opening > 0.) {
                if ( opening > this->fibreActivationOpening ) {
                    slip =  this->computeShearSlipOnCrack(gp, i);
                    gammaCrack = max(gammaCrack, slip / opening);
                }
            } // initiation condition
        } // active crack loop

        if ( fiberDamageType == FDAM_GammaCrackLin ) {
            omega = min(gammaCrack / this->gammaCrackFail, 1.);
        } else if ( fiberDamageType == FDAM_GammaCrackExp ) {
            omega = 1. - exp(-gammaCrack / this->gammaCrackFail);
        } else {
            OOFEM_ERROR("Unknown FiberDamageType");
        }

        omega = max( omega, status->giveDamage() );
        status->setTempDamage(omega);
    }

#if DEBUG
    if ( omega > 0.9 ) {
        OOFEM_WARNING( "High value of damage in Element %d", gp->giveElement()->giveNumber() );
    }
#endif

    return omega;
}


double
FRCFCM :: maxShearStress(GaussPoint *gp, int shearDirection)
{
    double maxTau_m;
    double minTau_f;
    double E = this->giveLinearElasticMaterial()->giveYoungsModulus();
    double crackStrain;
    double gamma_cr;
    double omega;
    double cos_theta = 0.5; // to reduce Vf for SRF to one half
    MaterialMode mMode = gp->giveMaterialMode();

    FRCFCMStatus *status = static_cast< FRCFCMStatus * >( this->giveStatus(gp) );

    int crackA, crackB;
    int icrack;

    // max shear in matrix
    maxTau_m =  ConcreteFCM :: maxShearStress(gp, shearDirection);
    maxTau_m *= ( 1. - this->Vf );

    // for now we simply compute the least allowable stress as a product of crack shear stiffness (fiber contribution) * max shear strain

    if ( shearDirection == 4 ) {
        crackA = 2;
        crackB = 3;
    } else if ( shearDirection == 5 ) {
        crackA = 1;
        crackB = 3;
    }  else if ( shearDirection == 6 ) {
        crackA = 1;
        crackB = 2;
    } else {
        OOFEM_ERROR("Unexpected value of index i (4, 5, 6 permitted only)");
    }

    minTau_f = E * fcm_BIGNUMBER;

    if ( mMode == _PlaneStress ) {
        gamma_cr = fabs( status->giveTempMaxCrackStrain(3) );
    } else {
        gamma_cr = fabs( status->giveTempMaxCrackStrain(shearDirection) );
    }

    for ( int i = 1; i <= 2; i++ ) {
        if ( i == 1 ) {
            icrack = crackA;
        } else { // i == 2
            icrack = crackB;
        }

        crackStrain = status->giveTempMaxCrackStrain(icrack);

        if ( ( this->isIntact(gp, icrack) ) || ( crackStrain <= 0. ) ) {
            minTau_f = min(minTau_f, E * fcm_BIGNUMBER);
        } else {
            if ( ( this->fiberType == FT_CAF ) || ( this->fiberType == FT_SAF ) ) {
                cos_theta = fabs( cos( this->computeCrackFibreAngle(gp, icrack) ) );
            }

            omega = this->computeTempDamage(gp);

            minTau_f = min(minTau_f, gamma_cr * ( 1. - omega ) * this->Vf * cos_theta * this->kfib * this->Gfib / crackStrain);
        }
    }

    return maxTau_m + minTau_f;
}



// computes crack spacing at saturated state
// random tensile strength is not supported here
// effect of fiber inclination is not considered
double
FRCFCM :: computeCrackSpacing() {
    double x_CA, lambda;
    double x = 0.;

    x_CA = ( 1. - this->Vf ) * this->Ft * this->Df / ( 4. * this->Vf * this->tau_0 );

    if ( fiberType == FT_CAF ) { // continuous aligned fibers
        x = x_CA;
    } else if ( fiberType == FT_SAF ) { // short aligned fibers
        x = 0.5 * sqrt(this->Lf * this->Lf - 4. * this->Lf * x_CA);
    } else if ( fiberType == FT_SRF ) { // short random fibers
        lambda = ( 2. / M_PI ) * ( 4. + this->f * this->f ) / ( 1. + exp(M_PI * f / 2.) );
        x = 0.5 * ( this->Lf - sqrt(this->Lf * this->Lf - 2. * M_PI * this->Lf * lambda * x_CA) );
    } else {
        OOFEM_ERROR("Unknown fiber type");
    }

    return x;
}


bool
FRCFCM :: isStrengthExceeded(const FloatMatrix &base, GaussPoint *gp, TimeStep *tStep, int iCrack, double trialStress)
{
    double Em, sigma_m;

    if ( trialStress <= 0. ) {
        return false;
    }

    Em = this->giveLinearElasticMaterial()->giveYoungsModulus();

    // matrix is stiffer -> carries higher stress
    if ( ( this->Ef <= Em ) && ( trialStress > this->giveTensileStrength(gp) ) ) {
        return true;
    } else {
        sigma_m = trialStress / ( 1. + this->Vf * ( this->Ef / Em - 1. ) );

        if ( sigma_m > this->giveTensileStrength(gp) ) {
            return true;
        } else {
            return false;
        }
    }
}


double
FRCFCM :: computeShearStiffnessRedistributionFactor(GaussPoint *gp, int ithCrackPlane, int jthCrackDirection)
{
    double factor_ij, D2_i, D2_j;

    // slightly modified version here. The problem is recursive calling of damage & slip evaluation for the FRCFCM material
    D2_i = this->estimateD2ModulusForCrack(gp, ithCrackPlane);
    D2_j = this->estimateD2ModulusForCrack(gp, jthCrackDirection);

    factor_ij = D2_j / ( D2_i + D2_j );

    return factor_ij;
}


int
FRCFCM :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FRCFCMStatus *status = static_cast< FRCFCMStatus * >( this->giveStatus(gp) );

    // nominal stress in fibers
    if  ( type == IST_FiberStressLocal ) {
        answer.resize(1);
        answer.at(1) = this->computeStressInFibersInCracked(gp, status->giveCrackStrain(1), 1);
        return 1;
    } else if  ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.at(1) = status->giveDamage();
        return 1;
    } else {
        return ConcreteFCM :: giveIPValue(answer, gp, type, tStep);
    }
}


// computes overall stiffness of the composite material: the main purpose of this method is to adjust the stiffness given by the linear elastic material which corresponds to the matrix. The same method is used by all fiber types.
double
FRCFCM :: computeOverallElasticStiffness(void) {
    double stiffness = 0.;

    double Em = this->giveLinearElasticMaterial()->giveYoungsModulus();

    if ( this->fiberType == FT_CAF ) { // continuous aligned fibers
        stiffness = this->Vf * this->Ef + ( 1. - this->Vf ) * Em;
    } else if ( this->fiberType == FT_SAF ) { // short aligned fibers
        stiffness = this->Vf * this->Ef + ( 1. - this->Vf ) * Em;
    } else if ( this->fiberType == FT_SRF ) { // short random fibers
        stiffness = this->Vf * this->Ef + ( 1. - this->Vf ) * Em;
    } else {
        OOFEM_ERROR("Unknown fiber type");
    }

    return stiffness;
}


///////////////////////////////////////////////////////////////////
//                      CONCRETE FCM STATUS                     ///
///////////////////////////////////////////////////////////////////


FRCFCMStatus :: FRCFCMStatus(int n, Domain *d, GaussPoint *gp) :
    ConcreteFCMStatus(n, d, gp)
{
    damage = tempDamage = 0.0;
}


FRCFCMStatus :: ~FRCFCMStatus()
{ }



void
FRCFCMStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    ConcreteFCMStatus :: printOutputAt(file, tStep);


    fprintf(file, "damage status { ");
    if ( this->damage > 0.0 ) {
        fprintf(file, "damage %f ", this->damage);
    }
    fprintf(file, "}\n");
}





void
FRCFCMStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
//
{
    ConcreteFCMStatus :: initTempStatus();

    this->tempDamage = this->damage;
}



void
FRCFCMStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables correspond to newly reched equilibrium.
//
{
    ConcreteFCMStatus :: updateYourself(tStep);

    this->damage = this->tempDamage;
}



contextIOResultType
FRCFCMStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = ConcreteFCMStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.write(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
FRCFCMStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = ConcreteFCMStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.read(damage) ) {
        return CIO_IOERR;
    }

    return CIO_OK; // return succes
}
} // end namespace oofem
