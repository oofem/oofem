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

#include "concretefcm.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"
#include "contextioerr.h"

namespace oofem {
REGISTER_Material(ConcreteFCM);

ConcreteFCM :: ConcreteFCM(int n, Domain *d) : FCMMaterial(n, d), RandomMaterialExtensionInterface()
{
    Gf = Ft = 0.;
    softType = ST_Unknown;
    shearType = SHR_Unknown;
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
}

IRResultType
ConcreteFCM :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = FCMMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    result = linearElasticMaterial->initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    result = RandomMaterialExtensionInterface :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    // type of TENSION SOFTENING
    int softening = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, softening, _IFT_ConcreteFCM_softType);
    switch ( softening ) {
    case 0:
        softType = ST_NONE;
        Ft = 1.e50;
        Gf = 1.e50;
        break;

    case 1:
        softType = ST_Exponential;
        IR_GIVE_FIELD(ir, Gf, _IFT_ConcreteFCM_gf);
        IR_GIVE_FIELD(ir, Ft, _IFT_ConcreteFCM_ft);
        break;

    case 2:
        softType = ST_Linear;
        IR_GIVE_FIELD(ir, Gf, _IFT_ConcreteFCM_gf);
        IR_GIVE_FIELD(ir, Ft, _IFT_ConcreteFCM_ft);
        break;

    case 3:
        softType = ST_Hordijk;
        IR_GIVE_FIELD(ir, Gf, _IFT_ConcreteFCM_gf);
        IR_GIVE_FIELD(ir, Ft, _IFT_ConcreteFCM_ft);
        break;

    case 4:
        softType = ST_UserDefinedCrack;
        IR_GIVE_FIELD(ir, Ft, _IFT_ConcreteFCM_ft);
        IR_GIVE_FIELD(ir, soft_w, _IFT_ConcreteFCM_soft_w);
        IR_GIVE_FIELD(ir, soft_function_w, _IFT_ConcreteFCM_soft_function_w);

        if ( soft_w.at(1) > 0. ) {
            OOFEM_WARNING("The first entry in soft_w must be 0.");
            return IRRT_BAD_FORMAT;
        }

        if ( soft_w.giveSize() != soft_function_w.giveSize() ) {
            OOFEM_WARNING("the size of 'soft_w' and 'soft(w)' must be the same");
            return IRRT_BAD_FORMAT;
        }
        break;

    case 5:
        softType = ST_LinearHardeningStrain;
        IR_GIVE_FIELD(ir, H, _IFT_ConcreteFCM_H);
        this->eps_f = 1.;
        IR_GIVE_OPTIONAL_FIELD(ir, eps_f, _IFT_ConcreteFCM_eps_f);
        IR_GIVE_FIELD(ir, Ft, _IFT_ConcreteFCM_ft);

        if ( this->giveCrackSpacing() > 0. ) {
            OOFEM_ERROR("crack spacing must not be defined in strain-defined materials (does not make sense)");
        }

        break;

    case 6:
        softType = ST_UserDefinedStrain;
        IR_GIVE_FIELD(ir, Ft, _IFT_ConcreteFCM_ft);
        IR_GIVE_FIELD(ir, soft_eps, _IFT_ConcreteFCM_soft_eps);
        IR_GIVE_FIELD(ir, soft_function_eps, _IFT_ConcreteFCM_soft_function_eps);

        if ( soft_eps.at(1) > 0. ) {
            OOFEM_WARNING("The first entry in soft_eps must be 0.");
            return IRRT_BAD_FORMAT;
        }

        if ( soft_eps.giveSize() != soft_function_eps.giveSize() ) {
            OOFEM_WARNING("the size of 'soft_eps' and 'soft(eps)' must be the same");
            return IRRT_BAD_FORMAT;
        }

        if ( this->giveCrackSpacing() > 0. ) {
            OOFEM_ERROR("crack spacing must not be defined in strain-defined materials (does not make sense)");
        }

        break;

    default:
        OOFEM_WARNING("Softening type number %d is unknown", softType);
        return IRRT_BAD_FORMAT;
    }


    // type of SHEAR STIFFNESS REDUCTION
    int shear = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, shear, _IFT_ConcreteFCM_shearType);
    switch ( shear ) {
    case 0:
        shearType = SHR_NONE;

        break;

    case 1:
        shearType = SHR_Const_ShearRetFactor;
        beta = 0.01;
        IR_GIVE_OPTIONAL_FIELD(ir, beta, _IFT_ConcreteFCM_beta);

        if ( beta >= 1. ) {
            OOFEM_ERROR("Parameter beta %f must be smaller than one", beta);
        }

        break;

    case 2:
        // todo: resolve problems with convergence when crack is wide open
        // this does not happen when beta is used instead sf
        shearType = SHR_Const_ShearFactorCoeff;
        sf = 20.;
        IR_GIVE_OPTIONAL_FIELD(ir, sf, _IFT_ConcreteFCM_sf);
        break;

    case 3: // w is the maximum reached crack opening, not the current value!
        shearType = SHR_UserDefined_ShearRetFactor;
        IR_GIVE_FIELD(ir, beta_w, _IFT_ConcreteFCM_beta_w);
        IR_GIVE_FIELD(ir, beta_function, _IFT_ConcreteFCM_beta_function);

        if ( beta_w.giveSize() != beta_function.giveSize() ) {
            OOFEM_WARNING("the size of 'beta_w' and 'beta(w)' must be the same");
            return IRRT_BAD_FORMAT;
        }
        break;


    default:
        OOFEM_WARNING("Shear stiffness reduction type number %d is unknown", softType);
        return IRRT_BAD_FORMAT;
    }


    // type of SHEAR STRENGTH CUTOFF
    shear = 0;

    IR_GIVE_OPTIONAL_FIELD(ir, shear, _IFT_ConcreteFCM_shearStrengthType);
    switch ( shear ) {
    case 0:
        shearStrengthType = SHS_NONE;
        break;

    case 1:
        shearStrengthType = SHS_Const_Ft;
        break;

    case 2:
        shearStrengthType = SHS_Collins_Interlock;

        IR_GIVE_FIELD(ir, fc, _IFT_ConcreteFCM_fc);
        IR_GIVE_FIELD(ir, ag, _IFT_ConcreteFCM_ag);
        IR_GIVE_FIELD(ir, lengthScale, _IFT_ConcreteFCM_lengthScale);

        break;

    default:
        OOFEM_WARNING("Shear stiffness reduction type number %d is unknown", softType);
        return IRRT_BAD_FORMAT;
    }

    return IRRT_OK;
}


double
ConcreteFCM :: giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp, int i)
//
// returns current cracking modulus according to crackStrain for i-th
// crackplane
//
{
    ConcreteFCMStatus *status = static_cast< ConcreteFCMStatus * >( this->giveStatus(gp) );

    double Cf = 0.;
    double ef, emax, c1, c2, Le, Ncr, w;

    double Ft = this->giveTensileStrength(gp);
    double Gf = this->giveFractureEnergy(gp);
    double ec =  status->giveTempCrackStrain(i);
    double E = this->computeOverallElasticStiffness();

    if ( status->giveTempCrackStatus(i) == pscm_NONE ) {
        return E * fcm_BIGNUMBER;
    }

    // very high stiffness also for compression
    if ( status->giveTempCrackStatus(i) == pscm_CLOSED ) {
        return E * fcm_BIGNUMBER;
    }

    emax = max(status->giveMaxCrackStrain(i), ec);

    Le = status->giveCharLength(i);

    // modification to take into account crack spacing - makes sense only in crack-opening-based materials
    Ncr = this->giveNumberOfCracksInDirection(gp, i);

    ec /= Ncr;
    emax /= Ncr;


    if ( ( rMode == TangentStiffness ) && ( ec >= emax ) ) { // tangent and NOT unloading
        if ( softType == ST_NONE ) {
            OOFEM_ERROR("For unknown reason the tensile strength has been exceeded and cracking has been initiated!");
        } else if ( softType == ST_Exponential ) {
            ef = Gf / ( Ft * Le );

            if ( emax == 0. ) { // just initiated and closing crack
                Cf = -Ft / ef;
            } else { // further softening - negative stiffness
                Cf = -Ft / ef *exp(-ec / ef);
            }
        } else if ( softType == ST_Linear ) {
            // fracturing strain at zero traction
            ef = 2. * Gf / ( Ft * Le );

            if ( ( ec >= ef )  || ( emax >= ef ) ) {
                // fully open crack - no stiffness
                Cf = 0.;
            } else { // further softening - negative stiffness
                Cf = -Ft / ef;
            }
        } else if ( softType == ST_Hordijk ) {
            c1 = 3.;
            c2 = 6.93;

            ef = 5.14 * Gf / Ft / Le;

            if ( ( ec >= ef )  || ( emax >= ef ) ) {
                Cf = 0.;
            }  else if ( emax == 0. ) { // just initiated and closing crack
                Cf = Ft * ( -c2 / ef - ( exp(-c2) * ( c1 * c1 * c1 + 1. ) ) / ef );
            } else { // further softening - negative stiffness
                Cf = Ft * ( ( 3. * c1 * c1 * c1 * ec * ec * exp(-( c2 * ec ) / ef) ) /  ( ef * ef * ef ) - ( c2 * exp(-( c2 * ec ) / ef) * ( c1 * c1 * c1 * ec * ec * ec + ef * ef * ef ) ) / ( ef * ef * ef * ef ) - ( exp(-c2) * ( c1 * c1 * c1 + 1. ) ) / ef );
            }
        } else if ( softType == ST_UserDefinedCrack ) {
            w = ec * Le;

            if ( w > soft_w.giveSize() ) {
                Cf = 0.;
                OOFEM_WARNING("Crack width is larger than the last user-defined value in the traction-opening law.");
            } else if ( emax == 0. ) {
                Cf = this->Ft * ( soft_function_w.at(2) - soft_function_w.at(1) ) / ( ( soft_w.at(2) - soft_w.at(1) ) / Le );
            } else { // softening
                for ( int i = 1; i <= soft_w.giveSize(); i++ ) {
                    if ( ( w - soft_w.at(i) ) < fcm_SMALL_STRAIN ) {
                        if ( i == 1 ) {
                            Cf =  this->Ft * ( soft_function_w.at(2) - soft_function_w.at(1) ) /  ( ( soft_w.at(2) - soft_w.at(1) ) / Le );
                            break;
                        } else {
                            Cf =  this->Ft * ( soft_function_w.at(i) - soft_function_w.at(i - 1) ) /  ( ( soft_w.at(i) - soft_w.at(i - 1) ) / Le );
                            break;
                        }
                    }
                }
            }
        } else if ( softType == ST_LinearHardeningStrain ) {
            if ( ( ec >= this->eps_f )  || ( emax >= this->eps_f ) ) {
                Cf = 0.;
            } else {
                // hardening
                Cf = this->H;
            }
        } else if ( softType == ST_UserDefinedStrain ) {
            if ( ec > soft_eps.giveSize() ) {
                Cf = 0.;
                OOFEM_WARNING("Crack strain is larger than the last user-defined value in the stress-crack-strain law.");
            } else if ( emax == 0. ) {
                Cf = this->Ft * ( soft_function_eps.at(2) - soft_function_eps.at(1) ) / ( soft_eps.at(2) - soft_eps.at(1) );
            } else { // softening
                for ( int i = 1; i <= soft_eps.giveSize(); i++ ) {
                    if ( ( ec - soft_eps.at(i) ) < fcm_SMALL_STRAIN ) {
                        if ( i == 1 ) {
                            Cf =  this->Ft * ( soft_function_eps.at(2) - soft_function_eps.at(1) ) /  ( soft_eps.at(2) - soft_eps.at(1) );
                            break;
                        } else {
                            Cf =  this->Ft * ( soft_function_eps.at(i) - soft_function_eps.at(i - 1) ) /  ( soft_eps.at(i) - soft_eps.at(i - 1) );
                            break;
                        }
                    }
                }
            }
        } else {
            OOFEM_ERROR("Unknown Softening Mode");
        }
    } else { // secant stiffness OR tangent stiffness and unloading
        if ( emax <= 0. ) { // just initiated crack that wants to close
            return E * fcm_BIGNUMBER;
        }

        Cf = ConcreteFCM :: giveNormalCrackingStress(gp, emax * Ncr, i); // gets stress
        Cf /= ( emax ); // converted into stiffness
    }

    // to take into account crack-spacing
    return Cf / Ncr;
}

double
ConcreteFCM :: giveNormalCrackingStress(GaussPoint *gp, double ec, int i)
//
// returns receivers Normal Stress in crack i  for given cracking strain
//
{
    ConcreteFCMStatus *status = static_cast< ConcreteFCMStatus * >( this->giveStatus(gp) );

    double traction;
    double Le, ef, emax, w, E, c1, c2;
    double Gf = this->giveFractureEnergy(gp);
    double Ft = this->giveTensileStrength(gp);

    if ( status->giveTempCrackStatus(i) == pscm_CLOSED ) {
        E = this->computeOverallElasticStiffness();
        return E * ec;
    }

    Le = status->giveCharLength(i);
    emax = max(status->giveMaxCrackStrain(i), ec);

    // to take into account crack spacing
    ec /= this->giveNumberOfCracksInDirection(gp, i);
    emax /= this->giveNumberOfCracksInDirection(gp, i);


    if ( softType == ST_NONE ) {
        OOFEM_ERROR("For unknown reason the tensile strength has been exceeded and cracking has been initiated!");
    } else if ( softType == ST_Exponential ) {
        ef = Gf / ( Ft * Le );
        if ( ec >= emax ) {
            // further softening
            traction = Ft * exp(-ec / ef);
        } else {
            // crack closing
            // or unloading or reloading regime
            traction = Ft * ec / emax *exp(-emax / ef);
        }
    } else if ( softType  == ST_Linear ) {
        // fracturing strain at zero traction
        ef = 2. * Gf / ( Ft * Le );

        if ( ( ec >= ef )  || ( emax >= ef ) ) {
            // fully open crack - no stiffness
            traction = 0.;
        } else if ( ec >= emax ) {
            // further softening
            traction = Ft - Ft * ec / ef;
        } else if ( ec <= 0. ) {
            traction = 0.;
        } else {
            traction = Ft * ec * ( ef - emax ) / ( emax * ef );
        }
    } else if ( softType == ST_Hordijk ) {
        c1 = 3.;
        c2 = 6.93;

        ef = 5.14 * Gf / Ft / Le;

        if ( ( ec >= ef )  || ( emax >= ef ) ) {
            // fully open crack - no stiffness
            traction = 0.;
        } else if ( ec >= emax ) {
            // further softening
            traction = Ft * ( ( 1. + pow( ( c1 * ec / ef ), 3. ) ) * exp(-c2 * ec / ef) - ec / ef * ( 1. + c1 * c1 * c1 ) * exp(-c2) );
        } else if ( ec <= 0. ) {
            traction = 0.;
        } else {
            traction = Ft * ec / emax * ( ( 1. + pow( ( c1 * emax / ef ), 3. ) ) * exp(-c2 * emax / ef) - emax / ef * ( 1. + c1 * c1 * c1 ) * exp(-c2) );
        }
    } else if ( softType == ST_UserDefinedCrack ) {
        w = ec * Le;


        if ( w > soft_w.giveSize() ) {
            traction = 0.;
            OOFEM_WARNING("Crack width is larger than the last user-defined value in the traction-opening law.");
        } else if ( ec >= emax ) { // softening
            for ( int i = 1; i <= soft_w.giveSize(); i++ ) {
                if ( ( w - soft_w.at(i) ) < fcm_SMALL_STRAIN ) {
                    if ( i == 1 ) { // eps \approx 0 bound
                        traction = soft_function_w.at(i) * this->Ft;
                        break;
                    } else {
                        traction = soft_function_w.at(i - 1) +  ( soft_function_w.at(i) - soft_function_w.at(i - 1) ) / ( soft_w.at(i) - soft_w.at(i - 1) ) * ( w - soft_w.at(i - 1) );
                        traction *= this->Ft;
                        break;
                    }
                }
            }
        } else if ( ec <= 0. ) { // closing
            traction = 0.;
        } else { //unlo-relo
            w = emax * Le; // lower stiffness

            for ( int i = 1; i <= soft_w.giveSize(); i++ ) {
                if ( ( w - soft_w.at(i) ) < fcm_SMALL_STRAIN ) {
                    if ( i == 1 ) { // eps \approx 0
                        traction = soft_function_w.at(i) * this->Ft * ec / emax;
                        break;
                    } else {
                        traction = soft_function_w.at(i - 1) +  ( soft_function_w.at(i) - soft_function_w.at(i - 1) ) / ( soft_w.at(i) - soft_w.at(i - 1) ) * ( w - soft_w.at(i - 1) );
                        traction *= this->Ft * ec / emax;
                        break;
                    }
                }
            }
        }
    } else if ( softType == ST_LinearHardeningStrain ) {
        if ( ( ec >= this->eps_f )  || ( emax >= this->eps_f ) ) {
            traction = 0.;
        } else if ( emax == 0. ) {
            traction = this->Ft;
        } else {
            // unloading or reloading regime
            // emax * Le = crack width
            //      traction = ( this->Ft + emax * Le * this->H ) * ec / emax;
            // written in terms of strain and not crack width
            traction = ( this->Ft + emax * this->H ) * ec / emax;
        }
    } else if ( softType == ST_UserDefinedStrain ) {
        if ( emax > soft_eps.giveSize() ) {
            traction = 0.;
            OOFEM_WARNING("Cracking strain is larger than the last user-defined value in the traction-cracking-strain law.");
        } else if ( ec >= emax ) { // softening
            for ( int i = 1; i <= soft_eps.giveSize(); i++ ) {
                if ( ( ec - soft_eps.at(i) ) < fcm_SMALL_STRAIN ) {
                    if ( i == 1 ) { // eps \approx 0 bound
                        traction = soft_function_eps.at(i) * this->Ft;
                        break;
                    } else {
                        traction = soft_function_eps.at(i - 1) +  ( soft_function_eps.at(i) - soft_function_eps.at(i - 1) ) / ( soft_eps.at(i) - soft_eps.at(i - 1) ) * ( ec - soft_eps.at(i - 1) );
                        traction *= this->Ft;
                        break;
                    }
                }
            }
        } else if ( ec <= 0. ) { // closing
            traction = 0.;
        } else { //unlo-relo
            for ( int i = 1; i <= soft_eps.giveSize(); i++ ) {
                if ( ( emax - soft_eps.at(i) ) < fcm_SMALL_STRAIN ) {
                    if ( i == 1 ) { // eps \approx 0
                        traction = soft_function_eps.at(i) * this->Ft * ec / emax;
                        break;
                    } else {
                        traction = soft_function_eps.at(i - 1) +  ( soft_function_eps.at(i) - soft_function_eps.at(i - 1) ) / ( soft_eps.at(i) - soft_eps.at(i - 1) ) * ( emax - soft_eps.at(i - 1) );
                        traction *= this->Ft * ec / emax;
                        break;
                    }
                }
            }
        }
    } else {
        OOFEM_ERROR("Unknown Softening Mode");
    }

    return traction;
}

double
ConcreteFCM :: computeEffectiveShearModulus(GaussPoint *gp, int shearDirection)

{
    double G, Geff, D2tot, N;
    int crackA, crackB;

    G = this->computeOverallElasticShearModulus();

    if ( this->isIntactForShear(gp, shearDirection) ) {
        Geff = G;
    } else {
        if ( this->shearType == SHR_NONE ) {
            Geff = G;
        } else if ( this->shearType == SHR_Const_ShearRetFactor ) {
            if ( multipleCrackShear ) {
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

                // total number of parallel cracks in both directions
                N = this->giveNumberOfCracksInDirection(gp, crackA) + this->giveNumberOfCracksInDirection(gp, crackB);
            } else {
                N = this->giveNumberOfCracksForShearDirection(gp, shearDirection);
            }

            Geff = G * this->beta / ( N - this->beta * ( N - 1 ) ); // for N = 1... Geff = G * beta
        } else if ( this->shearType == SHR_Const_ShearFactorCoeff ) {
            // contributions from shear terms in: De - De (De + Dcr)^-1 De
            D2tot = this->computeTotalD2Modulus(gp, shearDirection);
            // number of parallel cracks already taken into account
            Geff = G * D2tot / ( G + D2tot );
        } else if ( this->shearType == SHR_UserDefined_ShearRetFactor ) {
            // contributions from shear terms in: De - De (De + Dcr)^-1 De
            D2tot = this->computeTotalD2Modulus(gp, shearDirection);
            // number of parallel cracks already taken into account
            Geff = G * D2tot / ( G + D2tot );
        } else {
            OOFEM_ERROR("Unknown Shear Mode");
            Geff = 0.; // happy compiler
        }
    }

    return Geff;
}


double
ConcreteFCM :: computeD2ModulusForCrack(GaussPoint *gp, int icrack)
{
    double D2, N, G, E, w, beta_m;

    if ( icrack >= 4 ) {
        OOFEM_ERROR("Unexpected crack number");
    }

    N = this->giveNumberOfCracksInDirection(gp, icrack);

    if ( this->isIntact(gp, icrack) || ( this->shearType == SHR_NONE ) ) {
        E = this->computeOverallElasticStiffness();
        D2 = E * fcm_BIGNUMBER;
    } else {
        if ( shearType == SHR_Const_ShearRetFactor ) {
            G = this->computeOverallElasticShearModulus();
            D2 = G * this->beta / ( 1. - this->beta );
            D2 /= N;
        } else if ( shearType == SHR_Const_ShearFactorCoeff ) {
            // for the number of parallel cracks is taken care of in "giveCrackingModulus"
            // however, it is essential to reduce the final modulus due to the presence of multiple cracks
            D2 = this->sf * this->giveCrackingModulus(SecantStiffness, gp, icrack);
            D2 /= N;
        } else if ( this->shearType == SHR_UserDefined_ShearRetFactor ) {
            // for the number of parallel cracks is taken care of in the evaluation of normalcrackopening
            // however, it is essential to reduce the final modulus due to the presence of multiple cracks

            G = this->computeOverallElasticShearModulus();
            w = max( this->computeNormalCrackOpening(gp, icrack), this->computeMaxNormalCrackOpening(gp, icrack) );

            ///  shear retention factor = beta(w)

            beta_m = 0.;
            for ( int i = 1; i <= this->beta_w.giveSize(); i++ ) {
                if ( ( w - this->beta_w.at(i) ) < fcm_SMALL_STRAIN ) {
                    if ( i == 1 ) {
                        beta_m = this->beta_function.at(i);
                    } else {
                        beta_m = this->beta_function.at(i - 1) + ( this->beta_function.at(i) - this->beta_function.at(i - 1) ) / ( this->beta_w.at(i) - this->beta_w.at(i - 1) ) * ( w - this->beta_w.at(i - 1) );
                    }

                    beta_m = min(1., beta_m);
                    break;
                }
            }

            if ( beta_m >= 1. ) {
                D2 = G * fcm_BIGNUMBER;
            } else {
                D2 = G * beta_m / ( 1. - beta_m );
                D2 /= N;
            }
        } else {
            D2 = 0;
            OOFEM_ERROR("Unknown Softening Mode");
        }
    }

    return D2;
}


void
ConcreteFCM :: checkSnapBack(GaussPoint *gp, int i)
{
    ConcreteFCMStatus *status = static_cast< ConcreteFCMStatus * >( this->giveStatus(gp) );

    double E, Gf, Ft, Le;
    double ef = 0.;
    double slope;
    double efail, c1, c2;

    Le = status->giveCharLength(i);

    Gf = this->giveFractureEnergy(gp);
    Ft = this->giveTensileStrength(gp);

    E = this->computeOverallElasticStiffness();

    if ( softType == ST_NONE ) {
        OOFEM_ERROR("For unknown reason the tensile strength has been exceeded and cracking has been initiated!");
    } else if ( softType == ST_Exponential ) {
        ef = Gf / ( Ft * Le );
    } else if ( softType == ST_Linear ) {
        ef = 2. * Gf / ( Ft * Le );
    } else if ( softType == ST_Hordijk ) {
        efail =  5.14 * Gf / Ft / Le;
        c1 = 3.;
        c2 = 6.93;
        // ef = ft / slope of the softening branch at the peak stress
        ef = 1. / ( c2 / efail + ( exp(-c2) * ( c1 * c1 * c1 + 1. ) ) / efail );
    } else if ( softType == ST_UserDefinedCrack ) {
        ef = 1.e20;

        for ( int i = 1; i < soft_w.giveSize(); i++ ) {
            slope = ( soft_function_w.at(i + 1) - soft_function_w.at(i) ) / ( ( soft_w.at(i + 1) - soft_w.at(i) ) / Le );
            if ( slope < 0. ) {
                ef = min(fabs(1. / slope), ef);
            }
        }
    } else if ( softType == ST_LinearHardeningStrain ) {
        return;
    } else if ( softType == ST_UserDefinedStrain ) {
        ef = 1.e20;

        for ( int i = 1; i < soft_eps.giveSize(); i++ ) {
            slope = ( soft_function_eps.at(i + 1) - soft_function_eps.at(i) ) / ( soft_eps.at(i + 1) - soft_eps.at(i) );
            if ( slope < 0. ) {
                ef = min(fabs(1. / slope), ef);
            }
        }
    } else {
        OOFEM_ERROR("Unknown Softening Mode");
    }

    if ( ef <= Ft / E ) {
        OOFEM_ERROR("ef %e < e0 %e, this leads to material snapback in element %d, characteristic length %f", ef, Ft / E, gp->giveElement()->giveNumber(), Le);
    }
}


double
ConcreteFCM :: maxShearStress(GaussPoint *gp, int i) {
    int dir_1, dir_2;
    double maxTau, crackOpening, wmax, scale;

    if ( this->isIntactForShear(gp, i) ) {
        return Ft * fcm_BIGNUMBER;
    }

    if ( this->shearStrengthType == SHS_NONE ) {
        maxTau = Ft * fcm_BIGNUMBER;
    } else if ( this->shearStrengthType == SHS_Const_Ft ) {
        maxTau = Ft;
    } else if ( this->shearStrengthType == SHS_Collins_Interlock ) {
        if ( i == 4 ) { // y-z
            dir_1 = 2;
            dir_2 = 3;
        } else if ( i == 5 ) { // x-z
            dir_1 = 1;
            dir_2 = 3;
        } else if ( i == 6 ) { // x-y
            dir_1 = 1;
            dir_2 = 2;
        } else {
            OOFEM_ERROR("Unexpected number for shear stress (must be either 4, 5 or 6).");
            dir_1 = dir_2 = 0; // happy compiler
        }

        // from temporary strain
        crackOpening = max( this->computeNormalCrackOpening(gp, dir_1), this->computeNormalCrackOpening(gp, dir_2) );

        // from the entire history
        wmax = max( this->computeMaxNormalCrackOpening(gp, dir_1), this->computeMaxNormalCrackOpening(gp, dir_2) );

        crackOpening = max(wmax, crackOpening);

        scale = 1000. / lengthScale;
        maxTau = 0.18 * sqrt(this->fc) / ( 0.31 + 24. * crackOpening * scale / ( this->ag * scale + 16. ) );
    } else {
        OOFEM_ERROR("Unexpected shearStrengthType");
        maxTau = 0.; // happy compiler
    }

    return maxTau;
}


double
ConcreteFCM :: give(int aProperty, GaussPoint *gp)
{
    double answer;
    if ( RandomMaterialExtensionInterface :: give(aProperty, gp, answer) ) {
        return answer;
    } else if ( aProperty == gf_ID ) {
        return this->Gf;
    } else if ( aProperty == ft_strength ) {
        return this->Ft;
    } else {
        return FCMMaterial :: give(aProperty, gp);
    }
}

int
ConcreteFCM :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    return FCMMaterial :: giveIPValue(answer, gp, type, tStep);
}



MaterialStatus *
ConcreteFCM :: giveStatus(GaussPoint *gp) const
{
    MaterialStatus *status = static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
    if ( status == NULL ) {
        // create a new one
        status = this->CreateStatus(gp);

        if ( status != NULL ) {
            gp->setMaterialStatus( status, this->giveNumber() );
            this->_generateStatusVariables(gp);
        }
    }

    return status;
}





///////////////////////////////////////////////////////////////////
//                      CONCRETE FCM STATUS                     ///
///////////////////////////////////////////////////////////////////

ConcreteFCMStatus :: ConcreteFCMStatus(int n, Domain *d, GaussPoint *gp) :
    FCMMaterialStatus(n, d, gp), RandomMaterialStatusExtensionInterface()
{}


ConcreteFCMStatus :: ~ConcreteFCMStatus()
{ }



void
ConcreteFCMStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    FCMMaterialStatus :: printOutputAt(file, tStep);
}





void
ConcreteFCMStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
//
{
    FCMMaterialStatus :: initTempStatus();
}



void
ConcreteFCMStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables correspond to newly reched equilibrium.
//
{
    FCMMaterialStatus :: updateYourself(tStep);
}

Interface *
ConcreteFCMStatus :: giveInterface(InterfaceType type)
{
    if ( type == RandomMaterialStatusExtensionInterfaceType ) {
        return static_cast< RandomMaterialStatusExtensionInterface * >(this);
    } else {
        return NULL;
    }
}



contextIOResultType
ConcreteFCMStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = FCMMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

contextIOResultType
ConcreteFCMStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = FCMMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK; // return succes
}
} // end namespace oofem
