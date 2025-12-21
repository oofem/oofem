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

#include "concretefcm.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"
#include "contextioerr.h"

namespace oofem {
REGISTER_Material(ConcreteFCM);

ConcreteFCM :: ConcreteFCM(int n, Domain *d) : FCMMaterial(n, d), RandomMaterialExtensionInterface()
{}

void
ConcreteFCM :: initializeFrom(InputRecord &ir)
{
    FCMMaterial :: initializeFrom(ir);
    RandomMaterialExtensionInterface :: initializeFrom(ir);

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
            throw ValueInputException(ir, _IFT_ConcreteFCM_soft_w, "The first entry in soft_w must be 0.");
        }

        if ( soft_w.giveSize() != soft_function_w.giveSize() ) {
            throw ValueInputException(ir, _IFT_ConcreteFCM_soft_function_w, "the size of 'soft_w' and 'soft(w)' must be the same");
        }
        break;

    case 5:
        softType = ST_LinearHardeningStrain;
        IR_GIVE_FIELD(ir, H, _IFT_ConcreteFCM_H);
        this->eps_f = 1.;
        IR_GIVE_OPTIONAL_FIELD(ir, eps_f, _IFT_ConcreteFCM_eps_f);
        IR_GIVE_FIELD(ir, Ft, _IFT_ConcreteFCM_ft);

        if ( this->giveCrackSpacing() > 0. ) {
            throw ValueInputException(ir, _IFT_FCM_crackSpacing, "crack spacing must not be defined in strain-defined materials (does not make sense)");
        }

        break;

    case 6:
        softType = ST_UserDefinedStrain;
        IR_GIVE_FIELD(ir, Ft, _IFT_ConcreteFCM_ft);
        IR_GIVE_FIELD(ir, soft_eps, _IFT_ConcreteFCM_soft_eps);
        IR_GIVE_FIELD(ir, soft_function_eps, _IFT_ConcreteFCM_soft_function_eps);

        if ( soft_eps.at(1) > 0. ) {
            throw ValueInputException(ir, _IFT_ConcreteFCM_soft_eps, "The first entry in soft_eps must be 0.");
        }

        if ( soft_eps.giveSize() != soft_function_eps.giveSize() ) {
            throw ValueInputException(ir, _IFT_ConcreteFCM_soft_function_eps, "the size of 'soft_eps' and 'soft(eps)' must be the same");
        }

        if ( this->giveCrackSpacing() > 0. ) {
            throw ValueInputException(ir, "foo", "crack spacing must not be defined in strain-defined materials (does not make sense)");
        }

        break;

    default:
        throw ValueInputException(ir, _IFT_ConcreteFCM_softType, "Unknown softening type");
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
            throw ValueInputException(ir, _IFT_ConcreteFCM_beta, "Parameter beta must be <= 1");
        }

        break;

    case 2:
        // todo: resolve problems with convergence when crack is wide open
        // this does not happen when beta is used instead sf
        shearType = SHR_Const_ShearFactorCoeff;
        sf = 20.;
        IR_GIVE_OPTIONAL_FIELD(ir, sf, _IFT_ConcreteFCM_sf);

        sf_numer = sf;
        IR_GIVE_OPTIONAL_FIELD(ir, sf_numer, _IFT_ConcreteFCM_sf_numer);
        break;

    case 3: // w is the maximum reached crack opening, not the current value!
        shearType = SHR_UserDefined_ShearRetFactor;
        IR_GIVE_FIELD(ir, beta_w, _IFT_ConcreteFCM_beta_w);
        IR_GIVE_FIELD(ir, beta_function, _IFT_ConcreteFCM_beta_function);

        if ( beta_w.giveSize() != beta_function.giveSize() ) {
            throw ValueInputException(ir, _IFT_ConcreteFCM_beta_function, "the size of 'beta_w' and 'beta(w)' must be the same");
        }
        break;


    default:
        throw ValueInputException(ir, _IFT_ConcreteFCM_shearType, "Shear stiffness reduction typ is unknown");
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

    case 3:
        shearStrengthType = SHS_Residual_Ft;
        break;

    default:
        throw ValueInputException(ir, _IFT_ConcreteFCM_shearStrengthType, "Shear stiffness reduction type is unknown");
    }
}


double
ConcreteFCM :: giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep, int i) const
//
// returns current cracking modulus according to crackStrain for i-th
// crackplane
//
{
    ConcreteFCMStatus *status = static_cast< ConcreteFCMStatus * >( this->giveStatus(gp) );
    double E = this->computeOverallElasticStiffness(gp, tStep);

    if ( status->giveTempCrackStatus(i) == pscm_NONE ) {
        return E * fcm_BIGNUMBER;
    }

    // very high stiffness also for compression
    if ( status->giveTempCrackStatus(i) == pscm_CLOSED ) {
        return E * fcm_BIGNUMBER;
    }

    return giveCrackingModulusInTension(rMode, gp, tStep, i);

}

  
double
ConcreteFCM :: giveCrackingModulusInTension(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep, int i) const
{

  ConcreteFCMStatus *status = static_cast< ConcreteFCMStatus * >( this->giveStatus(gp) );
  
  double Cf = 0.;
  double ef, emax, c1, c2, Le, Ncr, w;
  
  double ft = this->giveTensileStrength(gp, tStep);
  double Gf = this->giveFractureEnergy(gp, tStep);
  double ec =  status->giveTempCrackStrain(i);
  double E = this->computeOverallElasticStiffness(gp, tStep);
    
    
  
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
            ef = Gf / ( ft * Le );

            if ( emax == 0. ) { // just initiated and closing crack
                Cf = -ft / ef;
            } else { // further softening - negative stiffness
                Cf = -ft / ef *exp(-ec / ef);
            }
        } else if ( softType == ST_Linear ) {
            // fracturing strain at zero traction
            ef = 2. * Gf / ( ft * Le );

            if ( ( ec >= ef )  || ( emax >= ef ) ) {
                // fully open crack - no stiffness
                Cf = 0.;
            } else { // further softening - negative stiffness
                Cf = -ft / ef;
            }
        } else if ( softType == ST_Hordijk ) {
            c1 = 3.;
            c2 = 6.93;

            ef = 5.14 * Gf / ft / Le;

            if ( ( ec >= ef )  || ( emax >= ef ) ) {
                Cf = 0.;
            }  else if ( emax == 0. ) { // just initiated and closing crack
                Cf = ft * ( -c2 / ef - ( exp(-c2) * ( c1 * c1 * c1 + 1. ) ) / ef );
            } else { // further softening - negative stiffness
                Cf = ft * ( ( 3. * c1 * c1 * c1 * ec * ec * exp(-( c2 * ec ) / ef) ) /  ( ef * ef * ef ) - ( c2 * exp(-( c2 * ec ) / ef) * ( c1 * c1 * c1 * ec * ec * ec + ef * ef * ef ) ) / ( ef * ef * ef * ef ) - ( exp(-c2) * ( c1 * c1 * c1 + 1. ) ) / ef );
            }
        } else if ( softType == ST_UserDefinedCrack ) {
            w = ec * Le;

            if ( w > soft_w.at(soft_w.giveSize() ) ) {
                Cf = 0.;
#ifdef DEBUG
                OOFEM_WARNING("Crack width is larger than the last user-defined value in the traction-opening law.");
#endif
            } else if ( emax == 0. ) {
                Cf = ft * ( soft_function_w.at(2) - soft_function_w.at(1) ) / ( ( soft_w.at(2) - soft_w.at(1) ) / Le );
            } else { // softening
                for ( int i = 1; i <= soft_w.giveSize(); i++ ) {
                    if ( ( w - soft_w.at(i) ) < fcm_SMALL_STRAIN ) {
                        if ( i == 1 ) {
                            Cf =  ft * ( soft_function_w.at(2) - soft_function_w.at(1) ) /  ( ( soft_w.at(2) - soft_w.at(1) ) / Le );
                            break;
                        } else {
                            Cf =  ft * ( soft_function_w.at(i) - soft_function_w.at(i - 1) ) /  ( ( soft_w.at(i) - soft_w.at(i - 1) ) / Le );
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
		
		if ( ec > soft_eps.at(soft_eps.giveSize() ) ) {
                Cf = 0.;
#ifdef DEBUG
                OOFEM_WARNING("Crack strain is larger than the last user-defined value in the stress-crack-strain law.");
#endif
            } else if ( emax == 0. ) {
                Cf = ft * ( soft_function_eps.at(2) - soft_function_eps.at(1) ) / ( soft_eps.at(2) - soft_eps.at(1) );
            } else { // softening
                for ( int i = 1; i <= soft_eps.giveSize(); i++ ) {
                    if ( ( ec - soft_eps.at(i) ) < fcm_SMALL_STRAIN ) {
                        if ( i == 1 ) {
                            Cf =  ft * ( soft_function_eps.at(2) - soft_function_eps.at(1) ) /  ( soft_eps.at(2) - soft_eps.at(1) );
                            break;
                        } else {
                            Cf =  ft * ( soft_function_eps.at(i) - soft_function_eps.at(i - 1) ) /  ( soft_eps.at(i) - soft_eps.at(i - 1) );
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

        Cf = ConcreteFCM :: giveNormalCrackingStress(gp, tStep, emax * Ncr, i); // gets stress
        Cf /= ( emax ); // converted into stiffness
    }

    // to take into account crack-spacing
    return Cf / Ncr;
}


double
ConcreteFCM :: giveNormalCrackingStress(GaussPoint *gp, TimeStep *tStep, double ec, int i) const
//
// returns receivers Normal Stress in crack i  for given cracking strain
//
{
    ConcreteFCMStatus *status = static_cast< ConcreteFCMStatus * >( this->giveStatus(gp) );

    double traction;
    double Le, ef, emax, w, E, c1, c2;
    double Gf = this->giveFractureEnergy(gp, tStep);
    double ft = this->giveTensileStrength(gp, tStep);

    if ( status->giveTempCrackStatus(i) == pscm_CLOSED ) {
      E = this->computeOverallElasticStiffness(gp, tStep);
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
        ef = Gf / ( ft * Le );
        if ( ec >= emax ) {
            // further softening
            traction = ft * exp(-ec / ef);
        } else {
            // crack closing
            // or unloading or reloading regime
            traction = ft * ec / emax *exp(-emax / ef);
        }
    } else if ( softType  == ST_Linear ) {
        // fracturing strain at zero traction
        ef = 2. * Gf / ( ft * Le );

        if ( ( ec >= ef )  || ( emax >= ef ) ) {
            // fully open crack - no stiffness
            traction = 0.;
        } else if ( ec >= emax ) {
            // further softening
            traction = ft - ft * ec / ef;
        } else if ( ec <= 0. ) {
            traction = 0.;
        } else {
            traction = ft * ec * ( ef - emax ) / ( emax * ef );
        }
    } else if ( softType == ST_Hordijk ) {
        c1 = 3.;
        c2 = 6.93;

        ef = 5.14 * Gf / ft / Le;

        if ( ( ec >= ef )  || ( emax >= ef ) ) {
            // fully open crack - no stiffness
            traction = 0.;
        } else if ( ec >= emax ) {
            // further softening
            traction = ft * ( ( 1. + pow( ( c1 * ec / ef ), 3. ) ) * exp(-c2 * ec / ef) - ec / ef * ( 1. + c1 * c1 * c1 ) * exp(-c2) );
        } else if ( ec <= 0. ) {
            traction = 0.;
        } else {
            traction = ft * ec / emax * ( ( 1. + pow( ( c1 * emax / ef ), 3. ) ) * exp(-c2 * emax / ef) - emax / ef * ( 1. + c1 * c1 * c1 ) * exp(-c2) );
        }
    } else if ( softType == ST_UserDefinedCrack ) {
        w = ec * Le;

	traction = 0.;
	
        if ( w > soft_w.at(soft_w.giveSize() ) ) {
	  traction = 0.;
#ifdef DEBUG
            OOFEM_WARNING("Crack width is larger than the last user-defined value in the traction-opening law.");
#endif	    
        } else if ( ec >= emax ) { // softening
            for ( int i = 1; i <= soft_w.giveSize(); i++ ) {
                if ( ( w - soft_w.at(i) ) < fcm_SMALL_STRAIN ) {
                    if ( i == 1 ) { // eps \approx 0 bound
                        traction = soft_function_w.at(i) * ft;
                        break;
                    } else {
                        traction = soft_function_w.at(i - 1) +  ( soft_function_w.at(i) - soft_function_w.at(i - 1) ) / ( soft_w.at(i) - soft_w.at(i - 1) ) * ( w - soft_w.at(i - 1) );
                        traction *= ft;
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
                        traction = soft_function_w.at(i) * ft * ec / emax;
                        break;
                    } else {
                        traction = soft_function_w.at(i - 1) +  ( soft_function_w.at(i) - soft_function_w.at(i - 1) ) / ( soft_w.at(i) - soft_w.at(i - 1) ) * ( w - soft_w.at(i - 1) );
                        traction *= ft * ec / emax;
                        break;
                    }
                }
            }
        }
    } else if ( softType == ST_LinearHardeningStrain ) {
        if ( ( ec >= this->eps_f )  || ( emax >= this->eps_f ) ) {
            traction = 0.;
        } else if ( emax == 0. ) {
            traction = ft;
        } else {
            // unloading or reloading regime
            // emax * Le = crack width
            //      traction = ( this->Ft + emax * Le * this->H ) * ec / emax;
            // written in terms of strain and not crack width
            traction = ( ft + emax * this->H ) * ec / emax;
        }
    } else if ( softType == ST_UserDefinedStrain ) {

	    traction = 0.;
	    
	    if ( emax > soft_eps.at(soft_eps.giveSize() ) ) {
	      traction = 0.;
#ifdef DEBUG
	      OOFEM_WARNING("Cracking strain is larger than the last user-defined value in the traction-cracking-strain law.");
#endif
	    } else if ( ec >= emax ) { // softening
	      for ( int i = 1; i <= soft_eps.giveSize(); i++ ) {
                if ( ( ec - soft_eps.at(i) ) < fcm_SMALL_STRAIN ) {
		  if ( i == 1 ) { // eps \approx 0 bound
		    traction = soft_function_eps.at(i) * ft;
		    break;
		  } else {
                        traction = soft_function_eps.at(i - 1) +  ( soft_function_eps.at(i) - soft_function_eps.at(i - 1) ) / ( soft_eps.at(i) - soft_eps.at(i - 1) ) * ( ec - soft_eps.at(i - 1) );
                        traction *= ft;
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
                        traction = soft_function_eps.at(i) * ft * ec / emax;
                        break;
                    } else {
                        traction = soft_function_eps.at(i - 1) +  ( soft_function_eps.at(i) - soft_function_eps.at(i - 1) ) / ( soft_eps.at(i) - soft_eps.at(i - 1) ) * ( emax - soft_eps.at(i - 1) );
                        traction *= ft * ec / emax;
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
ConcreteFCM :: computeEffectiveShearModulus(GaussPoint *gp, TimeStep *tStep, int shearDirection) const

{
    double G, Geff, D2tot, N;
    int crackA, crackB;

    G = this->computeOverallElasticShearModulus(gp, tStep);

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
                }

                // total number of parallel cracks in both directions
                N = this->giveNumberOfCracksInDirection(gp, crackA) + this->giveNumberOfCracksInDirection(gp, crackB);
            } else {
                N = this->giveNumberOfCracksForShearDirection(gp, shearDirection);
            }

            Geff = G * this->beta / ( N - this->beta * ( N - 1 ) ); // for N = 1... Geff = G * beta
        } else if ( this->shearType == SHR_Const_ShearFactorCoeff ) {
            // contributions from shear terms in: De - De (De + Dcr)^-1 De
	  D2tot = this->computeNumerD2Modulus(gp, tStep, shearDirection);
            // number of parallel cracks already taken into account
            Geff = G * D2tot / ( G + D2tot );
        } else if ( this->shearType == SHR_UserDefined_ShearRetFactor ) {
            // contributions from shear terms in: De - De (De + Dcr)^-1 De
	  D2tot = this->computeNumerD2Modulus(gp, tStep, shearDirection);
            // number of parallel cracks already taken into account
            Geff = G * D2tot / ( G + D2tot );
        } else {
            OOFEM_ERROR("Unknown Shear Mode");
        }
    }

    return Geff;
}


double
ConcreteFCM :: computeD2ModulusForCrack(GaussPoint *gp, TimeStep *tStep, int icrack) const
{
    double D2, N, G, E, w, beta_m;

    if ( icrack >= 4 ) {
        OOFEM_ERROR("Unexpected crack number");
    }

    N = this->giveNumberOfCracksInDirection(gp, icrack);

    if ( this->isIntact(gp, icrack) || ( this->shearType == SHR_NONE ) ) {
      E = this->computeOverallElasticStiffness(gp, tStep);
        D2 = E * fcm_BIGNUMBER;
    } else {
        if ( shearType == SHR_Const_ShearRetFactor ) {
	  G = this->computeOverallElasticShearModulus(gp, tStep);
            D2 = G * this->beta / ( 1. - this->beta );
            D2 /= N;
        } else if ( shearType == SHR_Const_ShearFactorCoeff ) {
            // for the number of parallel cracks is taken care of in "giveCrackingModulus"
            // however, it is essential to reduce the final modulus due to the presence of multiple cracks
	  D2 = this->sf * ConcreteFCM::giveCrackingModulusInTension(SecantStiffness, gp, tStep, icrack);
            D2 /= N;
        } else if ( this->shearType == SHR_UserDefined_ShearRetFactor ) {
            // for the number of parallel cracks is taken care of in the evaluation of normalcrackopening
            // however, it is essential to reduce the final modulus due to the presence of multiple cracks

	  G = this->computeOverallElasticShearModulus(gp, tStep);
	  w = max( this->computeNormalCrackOpening(gp, icrack), this->computeMaxNormalCrackOpening(gp, tStep, icrack) );

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



double
ConcreteFCM :: computeNumerD2ModulusForCrack(GaussPoint *gp, TimeStep *tStep, int icrack) const
{
    double D2, N, E;

    if ( icrack >= 4 ) {
        OOFEM_ERROR("Unexpected crack number");
    }

    N = this->giveNumberOfCracksInDirection(gp, icrack);

    if ( this->isIntact(gp, icrack) || ( this->shearType == SHR_NONE ) ) {
      E = this->computeOverallElasticStiffness(gp, tStep);
        D2 = E * fcm_BIGNUMBER;
    } else {
        if ( shearType == SHR_Const_ShearRetFactor ) {
	  D2 = ConcreteFCM :: computeD2ModulusForCrack(gp, tStep, icrack);
        } else if ( shearType == SHR_Const_ShearFactorCoeff ) {
            // for the number of parallel cracks is taken care of in "giveCrackingModulus"
            // however, it is essential to reduce the final modulus due to the presence of multiple cracks
	  D2 = this->sf_numer * ConcreteFCM::giveCrackingModulusInTension(SecantStiffness, gp, tStep, icrack);
	  D2 /= N;
        } else if ( this->shearType == SHR_UserDefined_ShearRetFactor ) {
	  D2 = ConcreteFCM :: computeD2ModulusForCrack(gp, tStep, icrack);
        } else {
            D2 = 0;
            OOFEM_ERROR("Unknown Softening Mode");
        }
    }

    return D2;
}  

  

void
ConcreteFCM :: checkSnapBack(GaussPoint *gp, TimeStep *tStep, int i) const
{
    ConcreteFCMStatus *status = static_cast< ConcreteFCMStatus * >( this->giveStatus(gp) );

    double E, Gf, ft, Le;
    double ef = 0.;
    double slope;
    double efail, c1, c2;

    Le = status->giveCharLength(i);

    Gf = this->giveFractureEnergy(gp, tStep);
    ft = this->giveTensileStrength(gp, tStep);

    E = this->computeOverallElasticStiffness(gp, tStep);

    if ( softType == ST_NONE ) {
        OOFEM_ERROR("For unknown reason the tensile strength has been exceeded and cracking has been initiated!");
    } else if ( softType == ST_Exponential ) {
        ef = Gf / ( ft * Le );
    } else if ( softType == ST_Linear ) {
        ef = 2. * Gf / ( ft * Le );
    } else if ( softType == ST_Hordijk ) {
        efail =  5.14 * Gf / ft / Le;
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

    if ( ef <= ft / E ) {
        OOFEM_WARNING("ef %e < e0 %e, this leads to material snapback in element %d, characteristic length %f (E=%e, ft=%e)", ef, ft / E, gp->giveElement()->giveGlobalNumber(), Le, E, ft);
    }
}


double
ConcreteFCM :: maxShearStress(GaussPoint *gp, TimeStep *tStep, int i) const {
    int dir_1, dir_2;
    double maxTau, crackOpening, wmax, scale;

    if ( this->isIntactForShear(gp, i) ) {
      return this->giveTensileStrength(gp, tStep) * fcm_BIGNUMBER;
    }

    if ( this->shearStrengthType == SHS_NONE ) {
        maxTau = this->giveTensileStrength(gp, tStep)  * fcm_BIGNUMBER;
    } else if ( this->shearStrengthType == SHS_Const_Ft ) {
        maxTau = this->giveTensileStrength(gp, tStep);
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
        }

        // from temporary strain
        crackOpening = max( this->computeNormalCrackOpening(gp, dir_1), this->computeNormalCrackOpening(gp, dir_2) );

        // from the entire history
        wmax = max( this->computeMaxNormalCrackOpening(gp, tStep, dir_1), this->computeMaxNormalCrackOpening(gp, tStep, dir_2) );

        crackOpening = max(wmax, crackOpening);

        scale = 1000. / lengthScale;
        maxTau = 0.18 * sqrt(this->fc) / ( 0.31 + 24. * crackOpening * scale / ( this->ag * scale + 16. ) );

    } else if ( this->shearStrengthType == SHS_Residual_Ft ) {
        maxTau = this->computeResidualTensileStrength(gp, tStep);

    } else {
        OOFEM_ERROR("Unexpected shearStrengthType");
    }

    return maxTau;
}

double
ConcreteFCM :: computeResidualTensileStrength(GaussPoint *gp, TimeStep *tStep) const {

  ConcreteFCMStatus *status = static_cast< ConcreteFCMStatus * >( this->giveStatus(gp) );
  
  double sigma;
  int nCracks;
  double emax;

  nCracks = status->giveNumberOfTempCracks();
  
  sigma = this->giveTensileStrength(gp, tStep);
  
  for ( int i = 1; i <= nCracks; i++ ) {
    
    emax = status->giveMaxCrackStrain(i);
    
    if (emax > 0.) {
      emax /= this->giveNumberOfCracksInDirection(gp, i);
      
      sigma = ConcreteFCM :: giveNormalCrackingStress(gp, tStep, emax, i);
      sigma = min( sigma, this->giveTensileStrength(gp, tStep) );
    }
  }

  return sigma;
}


double
ConcreteFCM :: give(int aProperty, GaussPoint *gp) const
{
    this->giveStatus(gp);
    
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

  if ( type == IST_TensileStrength ) {
    answer.resize(1);
    answer.at(1) = this->giveTensileStrength(gp, tStep);   
    return 1;
    
  } else if ( type == IST_ResidualTensileStrength ) {

    answer.resize(1);
    answer.zero();
    answer.at(1) = this->computeResidualTensileStrength(gp, tStep);
    
    return 1;
  }

    return FCMMaterial :: giveIPValue(answer, gp, type, tStep);
}


MaterialStatus *
ConcreteFCM :: giveStatus(GaussPoint *gp) const
{
    if (gp->hasMaterialStatus()) {
        return static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
    } else {
        MaterialStatus *status = static_cast<MaterialStatus*> (gp->setMaterialStatus (this->CreateStatus(gp)));
        this->_generateStatusVariables(gp);
        return status;
    }
}




///////////////////////////////////////////////////////////////////
//                      CONCRETE FCM STATUS                     ///
///////////////////////////////////////////////////////////////////

ConcreteFCMStatus :: ConcreteFCMStatus(GaussPoint *gp) :
    FCMMaterialStatus(gp), RandomMaterialStatusExtensionInterface()
{}



void
ConcreteFCMStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
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
        return nullptr;
    }
}


void
ConcreteFCMStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    FCMMaterialStatus :: saveContext(stream, mode);
}


void
ConcreteFCMStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    FCMMaterialStatus :: restoreContext(stream, mode);
}

} // end namespace oofem
