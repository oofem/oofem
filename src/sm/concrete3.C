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

#include "concrete3.h"
#include "gausspnt.h"
#include "mathfem.h"

namespace oofem {
Concrete3 :: Concrete3(int n, Domain *d) : RCM2Material(n, d)
    //
    // constructor
    //
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
}




MaterialStatus *
Concrete3 :: CreateStatus(GaussPoint *gp) const
/*
 * creates new  material status  corresponding to this class
 */
{
    // double Ee, Gf, Et, Etp, Ft, Le, minEffStrainForFullyOpenCrack;
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) RCM2Material :: CreateStatus(gp);
    return status;
}


int
Concrete3 :: checkSizeLimit(GaussPoint *gp, double charLength)
//
// checks if element size (charLength) is too big
// so that tension strength must be reduced followed
// by sudden stress drop
//
{
    double Ee, Gf, Ft, LeCrit;

    Ee = this->give(pscm_Ee, gp);
    Gf = this->give(pscm_Gf, gp);
    Ft = this->give(pscm_Ft, gp);

    LeCrit = 2.0 * Gf * Ee / ( Ft * Ft );
    return ( charLength < LeCrit );
}

double
Concrete3 :: computeStrength(GaussPoint *gp, double charLength)
//
// computes strength for given gp,
// which may be reduced according to length of "fracture process zone"
// to be energetically correct
//
{
    double Ee, Gf, Ft;

    Ee = this->give(pscm_Ee, gp);
    Gf = this->give(pscm_Gf, gp);
    Ft = this->give(pscm_Ft, gp);

    if ( this->checkSizeLimit(gp, charLength) ) {
        ;
    } else {
        // we reduce Ft and there is no softening but sudden drop
        Ft = sqrt(2. * Ee * Gf / charLength);
        //
        printf("\n Reducing Ft to %f in element %d, gp %d, Le %f",
               Ft, gp->giveElement()->giveNumber(), gp->giveNumber(), charLength);
        //
    }

    return Ft;
}

/*
 * double
 * Concrete3 ::giveMinCrackStrainsForFullyOpenCrack (GaussPoint* gp, int i)
 * //
 * // computes MinCrackStrainsForFullyOpenCrack for given gp and i-th crack
 * //
 * {
 * RCM2MaterialStatus *status = (RCM2MaterialStatus*) this -> giveStatus (gp);
 * double Le, Gf, Ft;
 *
 * Le = status -> giveCharLength (i);
 * Gf = this->give(pscm_Gf);
 * Ft = this->computeStrength (gp, Le);
 *
 * return 2.0*Gf/(Le*Ft);
 * }
 */


double
Concrete3 :: giveMinCrackStrainsForFullyOpenCrack(GaussPoint *gp, int i)
//
// computes MinCrackStrainsForFullyOpenCrack for given gp and i-th crack
//
{
    if ( softeningMode == linearSoftening ) {
        RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(gp);
        double Le, Gf, Ft;

        Le = status->giveCharLength(i);
        Gf = this->give(pscm_Gf, gp);
        Ft = this->computeStrength(gp, Le);

        return 2.0 * Gf / ( Le * Ft );
    } else {
        //return Gf/(Le*Ft); // Exponential softening
        return 1.e6; // Exponential softening
    }
}


/*
 * void
 * Concrete3 :: updateStatusForNewCrack (GaussPoint* gp, int i, double Le)
 * //
 * // updates gp status when new crack-plane i is formed with charLength Le
 * // updates Le and computes and sets minEffStrainForFullyOpenCrack
 * //
 * {
 * RCM2MaterialStatus *status = (RCM2MaterialStatus*) this -> giveStatus (gp);
 *
 * if (Le <= 0) {
 * char errMsg [80];
 * sprintf (errMsg,"Element %d returned zero char length",
 *    gp->giveElement()->giveNumber());
 * _error (errMsg);
 * }
 *
 * status -> setCharLength(i, Le);
 * status ->
 * setMinCrackStrainsForFullyOpenCrack (i, this->giveMinCrackStrainsForFullyOpenCrack(gp,i));
 * }
 */

/*
 * double
 * Concrete3 :: giveCrackingModulus (MatResponseMode rMode, GaussPoint* gp,
 *          double crackStrain,int i)
 * //
 * // returns current cracking modulus according to crackStrain for i-th
 * // crackplane
 * // now linear softening is implemented
 * // see also CreateStatus () function.
 * // softening modulus represents a relation between the normal crack strain
 * // rate and the normal stress rate.
 * //
 * {
 * double Ee, Gf, Et, Etp, Cf, Ft, Le, minEffStrainForFullyOpenCrack;
 * RCM2MaterialStatus *status = (RCM2MaterialStatus*) this -> giveStatus (gp);
 *
 * //
 * // now we have to set proper reduced strength and softening modulus Et
 * // in order to obtain energetically correct solution using the concept
 * // of fracture energy Gf, which is a material constant.
 * //
 * Ee = this->give(pscm_Ee);
 * Gf = this->give(pscm_Gf);
 * Ft = this->give(pscm_Ft);
 * Le = status->giveCharLength(i);
 *
 * Ft = this->computeStrength (gp, Le);
 * minEffStrainForFullyOpenCrack = this->giveMinCrackStrainsForFullyOpenCrack(gp,i);
 *
 * if (rMode == TangentStiffness) {
 * if (this->checkSizeLimit(gp, Le)) {
 * if ((crackStrain >= minEffStrainForFullyOpenCrack) ||
 *   (status->giveTempMaxCrackStrain(i) >= minEffStrainForFullyOpenCrack)) {
 *  // fully open crack - no stiffness
 *  Cf = 0.;
 * } else if (crackStrain >= status->giveTempMaxCrackStrain(i)) {
 *  // further softening
 *  Cf = -Ft/minEffStrainForFullyOpenCrack;
 * } else {
 *  // unloading or reloading regime
 *  Cf = Ft*(minEffStrainForFullyOpenCrack - status->giveTempMaxCrackStrain(i)) /
 *   (status->giveTempMaxCrackStrain(i) * minEffStrainForFullyOpenCrack);
 * }
 * } else {
 * Cf = 0.;
 * }
 * } else {
 * if (this->checkSizeLimit(gp, Le)) {
 * if ((crackStrain >= minEffStrainForFullyOpenCrack) ||
 *   (status->giveTempMaxCrackStrain(i) >= minEffStrainForFullyOpenCrack)) {
 *  // fully open crack - no stiffness
 *  Cf = 0.;
 * } else  {
 *  // unloading or reloading regime
 *  Cf = Ft*(minEffStrainForFullyOpenCrack - status->giveTempMaxCrackStrain(i)) /
 *   (status->giveTempMaxCrackStrain(i) * minEffStrainForFullyOpenCrack);
 * }
 * } else {
 * Cf = 0.;
 * }
 *
 * }
 * return Cf;
 * }
 */

double
Concrete3 :: giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp,
                                 double crackStrain, int i)
//
// returns current cracking modulus according to crackStrain for i-th
// crackplane
// now linear softening is implemented
// see also CreateStatus () function.
// softening modulus represents a relation between the normal crack strain
// rate and the normal stress rate.
//
{
    //double Ee, Gf;
    double Cf, Gf, Ft, Le, ef, minEffStrainForFullyOpenCrack;
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(gp);

    //
    // now we have to set proper reduced strength and softening modulus Et
    // in order to obtain energetically correct solution using the concept
    // of fracture energy Gf, which is a material constant.
    //
    //Ee = this->give(pscm_Ee);
    Gf = this->give(pscm_Gf, gp);
    Ft = this->give(pscm_Ft, gp);
    Le = status->giveCharLength(i);

    Ft = this->computeStrength(gp, Le);
    minEffStrainForFullyOpenCrack = this->giveMinCrackStrainsForFullyOpenCrack(gp, i);

    if ( rMode == TangentStiffness ) {
        if ( this->checkSizeLimit(gp, Le) ) {
            if ( softeningMode == linearSoftening ) {
                if ( ( crackStrain >= minEffStrainForFullyOpenCrack ) ||
                    ( status->giveTempMaxCrackStrain(i) >= minEffStrainForFullyOpenCrack ) ) {
                    // fully open crack - no stiffness
                    Cf = 0.;
                } else if ( crackStrain >= status->giveTempMaxCrackStrain(i) ) {
                    // further softening
                    Cf = -Ft / minEffStrainForFullyOpenCrack;
                } else {
                    // unloading or reloading regime
                    Cf = Ft * ( minEffStrainForFullyOpenCrack - status->giveTempMaxCrackStrain(i) ) /
                         ( status->giveTempMaxCrackStrain(i) * minEffStrainForFullyOpenCrack );
                }
            } else { // exponential softening
                ef = Gf / ( Le * Ft );
                if ( crackStrain >= status->giveTempMaxCrackStrain(i) ) {
                    // further softening
                    Cf = -Ft / ef *exp(-crackStrain / ef);
                } else {
                    // unloading or reloading regime
                    Cf = Ft / status->giveTempMaxCrackStrain(i) * exp(-status->giveTempMaxCrackStrain(i) / ef);
                }
            }
        } else {
            Cf = 0.;
        }
    } else { // secant stiffness
        if ( this->checkSizeLimit(gp, Le) ) {
            if ( softeningMode == linearSoftening ) {
                if ( ( crackStrain >= minEffStrainForFullyOpenCrack ) ||
                    ( status->giveTempMaxCrackStrain(i) >= minEffStrainForFullyOpenCrack ) ) {
                    // fully open crack - no stiffness
                    Cf = 0.;
                } else {
                    // unloading or reloading regime
                    Cf = Ft * ( minEffStrainForFullyOpenCrack - status->giveTempMaxCrackStrain(i) ) /
                         ( status->giveTempMaxCrackStrain(i) * minEffStrainForFullyOpenCrack );
                }
            } else { // exponential softening
                ef = Gf / ( Le * Ft );
                Cf = Ft / status->giveTempMaxCrackStrain(i) * exp(-status->giveTempMaxCrackStrain(i) / ef);
            }
        } else {
            Cf = 0.;
        }
    }

    return Cf;
}

/*
 *
 * double
 * Concrete3 :: giveNormalCrackingStress (GaussPoint*gp, double crackStrain, int i)
 * //
 * // returns receivers Normal Stress in crack i  for given cracking strain
 * //
 * {
 * double Cf, Ft, Le, answer, minEffStrainForFullyOpenCrack;
 * RCM2MaterialStatus *status = (RCM2MaterialStatus*) this->giveStatus (gp);
 * minEffStrainForFullyOpenCrack = this->giveMinCrackStrainsForFullyOpenCrack(gp,i);
 *
 * Cf = this -> giveCrackingModulus (TangentStiffness,gp,crackStrain,i);   // < 0
 * Le = status->giveCharLength(i);
 * Ft = this->computeStrength (gp, Le);
 *
 * if (this->checkSizeLimit(gp,Le)) {
 *
 * if ((crackStrain >= minEffStrainForFullyOpenCrack) ||
 *  (status->giveTempMaxCrackStrain(i) >= minEffStrainForFullyOpenCrack)) {
 * // fully open crack - no stiffness
 * answer = 0.;
 * } else if (crackStrain >= status->giveTempMaxCrackStrain(i)) {
 * // further softening
 * answer = Ft + Cf*crackStrain;
 * } else if (crackStrain <= 0.) {
 * // crack closing
 * // no stress due to cracking
 * answer = 0.;
 * } else {
 * // unloading or reloading regime
 * answer = crackStrain * Ft *
 *  (minEffStrainForFullyOpenCrack - status->giveTempMaxCrackStrain(i))/
 *   (status->giveTempMaxCrackStrain(i) * minEffStrainForFullyOpenCrack);
 * }
 * } else {
 * answer = 0.;
 * }
 *
 * return answer;
 * }
 */

double
Concrete3 :: giveNormalCrackingStress(GaussPoint *gp, double crackStrain, int i)
//
// returns receivers Normal Stress in crack i  for given cracking strain
//
{
    double Cf, Ft, Gf, Le, answer, ef, minEffStrainForFullyOpenCrack;
    RCM2MaterialStatus *status = ( RCM2MaterialStatus * ) this->giveStatus(gp);
    minEffStrainForFullyOpenCrack = this->giveMinCrackStrainsForFullyOpenCrack(gp, i);

    Cf = this->giveCrackingModulus(TangentStiffness, gp, crackStrain, i); // < 0
    Le = status->giveCharLength(i);
    Ft = this->computeStrength(gp, Le);
    Gf = this->give(pscm_Gf, gp);

    if ( this->checkSizeLimit(gp, Le) ) {
        if ( softeningMode == linearSoftening ) {
            if ( ( crackStrain >= minEffStrainForFullyOpenCrack ) ||
                ( status->giveTempMaxCrackStrain(i) >= minEffStrainForFullyOpenCrack ) ) {
                // fully open crack - no stiffness
                answer = 0.;
            } else if ( crackStrain >= status->giveTempMaxCrackStrain(i) ) {
                // further softening
                answer = Ft + Cf * crackStrain;
            } else if ( crackStrain <= 0. ) {
                // crack closing
                // no stress due to cracking
                answer = 0.;
            } else {
                // unloading or reloading regime
                answer = crackStrain * Ft *
                         ( minEffStrainForFullyOpenCrack - status->giveTempMaxCrackStrain(i) ) /
                         ( status->giveTempMaxCrackStrain(i) * minEffStrainForFullyOpenCrack );
            }
        } else { // Exponential softening
            ef = Gf / ( Le * Ft );
            if ( crackStrain >= status->giveTempMaxCrackStrain(i) ) {
                // further softening
                answer = Ft * exp(-crackStrain / ef);
            } else {
                // crack closing
                // or unloading or reloading regime
                answer = Ft * crackStrain / status->giveTempMaxCrackStrain(i) *
                         exp(-status->giveTempMaxCrackStrain(i) / ef);
            }
        }
    } else {
        answer = 0.;
    }

    return answer;
}

/*
 * double
 * Concrete3 :: giveShearRetentionFactor (GaussPoint* gp, double eps_cr, int i)
 * //
 * // Returns shear retention factor, according to given crack strain
 * // for i-th crack in gp.
 * //
 * {
 * double answer;
 * RCM2MaterialStatus *status = (RCM2MaterialStatus*) this->giveStatus (gp);
 * double maxEps_cr = status->giveMinCrackStrainsForFullyOpenCrack()->at(i);
 *
 * // lack of experimental evidence - just poor constant
 * answer = this-> shearRetFactor;
 *
 * //
 * // answer = 1. * (maxEps_cr - eps_cr) / (maxEps_cr);
 * // // optional - check the result
 * // if ((answer > 1) || (answer < 0)) error ("giveShearRetentionFactor: consistency error");
 * //
 * return answer;
 * }
 */


IRResultType
Concrete3 :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    //double value ;
    int exmode;

    RCM2Material :: initializeFrom(ir);
    exmode = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, exmode, IFT_Concrete3_exp_soft, "exp_soft"); // Macro
    if ( exmode ) {
        softeningMode = exponentialSoftening;
    } else {
        softeningMode = linearSoftening;
    }

    /*
     * value = readDouble (initString,"shearretfactor");
     * shearRetFactor = value;
     * if (shearRetFactor < 0.001) shearRetFactor = 0.001;
     * if (shearRetFactor > 0.50) shearRetFactor = 0.50;
     */

    return IRRT_OK;
}
} // end namespace oofem
