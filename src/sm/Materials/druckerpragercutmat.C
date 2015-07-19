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

#include "druckerpragercutmat.h"
#include "Materials/isolinearelasticmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(DruckerPragerCutMat);

// constructor
DruckerPragerCutMat :: DruckerPragerCutMat(int n, Domain *d) : MPlasticMaterial2(n, d)
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    this->nsurf = 4;
    this->rmType = mpm_CuttingPlane;
    //     this->rmType = mpm_ClosestPoint;

    this->plType = nonassociatedPT; // Rankine associated, DP nonassociated

    sigT = 0.;
    H = 0.;
    omegaCrit = 0.;
    a = 0.;
    yieldTol = 0.;
    newtonIter = 30;
    G = 0.;
    K = 0.;
}

// destructor
DruckerPragerCutMat :: ~DruckerPragerCutMat()
{ }

// specifies whether a given material mode is supported by this model
int
DruckerPragerCutMat :: hasMaterialModeCapability(MaterialMode mode)
{
    return mode == _3dMat || mode == _PlaneStrain;
}

// reads the model parameters from the input file
IRResultType
DruckerPragerCutMat :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // required by IR_GIVE_FIELD macro

    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;
    result = linearElasticMaterial->initializeFrom(ir); // takes care of elastic constants
    if ( result != IRRT_OK ) return result;

    G = static_cast< IsotropicLinearElasticMaterial * >(linearElasticMaterial)->giveShearModulus();
    K = static_cast< IsotropicLinearElasticMaterial * >(linearElasticMaterial)->giveBulkModulus();

    IR_GIVE_FIELD(ir, tau0, _IFT_DruckerPragerCutMat_tau0); // initial yield stress under pure shear (DP model)
    IR_GIVE_FIELD(ir, sigT, _IFT_DruckerPragerCutMat_sigT); // uniaxial tensile strength for cut-off, (Rankine plasticity model)
    IR_GIVE_FIELD(ir, alpha, _IFT_DruckerPragerCutMat_alpha); // friction coefficient (DP model)
    alphaPsi = alpha;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaPsi, _IFT_DruckerPragerCutMat_alphapsi); //dilatancy coefficient (DP model)
    IR_GIVE_OPTIONAL_FIELD(ir, H, _IFT_DruckerPragerCutMat_h); // hardening modulus  (DP model)
    IR_GIVE_OPTIONAL_FIELD(ir, omegaCrit, _IFT_DruckerPragerCutMat_omegaCrit); // critical damage
    IR_GIVE_OPTIONAL_FIELD(ir, a, _IFT_DruckerPragerCutMat_a); // exponent in damage law
    IR_GIVE_OPTIONAL_FIELD(ir, yieldTol, _IFT_DruckerPragerCutMat_yieldTol); //tolerance of the error in the yield criterion
    IR_GIVE_OPTIONAL_FIELD(ir, newtonIter, _IFT_DruckerPragerCutMat_newtonIter); //Maximum number of iterations in lambda search

    return IRRT_OK;
}


MaterialStatus *
DruckerPragerCutMat :: CreateStatus(GaussPoint *gp) const
{
    return new MPlasticMaterial2Status(1, this->giveDomain(), gp, this->giveSizeOfReducedHardeningVarsVector(gp));
}

double
DruckerPragerCutMat :: computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector, const FloatArray &strainSpaceHardeningVariables)
{
    //strainSpaceHardeningVariables = kappa
    if ( isurf <= 3 ) { //Rankine, surfaces 1,2,3
        FloatArray princStress(3);
        this->computePrincipalValues(princStress, stressVector, principal_stress);
        return princStress.at(isurf) - this->sigT;
    } else { //Drucker-Prager, surface 4
        double volumetricStress;
        double DPYieldStressInShear = tau0 + H *strainSpaceHardeningVariables.at(4);
        double JTwo;
        FloatArray deviatoricStress;

        volumetricStress = this->computeDeviatoricVolumetricSplit(deviatoricStress, stressVector);
        JTwo = computeSecondStressInvariant(deviatoricStress);
        return 3. * alpha * volumetricStress + sqrt(JTwo) - DPYieldStressInShear;
    }
}

//associated and nonassociated flow rule
void
DruckerPragerCutMat :: computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector, const FloatArray &stressSpaceHardeningVars)
{
    FloatArray princStress(3);
    FloatMatrix t(3, 3);

    // compute principal stresses and their directions
    this->computePrincipalValDir(princStress, t, stressVector, principal_stress);

    // derivation through stress transformation matrix. The transformation matrix is stored in t columnwise
    answer.resize(6);
    if ( isurf <= 3 ) { //Rankine associated
        answer.at(1) = t.at(1, isurf) * t.at(1, isurf); //xx = 11
        answer.at(2) = t.at(2, isurf) * t.at(2, isurf); //yy = 22
        answer.at(3) = t.at(3, isurf) * t.at(3, isurf); //zz = 33
        answer.at(4) = t.at(2, isurf) * t.at(3, isurf); //yz = 23
        answer.at(5) = t.at(1, isurf) * t.at(3, isurf); //xz = 13
        answer.at(6) = t.at(1, isurf) * t.at(2, isurf); //xy = 12
    } else { //DP nonassociated
        double sqrtJTwo;
        FloatArray deviatoricStress;

        this->computeDeviatoricVolumetricSplit(deviatoricStress, stressVector);
        sqrtJTwo = sqrt( computeSecondStressInvariant(deviatoricStress) );

        answer.at(1) = alphaPsi + deviatoricStress.at(1) / 2. / sqrtJTwo;
        answer.at(2) = alphaPsi + deviatoricStress.at(2) / 2. / sqrtJTwo;
        answer.at(3) = alphaPsi + deviatoricStress.at(3) / 2. / sqrtJTwo;
        answer.at(4) = stressVector.at(4) / sqrtJTwo;
        answer.at(5) = stressVector.at(5) / sqrtJTwo;
        answer.at(6) = stressVector.at(6) / sqrtJTwo;
    }
}

//necesarry only for mpm_ClosestPoint, see Jirasek: Inelastic analysis of structures, pp. 411.
//Hessian matrix
void
DruckerPragerCutMat :: computeReducedSSGradientMatrix(FloatMatrix &gradientMatrix,  int isurf, GaussPoint *gp, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables)
{
    switch ( gp->giveMaterialMode() ) {
    case _3dMat:
        gradientMatrix.resize(6, 6);
        break;
    case _PlaneStrain:
        gradientMatrix.resize(4, 4);
        break;
    default:
        OOFEM_ERROR("Unknown material mode (%s)", __MaterialModeToString( gp->giveMaterialMode() ) );
    }

    gradientMatrix.zero();

    if ( isurf == 4 ) {
        double c1 = 0.;
        FloatArray deviatoricStress;
        double JTwo, sqrtJTwo;

        computeDeviatoricVolumetricSplit(deviatoricStress, fullStressVector);
        JTwo = computeSecondStressInvariant(deviatoricStress);
        sqrtJTwo = sqrt(JTwo);

        if ( gp->giveMaterialMode() == _3dMat ) {
            for ( int i = 1; i <= 6; i++ ) {
                for ( int j = i; j <= 6; j++ ) {
                    if ( ( i == 1 && j == 1 ) || ( i == 2 && j == 2 ) || ( i == 3 && j == 3 ) ) {
                        c1 = 2 / 3.;
                    } else if ( ( i == 4 && j == 4 ) || ( i == 5 && j == 5 ) || ( i == 6 && j == 6 ) ) {
                        c1 = 1.;
                    } else if ( i <= 3 && j <= 3 ) {
                        c1 = -1 / 3.;
                    } else {
                        c1 = 0;
                    }
                    gradientMatrix.at(i, j) = 0.5 / JTwo * ( c1 * sqrtJTwo - deviatoricStress.at(i) * deviatoricStress.at(j) / 2. / sqrtJTwo );
                }
            }
            gradientMatrix.symmetrized();
        } else if ( gp->giveMaterialMode() == _PlaneStrain ) {
            for ( int i = 1; i <= 4; i++ ) {
                for ( int j = i; j <= 4; j++ ) {
                    if ( ( i == 1 && j == 1 ) || ( i == 2 && j == 2 ) || ( i == 3 && j == 3 ) ) {
                        c1 = 2 / 3.;
                    } else if ( ( i == 4 && j == 4 ) ) {
                        c1 = 1.;
                    } else if ( i <= 3 && j <= 3 ) {
                        c1 = -1 / 3.;
                    } else {
                        c1 = 0;
                    }
                    gradientMatrix.at(i, j) = 0.5 / JTwo * ( c1 * sqrtJTwo - deviatoricStress.at(i) * deviatoricStress.at(j) / 2. / sqrtJTwo );
                }
            }
            gradientMatrix.symmetrized();
        } else {
            OOFEM_ERROR("Unknown material mode (%s)", __MaterialModeToString( gp->giveMaterialMode() ) );
        }
    }
}


void
DruckerPragerCutMat :: computeReducedElasticModuli(FloatMatrix &answer,
                                                   GaussPoint *gp,
                                                   TimeStep *tStep)
{  /* Returns elastic moduli in reduced stress-strain space*/
    this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, ElasticStiffness, gp, tStep);
}

//answer is dkappa (cumulative plastic strain), flow rule
void DruckerPragerCutMat :: computeStrainHardeningVarsIncrement(FloatArray &answer, GaussPoint *gp, const FloatArray &stress, const FloatArray &dlambda, const FloatArray &dplasticStrain, const IntArray &activeConditionMap)
{
    answer.resize(4);
    answer.zero();

    if ( dlambda.at(4) > 0. ) {
        answer.at(4) = dlambda.at(4) * sqrt(1. / 3. + 2 * alphaPsi * alphaPsi);
    }
}

// Computes the derivative of yield/loading function with respect to kappa_1, kappa_2 etc.
void DruckerPragerCutMat :: computeKGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables)
{
    answer.resize(1);  //1 hardening variable for DP model - kappa
    answer.zero();

    if ( isurf == 4 ) {
        answer.at(1) = this->H;
    }
}

//necesarry only for mpm_ClosestPoint
//Computes second mixed derivative of loading function with respect to stress and hardening vars.
void DruckerPragerCutMat :: computeReducedSKGradientMatrix(FloatMatrix &gradientMatrix, int isurf, GaussPoint *gp, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables)
{
    int size = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );
    gradientMatrix.resize(size, 1); //six stresses in 3D and one kappa
    gradientMatrix.zero();
}

// computes dKappa_i/dsig_j gradient matrix
void DruckerPragerCutMat :: computeReducedHardeningVarsSigmaGradient(FloatMatrix &answer, GaussPoint *gp, const IntArray &activeConditionMap, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVars, const FloatArray &dlambda)
{
    int size = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );
    answer.resize(1, size);
    answer.zero();
}

// computes dKappa_i/dLambda_j for one surface
void DruckerPragerCutMat :: computeReducedHardeningVarsLamGradient(FloatMatrix &answer, GaussPoint *gp, int actSurf, const IntArray &activeConditionMap, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVars, const FloatArray &dlambda)
{
    int indx;
    answer.resize(1, actSurf); //actSurf = number of active surfaces
    answer.zero();

    if ( ( indx = activeConditionMap.at(4) ) ) {
        if ( dlambda.at(4) > 0. ) {
            answer.at(1, indx) = sqrt(1. / 3. + 2 * alphaPsi * alphaPsi);
        }
    }
}


int
DruckerPragerCutMat :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    //MaterialStatus *status = this->giveStatus(gp);
    if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.zero();
        ///@todo Actually export the relevant damage value here!
        //answer.at(1) = status->giveDamage();
        return 1;
    } else if ( type == IST_DamageTensor ) {
        answer.resize(6);
        answer.zero();
        //answer.at(1) = answer.at(2) = answer.at(3) = status->giveDamage();
        return 1;
    } else {
        return MPlasticMaterial2 :: giveIPValue(answer, gp, type, tStep);
    }
}
} // end namespace oofem
