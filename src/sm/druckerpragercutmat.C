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

#include "druckerpragercutmat.h"
#include "isolinearelasticmaterial.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "gausspnt.h"
#include "stressvector.h"
#include "strainvector.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"

namespace oofem {
// constructor
DruckerPragerCutMat :: DruckerPragerCutMat(int n, Domain *d) : MPlasticMaterial2(n, d)
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    this->nsurf = 4;
//     this->rmType = mpm_CuttingPlane;
    this->rmType = mpm_ClosestPoint;
    
    this->plType = nonassociatedPT;// Rankine associated, DP nonassociated
    
    sigT = 0.;
    H = 0.;
    omegaCrit = 0.;
    a=0.;
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
    if ( ( mode == _3dMat ) || ( mode == _PlaneStrain ) ) {
        return 1;
    }

    return 0;
}

// reads the model parameters from the input file
IRResultType
DruckerPragerCutMat :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // required by IR_GIVE_FIELD macro
    IRResultType result;                 // required by IR_GIVE_FIELD macro

    StructuralMaterial :: initializeFrom(ir);
    linearElasticMaterial->initializeFrom(ir); // takes care of elastic constants
    G = ( ( IsotropicLinearElasticMaterial * ) linearElasticMaterial )->giveShearModulus();
    K = ( ( IsotropicLinearElasticMaterial * ) linearElasticMaterial )->giveBulkModulus();

    IR_GIVE_FIELD(ir, tau0, IFT_DruckerPragerCutMat_tau0, "tau0"); // initial yield stress under pure shear (DP model)
    IR_GIVE_FIELD(ir, sigT, IFT_DruckerPragerCutMat_sigT, "sigt"); // uniaxial tensile strength for cut-off, (Rankine plasticity model)
    IR_GIVE_FIELD(ir, alpha, IFT_DruckerPragerCutMat_alpha, "alpha"); // friction coefficient (DP model)
    alphaPsi=alpha;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaPsi, IFT_DruckerPragerCutMat_alphapsi, "alphapsi"); //dilatancy coefficient (DP model)
    IR_GIVE_OPTIONAL_FIELD(ir, H, IFT_DruckerPragerCutMat_h, "h"); // hardening modulus  (DP model)
    IR_GIVE_OPTIONAL_FIELD(ir, omegaCrit, IFT_DruckerPragerCutMat_omegaCrit, "omega_crit"); // critical damage
    IR_GIVE_OPTIONAL_FIELD(ir, a, IFT_DruckerPragerCutMat_a, "a"); // exponent in damage law
    IR_GIVE_OPTIONAL_FIELD(ir, yieldTol, IFT_DruckerPragerCutMat_yieldTol, "yieldtol"); //tolerance of the error in the yield criterion
    IR_GIVE_OPTIONAL_FIELD(ir, newtonIter, IFT_DruckerPragerCutMat_newtonIter, "newtoniter"); //Maximum number of iterations in lambda search
    
    return IRRT_OK;
}


// creates a new material status corresponding to this class
MaterialStatus *
DruckerPragerCutMat :: CreateStatus(GaussPoint *gp) const
{
    MPlasticMaterial2Status *status;
    status = new MPlasticMaterial2Status(1, this->giveDomain(), gp);
    
    return status;
}

double
DruckerPragerCutMat :: computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector, const FloatArray &strainSpaceHardeningVariables)
{
    //strainSpaceHardeningVariables = kappa
    if (isurf <=3){//Rankine, surfaces 1,2,3
        FloatArray princStress(3);
        this->computePrincipalValues(princStress, stressVector, principal_stress);
        return princStress.at(isurf) - this->sigT;
    } else {//Drucker-Prager, surface 4
        double volumetricStress;
        double DPYieldStressInShear = tau0 + H * strainSpaceHardeningVariables.at(4);
        double JTwo;
        MaterialMode mmode = gp->giveMaterialMode();
        StressVector stressVector1(stressVector, mmode);//convert from array
        StressVector deviatoricStress(mmode);
        
        stressVector1.computeDeviatoricVolumetricSplit(deviatoricStress,volumetricStress);
        JTwo = deviatoricStress.computeSecondInvariant();
        return 3.*alpha*volumetricStress + sqrt(JTwo) - DPYieldStressInShear;
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
    if (isurf <=3){ //Rankine associated
        answer.at(1) = t.at(1, isurf) * t.at(1, isurf);//xx = 11
        answer.at(2) = t.at(2, isurf) * t.at(2, isurf);//yy = 22
        answer.at(3) = t.at(3, isurf) * t.at(3, isurf);//zz = 33
        answer.at(4) = t.at(2, isurf) * t.at(3, isurf);//yz = 23
        answer.at(5) = t.at(1, isurf) * t.at(3, isurf);//xz = 13
        answer.at(6) = t.at(1, isurf) * t.at(2, isurf);//xy = 12
    } else { //DP nonassociated
        MaterialMode mmode = gp->giveMaterialMode();
        StressVector stressVector1(stressVector, mmode);//convert from array
        StressVector deviatoricStress(mmode);
        double sqrtJTwo, volumetricStress;
        
        stressVector1.computeDeviatoricVolumetricSplit(deviatoricStress,volumetricStress);
        sqrtJTwo = sqrt(deviatoricStress.computeSecondInvariant());
        
        answer.at(1) = alphaPsi + deviatoricStress.at(1) / 2. / sqrtJTwo;
        answer.at(2) = alphaPsi + deviatoricStress.at(2) / 2. / sqrtJTwo;
        answer.at(3) = alphaPsi + deviatoricStress.at(3) / 2. / sqrtJTwo;
        answer.at(4) = stressVector1.at(4) / sqrtJTwo;
        answer.at(5) = stressVector1.at(5) / sqrtJTwo;
        answer.at(6) = stressVector1.at(6) / sqrtJTwo;
    }
}

//necesarry only for mpm_ClosestPoint, see Jirasek: Inelastic analysis of structures, pp. 411.
//Hessian matrix
void
DruckerPragerCutMat :: computeReducedSSGradientMatrix(FloatMatrix &gradientMatrix,  int isurf, GaussPoint *gp, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables)
{
    gradientMatrix.resize(6, 6);
    gradientMatrix.zero();
    
    if ( isurf == 4 ) {
        int i,j;
        double c1=0.;
        MaterialMode mmode = gp->giveMaterialMode();
        StressVector stressVector1(fullStressVector, mmode);//convert from array
        StressVector deviatoricStress(mmode);
        double JTwo, sqrtJTwo, volumetricStress;
        
        stressVector1.computeDeviatoricVolumetricSplit(deviatoricStress,volumetricStress);
        JTwo = deviatoricStress.computeSecondInvariant();
        sqrtJTwo = sqrt(JTwo);
        
        for(i=1; i<=6;i++){
            for(j=i; j<=6;j++){
                if( (i==1 && j==1) || (i==2 && j==2) || (i==3 && j==3) ){
                    c1 = 2/3.;
                } else if ( (i==4 && j==4) || (i==5 && j==5) || (i==6 && j==6) ){
                    c1 = 1.;
                } else if (i<=3 && j<=3){
                    c1 = -1/3.;
                } else{
                    c1=0;
                }
                
                gradientMatrix.at(i, j) = 0.5/JTwo * ( c1*sqrtJTwo - deviatoricStress.at(i)*deviatoricStress.at(j)/2./sqrtJTwo );
            }
        }
        gradientMatrix.symmetrized();
    }
}


void
DruckerPragerCutMat :: computeReducedElasticModuli(FloatMatrix &answer,
                                         GaussPoint *gp,
                                         TimeStep *atTime)
{  /* Returns elastic moduli in reduced stress-strain space*/
    //MaterialMode mode = gp->giveMaterialMode();
    
    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, ReducedForm,
                                                                    ElasticStiffness,
                                                                    gp, atTime);
}

//answer is dkappa (cumulative plastic strain), flow rule
void DruckerPragerCutMat :: computeStrainHardeningVarsIncrement(FloatArray &answer, GaussPoint *gp, const FloatArray &stress, const FloatArray &dlambda, const FloatArray &dplasticStrain, const IntArray &activeConditionMap){
    answer.resize(4);
    answer.zero();
    
    if( dlambda.at(4) > 0. ){
        answer.at(4) = dlambda.at(4) * sqrt( 1./3. + 2*alphaPsi*alphaPsi );
    }
}

// Computes the derivative of yield/loading function with respect to kappa_1, kappa_2 etc.
void DruckerPragerCutMat :: computeKGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables){
    answer.resize( 1 );//1 hardening variable for DP model - kappa
    answer.zero();
    
    if (isurf == 4){
        answer.at(1) = this->H;
    }
}

//necesarry only for mpm_ClosestPoint
//Computes second mixed derivative of loading function with respect to stress and hardening vars. 
void DruckerPragerCutMat :: computeReducedSKGradientMatrix(FloatMatrix &gradientMatrix, int isurf, GaussPoint *gp, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables){
    gradientMatrix.resize(6, 1);//six stresses and one kappa
    gradientMatrix.zero();
}

// computes dKappa_i/dsig_j gradient matrix
void DruckerPragerCutMat :: computeReducedHardeningVarsSigmaGradient(FloatMatrix &answer, GaussPoint *gp, const IntArray &activeConditionMap, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVars, const FloatArray &dlambda){
    answer.resize(1, 6);
    answer.zero();
}

// computes dKappa_i/dLambda_j for one surface
void DruckerPragerCutMat :: computeReducedHardeningVarsLamGradient(FloatMatrix &answer, GaussPoint *gp, int actSurf, const IntArray &activeConditionMap, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVars, const FloatArray &dlambda){
    int indx;
    answer.resize(1, actSurf);//actSurf = number of active surfaces
    answer.zero();
    
     if ( ( indx = activeConditionMap.at(4) ) ) {
        if ( dlambda.at(4) > 0. ) {
            answer.at(1, indx) = sqrt( 1./3. + 2*alphaPsi*alphaPsi );
        }
    }
}


int
DruckerPragerCutMat :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    //MaterialStatus *status = this->giveStatus(aGaussPoint);
    if ( ( type == IST_DamageScalar ) || (type == IST_DamageTensor ) ) {
        answer.resize(1);
        answer.zero();
        //answer.at(1) = status->giveDamage();
        return 1;
    } else {
        return MPlasticMaterial2 :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}



InternalStateValueType
DruckerPragerCutMat :: giveIPValueType(InternalStateType type)
{
    if ( type == IST_PlasticStrainTensor ) {
        return ISVT_TENSOR_S3;
    } else if ( ( type == IST_MaxEquivalentStrainLevel ) || ( type == IST_DamageScalar ) ) {
        return ISVT_SCALAR;
    } else {
        return StructuralMaterial :: giveIPValueType(type);
    }
}


int
DruckerPragerCutMat :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( type == IST_PlasticStrainTensor ) {
        if ( ( mmode == _3dMat ) || ( mmode == _3dMat_F ) ) {
            answer.resize(6);
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 3;
            answer.at(4) = 4;
            answer.at(5) = 5;
            answer.at(6) = 6;
            return 1;
        } else if ( mmode == _PlaneStrain ) {
            answer.resize(6);
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 3;
            answer.at(6) = 4;
            return 1;
        } else if ( mmode == _1dMat ) {
            answer.resize(1);
            answer.at(1) = 1;
            return 1;
        }
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else if ( ( type == IST_DamageScalar ) || ( type == IST_DamageTensor) ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    }

    return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
}


int
DruckerPragerCutMat :: giveIPValueSize(InternalStateType type, GaussPoint *gp)
{
    if ( type == IST_PlasticStrainTensor ) {
        MaterialMode mode = gp->giveMaterialMode();
        if (mode == _3dMat || mode == _3dMat_F)
            return 6;
        else if (mode == _PlaneStrain)
            return 4;
        else if (mode == _PlaneStress)
            return 3;
        else if (mode == _1dMat)
            return 1;
        else
            return 0;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        return 1;
    } else if ( type == IST_DamageScalar || type == IST_DamageTensor) {
        return 1;
    } else {
        return StructuralMaterial :: giveIPValueSize(type, gp);
    }
}

} // end namespace oofem
