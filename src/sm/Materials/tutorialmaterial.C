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

#include "tutorialmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(TutorialMaterial);

TutorialMaterial :: TutorialMaterial(int n, Domain *d) :StructuralMaterial(n, d),  D(n, d)
{}

TutorialMaterial :: ~TutorialMaterial()
{}

IRResultType
TutorialMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    StructuralMaterial :: initializeFrom(ir);
    D.initializeFrom(ir);
    
    IR_GIVE_FIELD(ir, this->sig0, _IFT_TutorialMaterial_yieldstress);

    IR_GIVE_FIELD(ir, this->H, _IFT_TutorialMaterial_hardeningmoduli);

    return IRRT_OK;
}


void
TutorialMaterial :: giveInputRecord(DynamicInputRecord &ir)
{
    StructuralMaterial :: giveInputRecord(ir);
    D.giveInputRecord(ir);
    
    ir.setField(this->sig0, _IFT_TutorialMaterial_yieldstress);
    ir.setField(this->H, _IFT_TutorialMaterial_hardeningmoduli);
}


MaterialStatus *
TutorialMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new TutorialMaterialStatus(1, this->giveDomain(), gp);
}


void
TutorialMaterial :: giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                 const FloatArray &totalStrain, TimeStep *tStep)
{
    TutorialMaterialStatus *status = static_cast< TutorialMaterialStatus * >( this->giveStatus(gp) );
     
    FloatArray deltaStrain;
    deltaStrain.beDifferenceOf(totalStrain, status->giveStrainVector());  
    
    // Compute trial stress = sig_tr = sig_old + E*delta_eps
    FloatMatrix elasticStiffness;
    D.give3dMaterialStiffnessMatrix(elasticStiffness, TangentStiffness, gp, tStep);
    FloatArray trialStress, deltaStress;
    deltaStress.beProductOf(elasticStiffness, deltaStrain);   // linear increment
    trialStress = status->giveStressVector(); // last equilibriated stress
    trialStress.add(deltaStress);

    FloatArray sphTrialStress, devTrialStress;
    computeSphDevPartOf(trialStress, sphTrialStress, devTrialStress);
    
    double J2 = this->computeJ2InvariantOf(devTrialStress);
    double effectiveTrialStress = sqrt(3 * J2);
    
    // evaluate the yield surface
    double kappa = status->giveKappa();
    double phiTrial = effectiveTrialStress - ( this->sig0 + kappa );
    
    if ( phiTrial < 0.0 ) { // elastic
        answer = trialStress;
    } else { // plastic loading
        double G = D.giveShearModulus();
        double mu = phiTrial / ( 3.0 * G + H ); // plastic multiplier
        answer = devTrialStress * ( 1.0 - 3.0*G*mu/effectiveTrialStress); // radial return
        answer.add(sphTrialStress);
        kappa += H*mu; 
    }
    
    // Store the temporary values for the given iteration
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    status->letTempKappaBe(kappa);
    status->letTempDevTrialStressBe(devTrialStress);
}


void
TutorialMaterial :: computeSphDevPartOf(const FloatArray &sigV, FloatArray &sigSph, FloatArray &sigDev)
{
    
    double volumetricPart = ( sigV.at(1) + sigV.at(2) + sigV.at(3) ) / 3.0;
    sigSph = {volumetricPart, volumetricPart, volumetricPart, 0.0, 0.0, 0.0};
   
    sigDev.beDifferenceOf(sigV, sigSph);
}


double
TutorialMaterial :: computeJ2InvariantOf(const FloatArray &sigV)
{
    // Computes the second invariant J2 of the stress tensor (in Voigt format) 
    // J2 = 0.5 * sigma : sigma = 0.5 * |sig|^2

    double normSig = sigV.at(1)*sigV.at(1) + sigV.at(2)*sigV.at(2) + sigV.at(3)*sigV.at(3) +
             2.0 * ( sigV.at(4)*sigV.at(4) + sigV.at(5)*sigV.at(5) + sigV.at(6)*sigV.at(6) );
    
    return 0.5 * normSig;
}


void 
TutorialMaterial :: giveDeviatoricProjectionMatrix(FloatMatrix &answer)
{
    answer.resize(6,6);
    answer.zero();
    answer.at(1,1) = answer.at(2,2) = answer.at(3,3) =  2.0/3.0;
    answer.at(1,2) = answer.at(1,3) = answer.at(2,3) = -1.0/3.0;
    answer.at(2,1) = answer.at(3,1) = answer.at(3,2) = -1.0/3.0;
    answer.at(4,4) = answer.at(5,5) = answer.at(6,6) =  1.0/2.0;
}


void
TutorialMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
 
    TutorialMaterialStatus *status = static_cast< TutorialMaterialStatus * >( this->giveStatus(gp) );
    FloatArray devTrialStress = status->giveTempDevTrialStress();

    double J2 = this->computeJ2InvariantOf(devTrialStress);
    double effectiveTrialStress = sqrt(3) * sqrt(J2);
    
    // evaluate the yield surface
    double kappa = status->giveKappa();
    double phiTrial = effectiveTrialStress - ( this->sig0 + kappa );
    
    FloatMatrix elasticStiffness;
    D.give3dMaterialStiffnessMatrix(elasticStiffness, mode, gp, tStep);
    
    if ( phiTrial < 0.0 ) { // elastic
        answer = elasticStiffness;
    } else { // plastic loading
        double G = D.giveShearModulus();
        // E_t = elasticStiffness - correction
        // correction =  2.0 * G * ( 2.0 * G / h *( sig0 + kappa ) / sigtre *openProd(nu,nu) + mu*3*G/sigtre *Idev);
        double h = 3.0 * G + H;
        double mu = phiTrial / h; // plasic multiplier
        
        FloatArray nu = (3.0/2.0 / effectiveTrialStress ) * devTrialStress; 
        FloatMatrix Idev, correction;
        giveDeviatoricProjectionMatrix(Idev);
        
        correction.plusDyadUnsym(nu, nu, 2.0 * G / h * ( this->sig0 + kappa ));
        correction.add(mu * 3.0 * G, Idev);
        correction.times(2.0 * G / effectiveTrialStress);
        answer = elasticStiffness;
        answer.subtract(correction);
    }
}


int
TutorialMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    TutorialMaterialStatus *status = static_cast< TutorialMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_PlasticStrainTensor ) {
        answer = status->givePlasticStrain();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}


//=============================================================================


TutorialMaterialStatus :: TutorialMaterialStatus(int n, Domain * d, GaussPoint * g) :
    StructuralMaterialStatus(n, d, g),
    tempPlasticStrain(6), plasticStrain(6),
    tempDevTrialStress(6),
    tempKappa(0.), kappa(0.)
{ }

void
TutorialMaterialStatus :: initTempStatus()
{
    //StructuralMaterialStatus :: initTempStatus();

    // reset temp vars.
    tempStressVector = stressVector;
    tempStrainVector = strainVector;
    tempPVector      = PVector;
    tempFVector      = FVector;


    tempPlasticStrain = plasticStrain;
    tempKappa = kappa;
    tempDevTrialStress.resize(6);
    tempDevTrialStress.zero();
}


void 
TutorialMaterialStatus :: updateYourself(TimeStep *tStep)
{
    // Copy the temporary values to the convered ones. 
    // This method is called after equilibrium has been reached and should 
    // set all...
    StructuralMaterialStatus :: updateYourself(tStep);

    plasticStrain = tempPlasticStrain;
    kappa = tempKappa;
    // deviatoric trial stress is not really a state variable and was used not to repeat some code...
}

} // end namespace oofem
