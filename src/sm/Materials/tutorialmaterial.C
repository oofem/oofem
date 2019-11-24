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
#include "sm/Elements/structuralelement.h"
#include "mathfem.h"

namespace oofem {
REGISTER_Material(TutorialMaterial);

TutorialMaterial :: TutorialMaterial(int n, Domain *d) : StructuralMaterial(n, d), D(n, d)
{}


void
TutorialMaterial :: initializeFrom(InputRecord &ir)
{
    StructuralMaterial :: initializeFrom(ir);

    D.initializeFrom(ir);
    IR_GIVE_FIELD(ir, this->sig0, _IFT_TutorialMaterial_yieldstress);
    IR_GIVE_FIELD(ir, this->H, _IFT_TutorialMaterial_hardeningmoduli);
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
    return new TutorialMaterialStatus(gp);
}


FloatArrayF<6>
TutorialMaterial :: giveRealStressVector_3d(const FloatArrayF<6> &totalStrain, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< TutorialMaterialStatus * >( this->giveStatus(gp) );

    // subtract stress thermal expansion
    auto thermalStrain = this->computeStressIndependentStrainVector_3d(gp, tStep, VM_Total);
    auto strain = totalStrain - thermalStrain;

    auto trialElastStrain = strain - status->givePlasticStrain();

    // Compute trial stress = sig_tr = sig_old + E*delta_eps
    const auto &elasticStiffness = D.giveTangent();
    auto trialStress = dot(elasticStiffness, trialElastStrain);

    //auto [devTrialStress, meanTrialStress] = computeDeviatoricVolumetricSplit(); // c++17
    auto tmp = computeDeviatoricVolumetricSplit(trialStress);
    auto devTrialStress = tmp.first;
    auto meanTrialStress = tmp.second;

    double J2 = this->computeSecondStressInvariant(devTrialStress);
    double effectiveTrialStress = sqrt(3 * J2);

    // evaluate the yield surface
    double k = status->giveK();
    double phiTrial = effectiveTrialStress - ( this->sig0 +  H * k );

    FloatArrayF<6> stress;
    if ( phiTrial < 0.0 ) { // elastic
        stress = trialStress;

        status->letTempPlasticStrainBe(status->givePlasticStrain());
    } else { // plastic loading
        double G = D.giveShearModulus();
        double mu = phiTrial / ( 3.0 * G + H ); // plastic multiplier
        // radial return
        auto devStress = ( 1.0 - 3.0*G*mu/effectiveTrialStress) * devTrialStress;
        stress = computeDeviatoricVolumetricSum(devStress, meanTrialStress);
        k += mu;

        auto plasticStrain = status->givePlasticStrain();
        auto dPlStrain = applyDeviatoricElasticCompliance(devTrialStress, 0.5);
        plasticStrain += (mu * 3. / (2. * effectiveTrialStress)) * dPlStrain;
        status->letTempPlasticStrainBe(plasticStrain);
    }

    // Store the temporary values for the given iteration
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(stress);
    status->letTempKBe(k);
    status->letTempDevTrialStressBe(devTrialStress);
    return stress;
}


FloatMatrixF<6,6>
TutorialMaterial :: give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< TutorialMaterialStatus * >( this->giveStatus(gp) );
    const auto &devTrialStress = status->giveTempDevTrialStress();

    double J2 = this->computeSecondStressInvariant(devTrialStress);
    double effectiveTrialStress = sqrt(3 * J2);

    // evaluate the yield surface
    double k = status->giveK();
    double phiTrial = effectiveTrialStress - ( this->sig0 + H * k );
    
    const auto &elasticStiffness = D.giveTangent();
    
    if ( phiTrial < 0.0 ) { // elastic
        return elasticStiffness;
    } else { // plastic loading
        double G = D.giveShearModulus();
        // E_t = elasticStiffness - correction
        // correction =  2.0 * G * ( 2.0 * G / h *( sig0 + kappa ) / sigtre *openProd(nu,nu) + mu*3*G/sigtre *Idev);
        double h = 3.0 * G + H;
        double mu = phiTrial / h; // plasic multiplier
        auto nu = ( 3.0/2.0 / effectiveTrialStress ) * devTrialStress; 

        auto correction = dyad(nu, nu) * (2.0 * G / h * ( this->sig0 + H * k ));
        correction += (mu * 3.0 * G) * I_dev6;
        correction *= 2.0 * G / effectiveTrialStress;
        return elasticStiffness - correction;
    }
}


FloatArrayF<6>
TutorialMaterial :: giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const
{
    return D.giveAlpha();
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


TutorialMaterialStatus :: TutorialMaterialStatus(GaussPoint * g) :
    StructuralMaterialStatus(g)
{
    strainVector.resize(6);
    stressVector.resize(6);
    tempStressVector = stressVector;
    tempStrainVector = strainVector;
}

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
    tempK = k;
    tempDevTrialStress = zeros<6>();
}


void 
TutorialMaterialStatus :: updateYourself(TimeStep *tStep)
{
    // Copy the temporary values to the convered ones. 
    // This method is called after equilibrium has been reached and should 
    // set all...
    StructuralMaterialStatus :: updateYourself(tStep);

    plasticStrain = tempPlasticStrain;
    k = tempK;
    // deviatoric trial stress is not really a state variable and was used not to repeat some code...
}

} // end namespace oofem
