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
 *               Copyright (C) 1993 - 2021   Borek Patzak
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

#include "structuralslipfe2material.h"
#include "gausspoint.h"
#include "engngm.h"
#include "oofemtxtdatareader.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "util.h"
#include "contextioerr.h"
#include "generalboundarycondition.h"
#include "prescribeddispsliphomogenization.h"
#include "exportmodulemanager.h"
#include "vtkxmlexportmodule.h"
#include "nummet.h"
#include "sm/EngineeringModels/xfemsolverinterface.h"
#include "sm/EngineeringModels/staticstructural.h"
#include "unknownnumberingscheme.h"
#include "xfem/xfemstructuremanager.h"
#include "mathfem.h"

#include "dynamicdatareader.h"

#include <sstream>

namespace oofem {
REGISTER_Material(StructuralSlipFE2Material);

StructuralSlipFE2Material :: StructuralSlipFE2Material(int n, Domain *d) : StructuralFE2Material(n, d)
{}


void
StructuralSlipFE2Material :: initializeFrom(InputRecord &ir)
{
    StructuralFE2Material :: initializeFrom(ir);

    useExtStiff = ir.hasField(_IFT_StructuralSlipFE2Material_useExternalStiffness);
    allGPRes = ir.hasField(_IFT_StructuralSlipFE2Material_allGPResults);
    IR_GIVE_OPTIONAL_FIELD(ir, outputSelected, _IFT_StructuralSlipFE2Material_outputSelectedResults);
    IR_GIVE_OPTIONAL_FIELD(ir, givendStressdEpsTangent, _IFT_StructuralSlipFE2Material_dStressdEps);
    IR_GIVE_OPTIONAL_FIELD(ir, givendStressdEpsTangent, _IFT_StructuralSlipFE2Material_dStressdEps);
    IR_GIVE_OPTIONAL_FIELD(ir, givendBStressdEpsTangent, _IFT_StructuralSlipFE2Material_dBStressdEps);
    IR_GIVE_OPTIONAL_FIELD(ir, givendRStressdEpsTangent, _IFT_StructuralSlipFE2Material_dRStressdEps);
    IR_GIVE_OPTIONAL_FIELD(ir, givendStressdSTangent, _IFT_StructuralSlipFE2Material_dStressdS);
    IR_GIVE_OPTIONAL_FIELD(ir, givendBStressdSTangent, _IFT_StructuralSlipFE2Material_dBStressdS);
    IR_GIVE_OPTIONAL_FIELD(ir, givendRStressdSTangent, _IFT_StructuralSlipFE2Material_dRStressdS);
    IR_GIVE_OPTIONAL_FIELD(ir, givendStressdGTangent, _IFT_StructuralSlipFE2Material_dStressdG);
    IR_GIVE_OPTIONAL_FIELD(ir, givendBStressdGTangent, _IFT_StructuralSlipFE2Material_dBStressdG);
    IR_GIVE_OPTIONAL_FIELD(ir, givendRStressdGTangent, _IFT_StructuralSlipFE2Material_dRStressdG);
}


void
StructuralSlipFE2Material :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralFE2Material :: giveInputRecord(input);

    if ( useExtStiff ) {
        input.setField(_IFT_StructuralSlipFE2Material_useExternalStiffness);
        input.setField(givendStressdEpsTangent, _IFT_StructuralSlipFE2Material_dStressdEps);
    }

    input.setField(allGPRes, _IFT_StructuralSlipFE2Material_allGPResults);
    input.setField(outputSelected, _IFT_StructuralSlipFE2Material_outputSelectedResults);
}


MaterialStatus *
StructuralSlipFE2Material :: CreateStatus(GaussPoint *gp) const
{
    int rank = -1;
    auto emodel = this->domain->giveEngngModel();
    if ( emodel->isParallel() && emodel->giveNumberOfProcesses() > 1 ) {
        rank = emodel->giveRank();
    }

    if ( !(outputSelected.giveSize()==0) ) {
        std::vector<GaussPoint*> gpArray;
        gpArray.resize(outputSelected.giveSize()/2);
        for (int i = 1; i <= outputSelected.giveSize()/2; i++) {
            GaussPoint* gpOut = domain->giveGlobalElement(outputSelected.at(i*2-1))->giveIntegrationRule(0)->getIntegrationPoint(outputSelected.at(i*2) - 1);
            gpArray[i - 1] = gpOut;
        }

        if ( std::find(gpArray.begin(), gpArray.end(), gp) != gpArray.end() ) {
            int nel = gp->giveElement()->giveGlobalNumber();
            int gpn = gp->giveNumber();
            return new StructuralSlipFE2MaterialStatus(rank, gp, this->inputfile, nel, gpn);
        } else {
            return new StructuralSlipFE2MaterialStatus(rank, gp, this->inputfile, 1, 1);
        }
    } else if ( allGPRes ) {
        int nel = gp->giveElement()->giveGlobalNumber();
        int gpn = gp->giveNumber();
        return new StructuralSlipFE2MaterialStatus(rank, gp, this->inputfile, nel, gpn);
    } else {
        return new StructuralSlipFE2MaterialStatus(rank, gp, this->inputfile, 1, 1);
    }

}

FloatMatrixF<3, 3>
StructuralSlipFE2Material::givePlaneStressStiffMtrx( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const
{
    if (useExtStiff) {
        auto status = static_cast<StructuralSlipFE2MaterialStatus*>(this->giveStatus(gp));
        FloatMatrixF<3,3> answer;
        answer = this->givendStressdEpsTangent;
        status->setdStressdEpsTangent(answer);
        return answer;
    } else if (useNumTangent) {
        return StructuralFE2Material::givePlaneStressStiffMtrx(rMode, gp, tStep);
    } else {
        auto status = static_cast<StructuralSlipFE2MaterialStatus*>(this->giveStatus(gp));
        FloatMatrix tangent;
        status->computeTangent(tStep); //implemented by BC
        tangent.beSubMatrixOf(status->givedStressdEpsTangent(), {1,2,3}, {1,2,3});
        FloatMatrixF<3,3> answer;
        answer = tangent;
        return answer;

    }
}


FloatArrayF<3> StructuralSlipFE2Material::giveRealStressVector_PlaneStress( const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep ) const
{
    auto ms = static_cast<StructuralSlipFE2MaterialStatus*>( this->giveStatus(gp) );

    ms->setTimeStep(tStep);
    ms->giveBC()->setDispGradient(strain);
    // Solve subscale problem
    ms->giveRVE()->solveYourselfAt(tStep);
    // Post-process the stress
    FloatArray stress;
    ms->giveBC()->computeStress(stress, tStep);

    FloatArrayF<6> updateStress;
    updateStress = {stress[0], stress[1], 0., 0., 0., 0.5*(stress[2]+stress[3])};

    FloatArrayF<6> updateStrain;
    updateStrain = {strain[0], strain[1], 0., 0., 0., strain[2]};

    FloatArrayF<3> answer;
    answer = {stress[0], stress[1], 0.5*(stress[2]+stress[3])};

    // Update the material status variables
    ms->letTempStressVectorBe(updateStress);
    ms->letTempStrainVectorBe(updateStrain);
    ms->markOldTangent(); // Mark this so that tangent is reevaluated if they are needed.

    return answer;
}


void StructuralSlipFE2Material::giveHomogenizedFields( FloatArray &stress, FloatArray &bStress, FloatArray &rStress, const FloatArray &strain, const FloatArray &slip, const FloatArray &slipGradient, GaussPoint *gp, TimeStep *tStep )
{
    StructuralSlipFE2MaterialStatus *ms = static_cast<StructuralSlipFE2MaterialStatus *>( this->giveStatus(gp) );

    ms->setTimeStep(tStep);
    //Set input
    ms->giveBC()->setDispGradient(strain);
    ms->giveBC()->setSlipField(slip);
    ms->giveBC()->setSlipGradient(slipGradient);
//     Solve subscale problem
//     OOFEM_LOG_INFO("Solving subscale problem at element %d, gp %d.\n Applied strain is %10.6e, %10.6e, %10.6e \n Applied slip is %10.6e, %10.6e \n Applied slip gradient is %10.6e, %10.6e, %10.6e, %10.6e \n", gp->giveElement()->giveGlobalNumber(), gp->giveNumber(), strain.at(1), strain.at(2), strain.at(3), slip.at(1), slip.at(2), slipGradient.at(1), slipGradient.at(2), slipGradient.at(3), slipGradient.at(4) );
    ms->giveRVE()->solveYourselfAt(tStep);
//     OOFEM_LOG_INFO("Solution of RVE problem @ element %d, gp %d OK. \n", gp->giveElement()->giveGlobalNumber(), gp->giveNumber() );
    //Homogenize fields
    FloatArray stress4;
    ms->giveBC()->computeStress(stress4, tStep);
    ms->giveBC()->computeTransferStress(bStress, tStep);
    ms->giveBC()->computeReinfStress(rStress, tStep);
//    OOFEM_LOG_INFO("Stress is  %10.6e, %10.6e, %10.6e \n Transfer stress is %10.6e, %10.6e \n Reinforcement membrane stress is %10.6e, %10.6e, %10.6e, %10.6e \n",  stress4.at(1), stress4.at(2), stress4.at(3), bStress.at(1), bStress.at(2), rStress.at(1), rStress.at(2), rStress.at(3), rStress.at(4) );

    if (stress4.giveSize() == 4 ) {
        stress = {stress4[0], stress4[1], 0.5*(stress4[2]+stress4[3])};
    } else {
        OOFEM_ERROR("Only 2D plane stress mode supported");
    }

    // Update the material status variables
    ms->letTempStressVectorBe(stress);
    ms->letTempStrainVectorBe(strain);
    ms->letTempTransferStressVectorBe(bStress);
    ms->letTempSlipVectorBe(slip);
    ms->letTempReinfStressVectorBe(rStress);
    ms->letTempSlipGradVectorBe(slipGradient);
    ms->markOldTangent(); // Mark this so that tangent is reevaluated if they are needed.
}


void StructuralSlipFE2Material::giveSensitivities( FloatMatrix &dStressdEps, FloatMatrix &dStressdS, FloatMatrix &dStressdG, FloatMatrix &dBStressdEps, FloatMatrix &dBStressdS,
    FloatMatrix &dBStressdG, FloatMatrix &dRStressdEps, FloatMatrix &dRStressdS, FloatMatrix &dRStressdG, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
    if ( useNumTangent ) {
        //Numerical tangent
        StructuralSlipFE2MaterialStatus *status = static_cast<StructuralSlipFE2MaterialStatus*>( this->giveStatus( gp ) );
        double h = 1.0e-9;

        //get current values from material status
        const FloatArray &epsRed = status->giveTempStrainVector();
        const FloatArray &slipRed = status->giveTempSlipVector();
        const FloatArray &slipGradRed = status->giveTempSlipGradVector();

        FloatArray eps;
        if ( epsRed.giveSize() != 3 ) {
            eps = { epsRed[0], epsRed[1], epsRed[5] };
        } else {
            eps = epsRed;
        }
        FloatArray slip=slipRed;
        FloatArray slipGrad=slipGradRed;

        int dim1 = eps.giveSize();
        int dim2 = slip.giveSize();
        int dim3 = slipGrad.giveSize();
        dStressdEps.resize(dim1,dim1);  dStressdEps.zero();
        dBStressdEps.resize(dim2,dim1); dBStressdEps.zero();
        dRStressdEps.resize(dim3,dim1); dRStressdEps.zero();
        dStressdS.resize(dim1,dim2);    dStressdS.zero();
        dBStressdS.resize(dim2,dim2);   dBStressdS.zero();
        dRStressdS.resize(dim3,dim2);   dRStressdS.zero();
        dStressdG.resize(dim1,dim3);    dStressdG.zero();
        dBStressdG.resize(dim2,dim3);   dBStressdG.zero();
        dRStressdG.resize(dim3,dim3);   dRStressdG.zero();

        FloatArray sig, sigPert, epsPert;
        FloatArray bsig, bsigPert, slipPert;
        FloatArray rsig, rsigPert, slipGradPert;

        for(int i = 1; i <= dim1; i++) {
            // Add a small perturbation to the strain
            epsPert = eps;
            epsPert.at(i) += h;

            giveHomogenizedFields( sigPert, bsigPert, rsigPert, epsPert, slip, slipGrad, gp, tStep );
            dStressdEps.setColumn( sigPert, i);
            dBStressdEps.setColumn( bsigPert, i);
            dRStressdEps.setColumn( rsigPert, i);
        }

        for(int i = 1; i <= dim2; i++) {
            // Add a small perturbation to the slip
            slipPert = slip;
            slipPert.at(i) += h;

            giveHomogenizedFields( sigPert, bsigPert, rsigPert, eps, slipPert, slipGrad, gp, tStep );
            dStressdS.setColumn( sigPert, i );
            dBStressdS.setColumn( bsigPert, i);
            dRStressdS.setColumn( rsigPert, i);
        }

        for(int i = 1; i <= dim3; i++) {
            // Add a small perturbation to slip gradient
            slipGradPert = slipGrad;
            slipGradPert.at(i) += h;

            giveHomogenizedFields( sigPert, bsigPert, rsigPert, eps, slip, slipGradPert, gp, tStep );
            dStressdG.setColumn( sigPert, i);
            dBStressdG.setColumn( bsigPert, i);
            dRStressdG.setColumn( rsigPert, i);
        }

        giveHomogenizedFields( sig, bsig, rsig, eps, slip, slipGrad, gp, tStep);

        for(int i = 1; i <= dim1; i++) {
            for(int j = 1; j <= dim1; j++) {
                dStressdEps.at(j,i) -= sig.at(j);
                dStressdEps.at(j,i) /= h;
            }
            for(int j = 1; j <= dim2; j++) {
                dBStressdEps.at(j,i) -= bsig.at(j);
                dBStressdEps.at(j,i) /= h;
            }
            for(int j = 1; j <= dim3; j++) {
                dRStressdEps.at(j,i) -= rsig.at(j);
                dRStressdEps.at(j,i) /= h;
            }
        }

        for(int i = 1; i <= dim2; i++) {
            for(int j = 1; j <= dim1; j++) {
                dStressdS.at(j,i) -= sig.at(j);
                dStressdS.at(j,i) /= h;
            }
            for(int j = 1; j <= dim2; j++) {
                dBStressdS.at(j,i) -= bsig.at(j);
                dBStressdS.at(j,i) /= h;
            }
            for(int j = 1; j <= dim3; j++) {
                dRStressdS.at(j,i) -= rsig.at(j);
                dRStressdS.at(j,i) /= h;
            }
        }

        for(int i = 1; i <= dim3; i++) {
            for(int j = 1; j <= dim1; j++) {
                dStressdG.at(j,i) -= sig.at(j);
                dStressdG.at(j,i) /= h;
            }
            for(int j = 1; j <= dim2; j++) {
                dBStressdG.at(j,i) -= bsig.at(j);
                dBStressdG.at(j,i) /= h;
            }
            for(int j = 1; j <= dim3; j++) {
                dRStressdG.at(j,i) -= rsig.at(j);
                dRStressdG.at(j,i) /= h;
            }
        }

        status->setdStressdEpsTangent(dStressdEps);
        status->setdBStressdEpsTangent(dBStressdEps);
        status->setdRStressdEpsTangent(dRStressdEps);
        status->setdStressdSTangent(dStressdS);
        status->setdBStressdSTangent(dBStressdS);
        status->setdRStressdSTangent(dRStressdS);
        status->setdStressdGTangent(dStressdG);
        status->setdBStressdGTangent(dBStressdG);
        status->setdRStressdGTangent(dRStressdG);

    } else if ( useExtStiff) {
        StructuralSlipFE2MaterialStatus *ms = static_cast<StructuralSlipFE2MaterialStatus *>( this->giveStatus(gp) );

        dStressdEps = this->givendStressdEpsTangent;
        dBStressdEps = this->givendBStressdEpsTangent;
        dRStressdEps = this->givendRStressdEpsTangent;
        dStressdS = this->givendStressdSTangent;
        dBStressdS = this->givendBStressdSTangent;
        dRStressdS = this->givendRStressdSTangent;
        dStressdG = this->givendStressdGTangent;
        dBStressdG = this->givendBStressdGTangent;
        dRStressdG = this->givendRStressdGTangent;

        ms->setdStressdEpsTangent(dStressdEps);
        ms->setdBStressdEpsTangent(dBStressdEps);
        ms->setdRStressdEpsTangent(dRStressdEps);
        ms->setdStressdSTangent(dStressdS);
        ms->setdBStressdSTangent(dBStressdS);
        ms->setdRStressdSTangent(dRStressdS);
        ms->setdStressdGTangent(dStressdG);
        ms->setdBStressdGTangent(dBStressdG);
        ms->setdRStressdGTangent(dRStressdG);
    } else {
        OOFEM_ERROR("Exact tangent not implemented. Either use_num_tangent or use_ext_stiffness");
    }
}

//=============================================================================


StructuralSlipFE2MaterialStatus :: StructuralSlipFE2MaterialStatus( int rank, GaussPoint * g,  const std :: string & inputfile, int el, int gp ) :
    StructuralMaterialStatus(g),
    mInputFile(inputfile)
{
    if ( !this->createRVE(inputfile, rank, el, gp) ) {
        OOFEM_ERROR("Couldn't create RVE");
    }
    stressVector.resize(6);
    strainVector.resize(6);
    tempStressVector.resize(6);
    tempStrainVector.resize(6);

    //this is now hardcoded for 2d plane stress
    slipVector.resize(2);
    bStressVector.resize(2);
    tempSlipVector.resize(2);
    tempBStressVector.resize(2);

    slipGradVector.resize(4);
    rStressVector.resize(4);
    tempSlipGradVector.resize(4);
    tempRStressVector.resize(4);
}


PrescribedDispSlipHomogenization* StructuralSlipFE2MaterialStatus::giveBC()
{
    this->bc = dynamic_cast< PrescribedDispSlipHomogenization * >( this->rve->giveDomain(1)->giveBc(1) );
    return this->bc;
}


bool StructuralSlipFE2MaterialStatus :: createRVE(const std :: string &inputfile, int rank, int el, int gp)
{
    OOFEMTXTDataReader dr( inputfile );
    this->rve = InstanciateProblem(dr, _processor, 0); // Everything but nrsolver is updated.
    dr.finish();
    this->rve->setProblemScale(microScale);
    this->rve->checkProblemConsistency();
    this->rve->initMetaStepAttributes( this->rve->giveMetaStep(1) );
    this->rve->giveNextStep(); // Makes sure there is a timestep (which we will modify before solving a step)
    this->rve->init();

    std :: ostringstream name;
    name << this->rve->giveOutputBaseFileName() << "_el" << el << "_gp" << gp;
    if ( rank >= 0 ) {
        name << "." << rank;
    }

    this->rve->letOutputBaseFileNameBe( name.str() );

    bc = this->giveBC();
    if ( !bc ) {
        OOFEM_ERROR("RVE doesn't have necessary boundary condition; should have a type of PrescribedDispSlipHomogenization as first b.c.");
    }

    return true;
}

void StructuralSlipFE2MaterialStatus :: setTimeStep(TimeStep *tStep)
{
    TimeStep *rveTStep = this->rve->giveCurrentStep(); // Should i create a new one if it is empty?
    rveTStep->setNumber( tStep->giveNumber() );
    rveTStep->setTime( tStep->giveTargetTime() );
    rveTStep->setTimeIncrement( tStep->giveTimeIncrement() );
}

void StructuralSlipFE2MaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();

    // see if vectors describing reached equilibrium are defined
    if ( this->giveSlipVector().giveSize() == 0 ) {
        slipVector.resize(2);
    }

    if ( this->giveTransferStressVector().giveSize() == 0) {
        bStressVector.resize(2);
    }

    if (this->giveSlipGradVector().giveSize() == 0) {
        slipGradVector.resize(4);
    }

    if (this->giveReinfStressVector().giveSize() == 0) {
        rStressVector.resize(4);
    }

    // reset temp vars.
    tempSlipVector = slipVector;
    tempBStressVector = bStressVector;
    tempSlipGradVector = slipGradVector;
    tempRStressVector = rStressVector;
}


void StructuralSlipFE2MaterialStatus :: computeTangent(TimeStep *tStep)
{
    if ( !tStep->isTheCurrentTimeStep() ) {
        OOFEM_ERROR("Only current timestep supported.");
    }

    if ( this->olddSdETangent ) {
        bc->computeTangent(this->givedStressdEpsTangent(), tStep);
    }

    this->olddSdETangent = false;
}

void StructuralSlipFE2MaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    this->rve->updateYourself(tStep);
    this->rve->terminate(tStep);

    slipVector = tempSlipVector;
    bStressVector = tempBStressVector;
    slipGradVector = tempSlipGradVector;
    rStressVector = tempRStressVector;
}

void StructuralSlipFE2MaterialStatus::markOldTangent()
{
    this->olddSdETangent = true;
    this->olddBSdETangent = true;
    this->olddRSdETangent = true;
    this->olddSdSTangent = true;
    this->olddBSdSTangent = true;
    this->olddRSdSTangent = true;
    this->olddSdGTangent = true;
    this->olddBSdGTangent = true;
    this->olddRSdGTangent = true;
}

} // end namespace oofem
