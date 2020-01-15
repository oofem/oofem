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

#include "structuralfe2material.h"
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
#include "prescribedgradienthomogenization.h"
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
REGISTER_Material(StructuralFE2Material);

int StructuralFE2Material :: n = 1;

StructuralFE2Material :: StructuralFE2Material(int n, Domain *d) : StructuralMaterial(n, d)
{}


void
StructuralFE2Material :: initializeFrom(InputRecord &ir)
{
    StructuralMaterial :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, this->inputfile, _IFT_StructuralFE2Material_fileName);

    useNumTangent = ir.hasField(_IFT_StructuralFE2Material_useNumericalTangent);
}


void
StructuralFE2Material :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);
    input.setField(this->inputfile, _IFT_StructuralFE2Material_fileName);

    if ( useNumTangent ) {
        input.setField(_IFT_StructuralFE2Material_useNumericalTangent);
    }
}


MaterialStatus *
StructuralFE2Material :: CreateStatus(GaussPoint *gp) const
{
    int rank = -1;
    auto emodel = this->domain->giveEngngModel();
    if ( emodel->isParallel() && emodel->giveNumberOfProcesses() > 1 ) {
        rank = emodel->giveRank();
    }
    return new StructuralFE2MaterialStatus(rank, gp, this->inputfile);
}


FloatArrayF<6>
StructuralFE2Material :: giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    auto ms = static_cast< StructuralFE2MaterialStatus * >( this->giveStatus(gp) );

#if 0
    XfemStructureManager *xMan = dynamic_cast<XfemStructureManager*>( ms->giveRVE()->giveDomain(1)->giveXfemManager() );
    if(xMan) {
        printf("Total crack length in RVE: %e\n", xMan->computeTotalCrackLength() );
    }
#endif

    ms->setTimeStep(tStep);
    // Set input
    ms->giveBC()->setPrescribedGradientVoigt(strain);
    // Solve subscale problem
    ms->giveRVE()->solveYourselfAt(tStep);
    // Post-process the stress
    FloatArray stress;
    ms->giveBC()->computeField(stress, tStep);

    FloatArrayF<6> answer;
    if ( stress.giveSize() == 6 ) {
        answer = stress;
    } else if ( stress.giveSize() == 9 ) {
        answer = {stress[0], stress[1], stress[2], 0.5*(stress[3]+stress[6]), 0.5*(stress[4]+stress[7]), 0.5*(stress[5]+stress[8])};
    } else if ( stress.giveSize() == 4 ) {
        // 2D mode is a bit wonky in FE2 code right now. It doesn't actually deal with plane strain out-of-plane component correctly.
        // Instead gives just the 2D state, non-symmetrically.
        // This needs to be cleared up and made consistent, to take out the guesswork on what computeField() actually computes.
        answer = {stress[0], stress[1], 0., 0., 0., 0.5*(stress[2]+stress[3])};
    } else {
        answer = {stress[0], 0., 0., 0., 0., 0.};
    }

    // Update the material status variables
    ms->letTempStressVectorBe(answer);
    ms->letTempStrainVectorBe(strain);
    ms->markOldTangent(); // Mark this so that tangent is reevaluated if they are needed.
    return answer;
}


FloatMatrixF<6,6>
StructuralFE2Material :: give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast<StructuralFE2MaterialStatus*>( this->giveStatus( gp ) );
    if ( useNumTangent ) {
        // Numerical tangent
        double h = 1.0e-9;

        FloatArrayF<6> eps = status->giveTempStrainVector();
        FloatArrayF<6> sig = status->giveTempStressVector();

        FloatMatrixF<6,6> answer;
        for(int i = 0; i < 6; i++) {
            // Add a small perturbation to the strain
            auto epsPert = eps;
            epsPert[i] += h;

            auto sigPert = const_cast<StructuralFE2Material*>(this)->giveRealStressVector_3d(epsPert, gp, tStep);
            answer.setColumn((sigPert - sig) / h, i);
        }
        const_cast<StructuralFE2Material*>(this)->giveRealStressVector_3d(eps, gp, tStep);
        return answer;

    } else {
        status->computeTangent(tStep);
        const auto &ans9 = status->giveTangent();
        FloatMatrix answer;
        StructuralMaterial::giveReducedSymMatrixForm(answer, ans9, _3dMat);
        return answer;

//        printf("ans9: "); ans9.printYourself();
//
//        // Compute the (minor) symmetrized tangent:
//        answer.resize(6, 6);
//        for ( int i = 0; i < 6; ++i ) {
//            for ( int j = 0; j < 6; ++j ) {
//                answer(i, j) = ans9(i, j);
//            }
//        }
//        for ( int i = 0; i < 6; ++i ) {
//            for ( int j = 6; j < 9; ++j ) {
//                answer(i, j-3) += ans9(i, j);
//                answer(j-3, i) += ans9(j, i);
//            }
//        }
//        for ( int i = 6; i < 9; ++i ) {
//            for ( int j = 6; j < 9; ++j ) {
//                answer(j-3, i-3) += ans9(j, i);
//            }
//        }
//        for ( int i = 0; i < 6; ++i ) {
//            for ( int j = 3; j < 6; ++j ) {
//                answer(j, i) *= 0.5;
//                answer(i, j) *= 0.5;
//            }
//        }



#if 0
        // Numerical ATS for debugging
        FloatMatrix numericalATS(6, 6);
        FloatArray dsig;
        // Note! We need a copy of the temp strain, since the pertubations might change it.
        FloatArray tempStrain = ms->giveTempStrainVector();

        FloatArray sig, strain, sigPert;
        giveRealStressVector_3d(sig, gp, tempStrain, tStep);
        double hh = 1e-6;
        for ( int k = 1; k <= 6; ++k ) {
            strain = tempStrain;
            strain.at(k) += hh;
            giveRealStressVector_3d(sigPert, gp, strain, tStep);
            dsig.beDifferenceOf(sigPert, sig);
            numericalATS.setColumn(dsig, k);
        }
        numericalATS.times(1. / hh);
        giveRealStressVector_3d(sig, gp, tempStrain, tStep); // Reset

        //answer.printYourself("Analytical deviatoric tangent");
        //numericalATS.printYourself("Numerical deviatoric tangent");

        numericalATS.subtract(answer);
        double norm = numericalATS.computeFrobeniusNorm();
        if ( norm > answer.computeFrobeniusNorm() * 1e-3 && norm > 0.0 ) {
            OOFEM_ERROR("Error in deviatoric tangent");
        }
#endif
    }
}

FloatMatrixF<4,4>
StructuralFE2Material :: givePlaneStrainStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    /// The plane-strain version of the tangent is overloaded primarily for performance when using the numerical tangent.
    auto status = static_cast<StructuralFE2MaterialStatus*>( this->giveStatus( gp ) );
    if ( useNumTangent ) {
        // Numerical tangent
        double h = 1.0e-9;

        auto eps = FloatArrayF<6>(status->giveTempStrainVector())[{0, 1, 2, 5}];
        auto sig = FloatArrayF<6>(status->giveTempStressVector())[{0, 1, 2, 5}];
        

        FloatMatrixF<4,4> answer;
        for(int i = 0; i < 4; i++) {
            // Add a small perturbation to the strain
            auto epsPert = eps;
            epsPert[i] += h;

            auto sigPert = this->giveRealStressVector_PlaneStrain(epsPert, gp, tStep);
            answer.setColumn((sigPert - sig) / h, i);
        }
        this->giveRealStressVector_PlaneStrain(eps, gp, tStep);
        return answer;

    } else {
        status->computeTangent(tStep);
        return status->giveTangent();
    }
}


//=============================================================================


StructuralFE2MaterialStatus :: StructuralFE2MaterialStatus(int rank, GaussPoint * g,  const std :: string & inputfile) :
    StructuralMaterialStatus(g),
    mInputFile(inputfile)
{
    if ( !this->createRVE(1, inputfile, rank) ) { ///@TODO FIXME createRVE
        OOFEM_ERROR("Couldn't create RVE");
    }
    stressVector.resize(6);
    strainVector.resize(6);
    tempStressVector.resize(6);
    tempStrainVector.resize(6);
}

PrescribedGradientHomogenization* StructuralFE2MaterialStatus::giveBC()
{
    this->bc = dynamic_cast< PrescribedGradientHomogenization * >( this->rve->giveDomain(1)->giveBc(1) );
    return this->bc;
}


bool
StructuralFE2MaterialStatus :: createRVE(int n, const std :: string &inputfile, int rank)
{
    OOFEMTXTDataReader dr( inputfile.c_str() );
    this->rve = InstanciateProblem(dr, _processor, 0); // Everything but nrsolver is updated.
    dr.finish();
    this->rve->setProblemScale(microScale);
    this->rve->checkProblemConsistency();
    this->rve->initMetaStepAttributes( this->rve->giveMetaStep(1) );
    this->rve->giveNextStep(); // Makes sure there is a timestep (which we will modify before solving a step)
    this->rve->init();

    std :: ostringstream name;
    name << this->rve->giveOutputBaseFileName() << "-gp" << n;
    if ( rank >= 0 ) {
        name << "." << rank;
    }

    this->rve->letOutputBaseFileNameBe( name.str() );

    this->bc = dynamic_cast< PrescribedGradientHomogenization * >( this->rve->giveDomain(1)->giveBc(1) );
    if ( !this->bc ) {
        OOFEM_ERROR("RVE doesn't have necessary boundary condition; should have a type of PrescribedGradientHomogenization as first b.c.");
    }

    return true;
}

void
StructuralFE2MaterialStatus :: setTimeStep(TimeStep *tStep)
{
    TimeStep *rveTStep = this->rve->giveCurrentStep(); // Should i create a new one if it is empty?
    rveTStep->setNumber( tStep->giveNumber() );
    rveTStep->setTime( tStep->giveTargetTime() );
    rveTStep->setTimeIncrement( tStep->giveTimeIncrement() );
}

void
StructuralFE2MaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
}

void
StructuralFE2MaterialStatus :: markOldTangent() { this->oldTangent = true; }

void
StructuralFE2MaterialStatus :: computeTangent(TimeStep *tStep)
{
    if ( !tStep->isTheCurrentTimeStep() ) {
        OOFEM_ERROR("Only current timestep supported.");
    }

    if ( this->oldTangent ) {
        bc->computeTangent(this->giveTangent(), tStep);
    }

    this->oldTangent = false;
}

void 
StructuralFE2MaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    this->rve->updateYourself(tStep);
    this->rve->terminate(tStep);

    mNewlyInitialized = false;
}


void
StructuralFE2MaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus :: saveContext(stream, mode);
    this->rve->saveContext(stream, mode);
}


void
StructuralFE2MaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus :: restoreContext(stream, mode);
    this->rve->restoreContext(stream, mode);
}

double StructuralFE2MaterialStatus :: giveRveLength()
{
    return sqrt( bc->domainSize() );
}

void StructuralFE2MaterialStatus :: copyStateVariables(const MaterialStatus &iStatus)
{
    //static int num = 0;
    //printf("Entering StructuralFE2MaterialStatus :: copyStateVariables.\n");

    this->oldTangent = true;

//    if ( !this->createRVE(this->giveNumber(), gp, mInputFile) ) {
//        OOFEM_ERROR("Couldn't create RVE");
//    }


    StructuralMaterialStatus :: copyStateVariables(iStatus);


    //////////////////////////////
    const StructuralFE2MaterialStatus *fe2ms = dynamic_cast<const StructuralFE2MaterialStatus*>(&iStatus);

    if ( !fe2ms ) {
        OOFEM_ERROR("Failed to cast StructuralFE2MaterialStatus.")
    }


    this->mNewlyInitialized = fe2ms->mNewlyInitialized;

    // The proper way to do this would be to clone the RVE from iStatus.
    // However, this is a mess due to all pointers that need to be tracked.
    // Therefore, we consider a simplified version: copy only the enrichment items.

    Domain *ext_domain = fe2ms->giveRVE()->giveDomain(1);
    if ( ext_domain->hasXfemManager() ) {

        Domain *rve_domain = rve->giveDomain(1);

        XfemManager *ext_xMan = ext_domain->giveXfemManager();
        XfemManager *this_xMan = rve->giveDomain(1)->giveXfemManager();
        DynamicDataReader dataReader("fe2");
        if ( ext_xMan != NULL ) {

            std::vector<std::unique_ptr<EnrichmentItem>> eiList;

            //auto &xmanRec = std::make_unique<DynamicInputRecord>();
            //ext_xMan->giveInputRecord(* xmanRec);
            //dataReader.insertInputRecord(DataReader :: IR_xfemManRec, xmanRec);

            // Enrichment items
            int nEI = ext_xMan->giveNumberOfEnrichmentItems();
            for ( int i = 1; i <= nEI; i++ ) {
                EnrichmentItem *ext_ei = ext_xMan->giveEnrichmentItem(i);
                ext_ei->appendInputRecords(dataReader);


                auto &mir = dataReader.giveInputRecord(DataReader :: IR_enrichItemRec, i);
                std :: string name;
                mir.giveRecordKeywordField(name);

                std :: unique_ptr< EnrichmentItem >ei( classFactory.createEnrichmentItem( name.c_str(), i, this_xMan, rve_domain ) );
                if ( ei.get() == NULL ) {
                    OOFEM_ERROR( "unknown enrichment item (%s)", name.c_str() );
                }

                ei->initializeFrom(mir);
                ei->instanciateYourself(dataReader);
                eiList.push_back( std :: move(ei) );

            }

            this_xMan->clearEnrichmentItems();
            this_xMan->appendEnrichmentItems(eiList);

            rve_domain->postInitialize();
            rve->forceEquationNumbering();
        }

    }

    //printf("done.\n");

#if 0
    Domain *newDomain = fe2ms->giveRVE()->giveDomain(1)->Clone();
    newDomain->SetEngngModel(rve.get());
    bool deallocateOld = true;
    rve->setDomain(1, newDomain, deallocateOld);

    //rve->giveDomain(1)->postInitialize();
    rve->giveNumericalMethod(NULL)->setDomain(newDomain);

    rve->postInitialize();
    //rve->forceEquationNumbering();

    rve->initMetaStepAttributes( rve->giveMetaStep(1) );
    rve->giveNextStep(); // Makes sure there is a timestep (which we will modify before solving a step)
    rve->init();


//    std :: ostringstream name;
//    name << this->rve->giveOutputBaseFileName() << "-gp" << n;
//    this->rve->letOutputBaseFileNameBe( name.str() );
//    n++;

    double crackLength = 0.0;
    XfemStructureManager *xMan = dynamic_cast<XfemStructureManager*>( rve->giveDomain(1)->giveXfemManager() );
    if ( xMan ) {
        crackLength = xMan->computeTotalCrackLength();
    }

    std :: ostringstream name;
    name << this->rve->giveOutputBaseFileName() << "-gp" << num << "crackLength" << crackLength;
    if ( this->domain->giveEngngModel()->isParallel() && this->domain->giveEngngModel()->giveNumberOfProcesses() > 1 ) {
        name << "." << this->domain->giveEngngModel()->giveRank();
    }

    num++;

    this->rve->letOutputBaseFileNameBe( name.str() );


    // Update BC
    this->bc = dynamic_cast< PrescribedGradientHomogenization * >( this->rve->giveDomain(1)->giveBc(1) );

#if 1

    XfemSolverInterface *xfemSolInt = dynamic_cast<XfemSolverInterface*>(rve.get());
    StaticStructural *statStruct = dynamic_cast<StaticStructural*>(rve.get());
    if ( xfemSolInt && statStruct ) {
        //printf("Successfully casted to XfemSolverInterface.\n");

        TimeStep *tStep = rve->giveCurrentStep();

        EModelDefaultEquationNumbering num;
        int numDofsNew = rve->giveNumberOfDomainEquations( 1, num );
        FloatArray u;
        u.resize(numDofsNew);
        u.zero();

        xfemSolInt->xfemUpdatePrimaryField(*statStruct, tStep, u);

        // Set domain pointer to various components ...
        rve->giveNumericalMethod(NULL)->setDomain(newDomain);
        //ioEngngModel.nMethod->setDomain(domain);

    }

//    TimeStep *tStep = rve->giveNextStep();
//    setTimeStep(tStep);
//    rve->solveYourselfAt(tStep);


    int numExpModules = rve->giveExportModuleManager()->giveNumberOfModules();
    for ( int i = 1; i <= numExpModules; i++ ) {
        //  ... by diving deep into the hierarchies ... :-/
        VTKXMLExportModule *vtkxmlMod = dynamic_cast< VTKXMLExportModule * >( rve->giveExportModuleManager()->giveModule(i) );
        if ( vtkxmlMod != NULL ) {
            vtkxmlMod->giveSmoother()->setDomain(newDomain);
            vtkxmlMod->givePrimVarSmoother()->setDomain(newDomain);
        }
    }
#endif
#endif
}


} // end namespace oofem
