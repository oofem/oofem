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

// Implementation according paper:
// S.R. Idelsohn, E. Onate, F. Del Pin
// The particle finite element method: a powerful tool to solve
// incompressible flows with free-surface and breaking waves

#include "pfem.h"
#include "nummet.h"
#include "timestep.h"
#include "metastep.h"
#include "element.h"
#include "dofmanager.h"
#include "elementside.h"
#include "dof.h"
#include "initialcondition.h"
#include "verbose.h"
#include "connectivitytable.h"
#include "pfemelement.h"
#include "tr1_2d_pfem.h"
#include "delaunaytriangulator.h"
#include "load.h"
#include "pfemparticle.h"
#include "octreelocalizert.h" //changed from "octreelocalizer.h"
#include "spatiallocalizer.h"
#include "classfactory.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "leplic.h"


namespace oofem {
REGISTER_EngngModel(PFEM);


void PFEMPressureRhsAssembler :: vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const
{
    PFEMElement &pelem = static_cast< PFEMElement & >( element );

    FloatMatrix d;
    pelem.giveCharacteristicMatrix(d, DivergenceMatrix, tStep);
    FloatArray u_star;
    pelem.computeVectorOf(VM_Intermediate, tStep, u_star);
    vec.beProductOf(d, u_star);

    FloatArray reducedVec;
    pelem.computePrescribedRhsVector(reducedVec, tStep, mode);
    reducedVec.negated();
    vec.assemble( reducedVec, pelem.givePressureDofMask() );
}


void PFEMCorrectionRhsAssembler :: vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const
{
    PFEMElement &pelem = static_cast< PFEMElement & >( element );

    FloatMatrix g;
    FloatArray p;
    pelem.computeVectorOfPressures(VM_Total, tStep, p);
    pelem.computeGradientMatrix(g, tStep);
    vec.beProductOf(g, p);
    vec.times(-deltaT);

    FloatMatrix m;
    FloatArray u;
    pelem.computeDiagonalMassMtrx(m, tStep);
    pelem.computeVectorOfVelocities(VM_Intermediate, tStep, u);
    for ( int i = 0; i < m.giveNumberOfColumns(); ++i ) {
        vec [ i ] += m(i, i) * u [ i ];
    }
}

void PFEMCorrectionRhsAssembler :: locationFromElement(IntArray &loc, Element &element, const UnknownNumberingScheme &s, IntArray *dofIds) const
{
    // Note: For 2D, the V_w will just be ignored anyweay, so it's safe to use all three always.
    //element.giveLocationArray(loc, {V_u, V_v, V_w}, s, dofIds);
    element.giveLocationArray(loc, { V_u, V_v }, s, dofIds);
}


void PFEMLaplaceVelocityAssembler :: vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const
{
    PFEMElement &pelem = static_cast< PFEMElement & >( element );

    FloatMatrix l;
    FloatArray u;
    pelem.computeStiffnessMatrix(l, TangentStiffness, tStep);
    pelem.computeVectorOfVelocities(VM_Total, tStep, u);
    vec.beProductOf(l, u);
}

void PFEMLaplaceVelocityAssembler :: locationFromElement(IntArray &loc, Element &element, const UnknownNumberingScheme &s, IntArray *dofIds) const
{
    //element.giveLocationArray(loc, {V_u, V_v, V_w}, s, dofIds);
    element.giveLocationArray(loc, { V_u, V_v }, s, dofIds);
}


void PFEMMassVelocityAssembler :: vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const
{
    PFEMElement &pelem = static_cast< PFEMElement & >( element );

    FloatMatrix m;
    FloatArray u;
    pelem.computeDiagonalMassMtrx(m, tStep);
    pelem.computeVectorOfVelocities(VM_Total, tStep, u);
    vec.beProductOf(m, u);
}

void PFEMMassVelocityAssembler :: locationFromElement(IntArray &loc, Element &element, const UnknownNumberingScheme &s, IntArray *dofIds) const
{
    //element.giveLocationArray(loc, {V_u, V_v, V_w}, s, dofIds);
    element.giveLocationArray(loc, { V_u, V_v }, s, dofIds);
}


void PFEMPressureLaplacianAssembler :: matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const
{
    PFEMElement &pelem = static_cast< PFEMElement & >( element );
    pelem.computePressureLaplacianMatrix(mat, tStep);
}

void PFEMPressureLaplacianAssembler :: locationFromElement(IntArray &loc, Element &element, const UnknownNumberingScheme &s, IntArray *dofIds) const
{
    element.giveLocationArray(loc, { P_f }, s, dofIds);
}


NumericalMethod *PFEM :: giveNumericalMethod(MetaStep *mStep)
{
    if ( nMethod ) {
        return nMethod;
    }

    nMethod = classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this);
    if ( nMethod == NULL ) {
        OOFEM_ERROR("linear solver creation failed for lstype %d", solverType);
    }

    return nMethod;
}


int
PFEM :: forceEquationNumbering(int id)
{
    resetEquationNumberings();
    // forces equation renumbering for current time step
    // intended mainly for problems with changes of static system
    // during solution
    // OUTPUT:
    // sets this->numberOfEquations and this->numberOfPrescribedEquations and returns this value

    // first initialize default numbering (for velocity unknowns only)

    Domain *domain = this->giveDomain(id);
    TimeStep *currStep = this->giveCurrentStep();

    this->domainNeqs.at(id) = 0;
    this->domainPrescribedNeqs.at(id) = 0;

    // number velocities first
    for ( auto &dman : domain->giveDofManagers() ) {
        for ( Dof *jDof: *dman ) {
            DofIDItem type = jDof->giveDofID();
            if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                jDof->askNewEquationNumber(currStep);
            }
        }
    }

    // initialize the independent pressure and auxiliary velocity/equation numbering
    pns.init(domain, currStep);
    avns.init(domain);

    return domainNeqs.at(id);
}



IRResultType
PFEM :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                            // Required by IR_GIVE_FIELD macro

    EngngModel :: initializeFrom(ir);
    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_lstype);
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = ( SparseMtrxType ) val;

    IR_GIVE_FIELD(ir, deltaT, _IFT_PFEM_deltat);
    minDeltaT = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, minDeltaT, _IFT_PFEM_mindeltat);

    alphaShapeCoef = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaShapeCoef, _IFT_PFEM_alphashapecoef);

    maxiter = 50;
    IR_GIVE_OPTIONAL_FIELD(ir, maxiter, _IFT_PFEM_maxiter);

    rtolv = 1.e-8;
    IR_GIVE_OPTIONAL_FIELD(ir, rtolv, _IFT_PFEM_rtolv);

    rtolp = 1.e-8;
    IR_GIVE_OPTIONAL_FIELD(ir, rtolp, _IFT_PFEM_rtolp);

    particleRemovalRatio = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, particleRemovalRatio, _IFT_PFEM_particalRemovalRatio);

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_PFEM_printVolumeReport);
    printVolumeReport = ( val == 1 );

    IR_GIVE_OPTIONAL_FIELD(ir, discretizationScheme, _IFT_PFEM_discretizationScheme);

    IR_GIVE_FIELD(ir, associatedMaterial, _IFT_PFEM_associatedMaterial);
    IR_GIVE_FIELD(ir, associatedCrossSection, _IFT_PFEM_associatedCrossSection);
    IR_GIVE_FIELD(ir, associatedPressureBC, _IFT_PFEM_pressureBC);

    return IRRT_OK;
}



double
PFEM :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like pressure or velocity of dof
{
    if ( mode == VM_Intermediate ) {
        if ( AuxVelocity.isNotEmpty() ) {
            int index = avns.giveDofEquationNumber(dof);
            if ( index > 0 && index <= AuxVelocity.giveSize() ) {
                return AuxVelocity.at(index);
            } else {
                return 0.0;
            }
        } else {
            return 0.;
        }
    } else {
        if ( this->requiresUnknownsDictionaryUpdate() ) {
            int hash = this->giveUnknownDictHashIndx(mode, tStep);
            if ( dof->giveUnknowns()->includes(hash) ) {
                return dof->giveUnknowns()->at(hash);
            } else {
                OOFEM_ERROR( "giveUnknown:  Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode) );
                return 0.; // to make compiler happy
            }
        } else {
            return 0.0;
        }
    }
}


TimeStep *
PFEM :: giveSolutionStepWhenIcApply()
{
    if ( !stepWhenIcApply ) {
        stepWhenIcApply.reset( new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0, 0.0, deltaT, 0) );
    }

    return stepWhenIcApply.get();
}

void
PFEM :: preInitializeNextStep()
{
    Domain *domain = this->giveDomain(1);
    domain->clearElements();

    DelaunayTriangulator myMesher(domain, alphaShapeCoef);
    myMesher.generateMesh();

    for ( auto &dman : domain->giveDofManagers() ) {
        PFEMParticle *particle = dynamic_cast< PFEMParticle * >( dman.get() );
        particle->setFree();
    }

    if ( domain->giveNumberOfElements() > 0 ) {
        for ( auto &element : domain->giveElements() ) {
            element->checkConsistency();

            for ( int j = 1; j <= element->giveNumberOfDofManagers(); j++ ) {
                dynamic_cast< PFEMParticle * >( element->giveDofManager(j) )->setFree(false);
            }
        }
    } else {
        VERBOSE_PRINTS("Mesh generation failed", "0 elements created");
    }

    for ( auto &dman : domain->giveDofManagers() ) {
        for ( Dof *jDof: *dman ) {
            DofIDItem type = jDof->giveDofID();
            if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                if ( jDof->giveBcId() ) {
                    dynamic_cast< PFEMParticle * >( dman.get() )->setFree(false);
                    break;
                }
            }
        }
    }
}

TimeStep *
PFEM :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();
    double totalTime = 0;
    StateCounterType counter = 1;
    Domain *domain = this->giveDomain(1);

    if ( !currentStep ) {
        // first step -> generate initial step
        currentStep.reset( new TimeStep( * giveSolutionStepWhenIcApply() ) );
    } else {
        istep = currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = std :: move(currentStep);

    EngngModel :: forceEquationNumbering();

    double volume = 0.0;

    double ndt = dynamic_cast< PFEMElement * >( domain->giveElement(1) )->computeCriticalTimeStep( previousStep.get() );
    // check for critical time step
    for ( int i = 2; i <= domain->giveNumberOfElements(); i++ ) {
        TR1_2D_PFEM *ielem = dynamic_cast< TR1_2D_PFEM * >( domain->giveElement(i) );
        double idt = ielem->computeCriticalTimeStep( previousStep.get() );
        if ( idt < ndt ) {
            ndt = idt;
            OOFEM_LOG_INFO("Reducing time step due to element #%i \n", i);
        }
        volume += ielem->computeArea();
    }

    ndt = min(ndt, deltaT);

    totalTime = previousStep->giveTargetTime() + ndt;

    currentStep.reset( new TimeStep(istep, this, 1, totalTime, ndt, counter) );
    // time and dt variables are set eq to 0 for staics - has no meaning

    OOFEM_LOG_INFO( "SolutionStep %d : t = %e, dt = %e\n", istep, totalTime * this->giveVariableScale(VST_Time), ndt * this->giveVariableScale(VST_Time) );
    if ( printVolumeReport ) {
        OOFEM_LOG_INFO("Volume leakage: %.3f%%\n", ( 1.0 - ( volume / domainVolume ) ) * 100.0);
    }

    return currentStep.get();
}

void
PFEM :: solveYourselfAt(TimeStep *tStep)
{
    int auxmomneq = this->giveNumberOfDomainEquations(1, avns);
    int momneq =  this->giveNumberOfDomainEquations(1, vns);
    int presneq =  this->giveNumberOfDomainEquations(1, pns);

    double deltaT = tStep->giveTimeIncrement();

    FloatArray rhs(auxmomneq);


    double d_pnorm = 1.0;
    double d_vnorm = 1.0;

    AuxVelocity.resize(auxmomneq);

    avLhs.clear();
    avLhs.resize(auxmomneq);

    this->assembleVector( avLhs, tStep, LumpedMassVectorAssembler(), VM_Total, avns, this->giveDomain(1) );

    pLhs.reset( classFactory.createSparseMtrx(sparseMtrxType) );
    if ( !pLhs ) {
        OOFEM_ERROR("solveYourselfAt: sparse matrix creation failed");
    }

    pLhs->buildInternalStructure(this, 1, pns);

    this->assemble( * pLhs, tStep, PFEMPressureLaplacianAssembler(), pns, this->giveDomain(1) );

    pLhs->times(deltaT);

    //////////////////////////////////////////////////////////////////////////

    vLhs.clear();
    vLhs.resize(momneq);

    this->assembleVector( vLhs, tStep, LumpedMassVectorAssembler(), VM_Total, vns, this->giveDomain(1) );

    if ( tStep->giveNumber() == giveNumberOfFirstStep() ) {
        TimeStep *stepWhenIcApply = tStep->givePreviousStep();
        this->applyIC(stepWhenIcApply);
    }

    if ( VelocityField.giveActualStepNumber() < tStep->giveNumber() ) {
        VelocityField.advanceSolution(tStep);
    }

    if ( PressureField.giveActualStepNumber() < tStep->giveNumber() ) {
        PressureField.advanceSolution(tStep);
    }

    FloatArray *velocityVector = VelocityField.giveSolutionVector(tStep);
    FloatArray *pressureVector = PressureField.giveSolutionVector(tStep);

    FloatArray velocityVectorThisStep(* velocityVector);
    FloatArray velocityVectorLastStep(* velocityVector);

    FloatArray pressureVectorThisStep(* pressureVector);
    FloatArray pressureVectorLastStep(* pressureVector);

    if ( tStep->isTheFirstStep() ) {
        for ( int i = 1; i <= this->giveDomain(1)->giveNumberOfDofManagers(); i++ ) {
            this->updateDofUnknownsDictionary(this->giveDomain(1)->giveDofManager(i), tStep);
        }
    }

    FloatArray externalForces;
    externalForces.resize(auxmomneq);
    this->assembleVector( externalForces, tStep, ExternalForceAssembler(), VM_Total, avns, this->giveDomain(1) );

    int iteration = 0;
    do {
        iteration++;

        velocityVectorLastStep = velocityVectorThisStep;
        pressureVectorLastStep = pressureVectorThisStep;

        /************************** STEP 1 - calculates auxiliary velocities *************************/

        rhs.resize(auxmomneq);
        rhs.zero();

        if ( discretizationScheme == 1 ) { //implicit
            if ( iteration > 1 ) {
                this->assembleVectorFromElements( rhs, tStep, PFEMLaplaceVelocityAssembler(), VM_Total, avns, this->giveDomain(1) );
            }
        } else if ( discretizationScheme == 0 ) { // explicit
            this->assembleVectorFromElements( rhs, tStep->givePreviousStep(), PFEMLaplaceVelocityAssembler(), VM_Total, avns, this->giveDomain(1) );
        }

        rhs.negated();

        // constant member placed outside of the loop
        rhs.add(externalForces);

        rhs.times(deltaT);


        /// assumming the mass matrix is always lumped, PFEMMassVelocityAssembler could be replaced by:
        // for ( int i = 1; i <= momneq; i++ ) {
        //     rhs.at(i) +=  avLhs.at(i) * oldVelocityVector.at(i);
        // }

        if ( tStep->isTheFirstStep() == false ) {
            this->assembleVectorFromElements( rhs, tStep->givePreviousStep(), PFEMMassVelocityAssembler(), VM_Total, avns, this->giveDomain(1) );
        }

        this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );
        AuxVelocity.resize( rhs.giveSize() );
        for ( int i = 1; i <= auxmomneq; i++ ) {
            AuxVelocity.at(i) = rhs.at(i) / avLhs.at(i);
        }

        /************************* STEP 2 - calculates pressure ************************************/

        rhs.resize(presneq);
        rhs.zero();

        this->assembleVectorFromElements( rhs, tStep, PFEMPressureRhsAssembler(), VM_Total, pns, this->giveDomain(1) );

        this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );
        pressureVector->resize(presneq);
        nMethod->solve(* pLhs, rhs, * pressureVector);

        for ( auto &dman : this->giveDomain(1)->giveDofManagers() ) {
            this->updateDofUnknownsDictionaryPressure(dman.get(), tStep);
        }

        /**************************** STEP 3 - velocity correction step ********************************/



        rhs.resize(momneq);
        rhs.zero();
        velocityVector->resize(momneq);

        this->assembleVectorFromElements( rhs, tStep, PFEMCorrectionRhsAssembler(deltaT), VM_Total, vns, this->giveDomain(1) );

        this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );

        velocityVector->resize( rhs.giveSize() );
        for ( int i = 1; i <= momneq; i++ ) {
            velocityVector->at(i) = rhs.at(i) / vLhs.at(i);
        }

        for ( int i = 1; i <= this->giveDomain(1)->giveNumberOfDofManagers(); i++ ) {
            PFEMParticle *particle = dynamic_cast< PFEMParticle * >( this->giveDomain(1)->giveDofManager(i) );
            this->updateDofUnknownsDictionaryVelocities(particle, tStep);
        }

        velocityVectorThisStep = * velocityVector;
        pressureVectorThisStep = * pressureVector;

        if ( iteration == 1 ) {
            d_vnorm = velocityVectorThisStep.computeNorm();
            d_pnorm = pressureVectorThisStep.computeNorm();
        } else {
            FloatArray diffVelocity = velocityVectorThisStep;
            diffVelocity.subtract(velocityVectorLastStep);

            d_vnorm = 0.0;
            for ( int i = 1; i <= diffVelocity.giveSize(); i++ ) {
                d_vnorm = max(d_vnorm, fabs( diffVelocity.at(i) ) > 1.e-6 ? fabs( diffVelocity.at(i) / velocityVectorLastStep.at(i) ) : 0);
            }

            FloatArray diffPressure = pressureVectorThisStep;
            diffPressure.subtract(pressureVectorLastStep);
            d_pnorm = 0.0;
            for ( int i = 1; i <= diffPressure.giveSize(); i++ ) {
                // TODO : what about divide by zero ????
                d_pnorm = max(d_pnorm, fabs( diffPressure.at(i) ) > 1.e-6 ? fabs( diffPressure.at(i) / pressureVectorLastStep.at(i) ) : 0);
            }
        }
    } while ( discretizationScheme == 1 && ( d_vnorm > rtolv || d_pnorm > rtolp ) && iteration < maxiter );

    if ( iteration > maxiter ) {
        OOFEM_ERROR("Maximal iteration count exceded");
    } else {
        OOFEM_LOG_INFO("\n %i iterations performed.\n", iteration);
    }



    Domain *d = this->giveDomain(1);
    for ( auto &dman : d->giveDofManagers() ) {
        PFEMParticle *particle = dynamic_cast< PFEMParticle * >( dman.get() );

        if ( particle->isFree() && particle->isActive() ) {
            for ( Dof *dof : *particle ) {
                DofIDItem type = dof->giveDofID();
                if ( type == V_u || type == V_v || type == V_w ) {
                    int eqnum = dof->giveEquationNumber(vns);
                    if ( eqnum ) {
                        double previousValue = dof->giveUnknown( VM_Total, tStep->givePreviousStep() );
                        Load *load = d->giveLoad(3);
                        FloatArray gVector;
                        load->computeComponentArrayAt(gVector, tStep, VM_Total);
                        previousValue += gVector [ type - V_w ] * deltaT; // Get the corresponding coordinate index for the velocities.
                        velocityVector->at(eqnum) = previousValue;
                    }
                }
            }
        }
    }
    tStep->incrementStateCounter();
}


void
PFEM :: updateYourself(TimeStep *stepN)
{
    this->updateInternalState(stepN);
    EngngModel :: updateYourself(stepN);

    this->deactivateTooCloseParticles();
}



void
PFEM :: updateInternalState(TimeStep *stepN)
{
    for ( auto &domain : this->domainList ) {
        if ( requiresUnknownsDictionaryUpdate() ) {
            for ( auto &dman : domain->giveDofManagers() ) {
                this->updateDofUnknownsDictionary(dman.get(), stepN);
            }
        }

        for ( auto &elem : domain->giveElements() ) {
            elem->updateInternalState(stepN);
        }
    }
}

int
PFEM :: giveUnknownDictHashIndx(ValueModeType mode, TimeStep *stepN)
{
    if ( ( stepN == this->giveCurrentStep() ) || ( stepN == this->givePreviousStep() ) ) {
        int index = ( stepN->giveNumber() % 2 ) * 100 + mode;
        return index;
    } else {
        OOFEM_ERROR("giveUnknownDictHashIndx: unsupported solution step");
    }

    return 0;
}

void
PFEM :: updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep)
{
    for ( Dof *iDof: *inode ) {
        int eqNum = 0;
        double val = 0;
        if ( iDof->hasBc(tStep) ) {                 // boundary condition
            val = iDof->giveBcValue(VM_Total, tStep);
        } else {
            if ( iDof->giveDofID() == P_f ) {
                eqNum = pns.giveDofEquationNumber(iDof);
                FloatArray *vect = PressureField.giveSolutionVector(tStep);
                if ( vect->giveSize() > 0 ) {                     // in the first step -> zero will be set
                    val = vect->at(eqNum);
                }
            } else {                 // velocities
                eqNum = iDof->__giveEquationNumber();
                FloatArray *vect = VelocityField.giveSolutionVector(tStep);
                if ( vect->giveSize() > 0 ) {
                    val = vect->at(eqNum);
                }
            }
        }

        iDof->updateUnknownsDictionary(tStep, VM_Total, val);
    }
}

void
PFEM :: updateDofUnknownsDictionaryPressure(DofManager *inode, TimeStep *tStep)
{
    Dof *iDof = inode->giveDofWithID(P_f);
    double val = 0;
    if ( iDof->hasBc(tStep) ) {
        val = iDof->giveBcValue(VM_Total, tStep);
    } else {
        int eqNum = pns.giveDofEquationNumber(iDof);
        FloatArray *vect = PressureField.giveSolutionVector(tStep);
        val = vect->at(eqNum);
    }

    iDof->updateUnknownsDictionary(tStep, VM_Total, val);
}

void
PFEM :: updateDofUnknownsDictionaryVelocities(DofManager *inode, TimeStep *tStep)
{
    for ( Dof *iDof : *inode ) {
        DofIDItem type = iDof->giveDofID();
        if ( type == V_u || type == V_v || type == V_w ) {
            int eqNum = 0;
            double val = 0;
            if ( iDof->hasBc(tStep) ) {             // boundary condition
                val = iDof->giveBcValue(VM_Total, tStep);
            } else {
                eqNum = iDof->__giveEquationNumber();
                FloatArray *vect = VelocityField.giveSolutionVector(tStep);
                if ( vect->giveSize() > 0 ) {
                    val = vect->at(eqNum);
                }
            }
            iDof->updateUnknownsDictionary(tStep, VM_Total, val);
        }
    }
}


// NOT ACTIVE
contextIOResultType
PFEM :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves state variable - displacement vector
//
{
    contextIOResultType iores;
    int closeFlag = 0;
    FILE *file = NULL;

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, this->giveCurrentStep()->giveNumber(),
                                    this->giveCurrentStep()->giveVersion(), contextMode_write) ) {
            THROW_CIOERR(CIO_IOERR);                             // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = EngngModel :: saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = PressureField.saveContext(* stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = VelocityField.saveContext(* stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }

    return CIO_OK;
}


// NOT ACTIVE
contextIOResultType
PFEM :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restore state variable - displacement vector
//
{
    contextIOResultType iores;
    int closeFlag = 0;
    int istep, iversion;
    FILE *file = NULL;

    this->resolveCorrespondingStepNumber(istep, iversion, obj);

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, istep, iversion, contextMode_read) ) {
            THROW_CIOERR(CIO_IOERR);                     // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = EngngModel :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = PressureField.restoreContext(* stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = VelocityField.restoreContext(* stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }

    return CIO_OK;
}


int
PFEM :: checkConsistency()
{
    Domain *domain = this->giveDomain(1);

    // check for proper element type
    for ( auto &elem : domain->giveElements() ) {
        if ( !dynamic_cast< PFEMElement * >( elem.get() ) ) {
            OOFEM_WARNING( "Element %d has no PFEM base", elem->giveLabel() );
            return 0;
        }
    }

    return EngngModel :: checkConsistency();
}


void
PFEM :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    DofIDItem type = iDof->giveDofID();
    if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
        iDof->printSingleOutputAt(stream, atTime, '*', VM_Intermediate);
        iDof->printSingleOutputAt(stream, atTime, 'u', VM_Total);

        // printing coordinate in DofMan Output
        DofManager *dman = iDof->giveDofManager();
        double coordinate = 0.0;
        int dofNumber = 0;
        switch ( type ) {
        case V_u: dofNumber = 1;
            break;
        case V_v: dofNumber = 2;
            break;
        case V_w: dofNumber = 3;
            break;
        default:;
        }
        coordinate = dman->giveCoordinate(dofNumber);
        fprintf(stream, "  dof %d   c % .8e\n", dofNumber, coordinate);
    } else if ( ( type == P_f ) ) {
        iDof->printSingleOutputAt(stream, atTime, 'p', VM_Total);
    } else {
        OOFEM_ERROR("printDofOutputAt: unsupported dof type");
    }
}

void
PFEM :: applyIC(TimeStep *stepWhenIcApply)
{
    Domain *domain = this->giveDomain(1);
    int mbneq =  this->giveNumberOfDomainEquations(1, vns);
    int pdneq =  this->giveNumberOfDomainEquations(1, pns);
    FloatArray *velocityVector, *pressureVector;

#ifdef VERBOSE
    OOFEM_LOG_INFO("Applying initial conditions\n");
#endif

    int velocityFieldStepNumber = VelocityField.giveActualStepNumber();

    if ( velocityFieldStepNumber < stepWhenIcApply->giveNumber() ) {
        VelocityField.advanceSolution(stepWhenIcApply);
    }
    velocityVector = VelocityField.giveSolutionVector(stepWhenIcApply);
    velocityVector->resize(mbneq);
    velocityVector->zero();

    int pressureFieldStepNumber = PressureField.giveActualStepNumber();
    if ( pressureFieldStepNumber < stepWhenIcApply->giveNumber() ) {
        PressureField.advanceSolution(stepWhenIcApply);
    }
    pressureVector = PressureField.giveSolutionVector(stepWhenIcApply);
    pressureVector->resize(pdneq);
    pressureVector->zero();


    for ( auto &node : domain->giveDofManagers() ) {
        for ( Dof *iDof: *node ) {
            // ask for initial values obtained from
            // bc (boundary conditions) and ic (initial conditions)
            //iDof  =  node->giveDof(k);
            if ( !iDof->isPrimaryDof() ) {
                continue;
            }

            int jj = iDof->__giveEquationNumber();
            DofIDItem type = iDof->giveDofID();

            if ( jj ) {
                if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                    velocityVector->at(jj) = iDof->giveUnknown(VM_Total, stepWhenIcApply);
                } else {
                    pressureVector->at(jj) = iDof->giveUnknown(VM_Total, stepWhenIcApply);
                }
            }
        }
    }

    // update element state according to given ic
    for ( auto &elem : domain->giveElements() ) {
        PFEMElement *element = static_cast< PFEMElement * >( elem.get() );
        element->updateInternalState(stepWhenIcApply);
        element->updateYourself(stepWhenIcApply);
        domainVolume += element->computeArea();
    }
}

int
PFEM :: giveNewEquationNumber(int domain, DofIDItem id)
{
    if ( ( id == V_u ) || ( id == V_v ) || ( id == V_w ) ) {
        return this->vns.askNewEquationNumber();
    } else {
        OOFEM_ERROR("giveNewEquationNumber:: Unknown DofIDItem");
    }

    return 0;
}

int
PFEM :: giveNewPrescribedEquationNumber(int domain, DofIDItem id)
{
    if ( ( id == V_u ) || ( id == V_v ) || ( id == V_w ) ) {
        return prescribedVns.askNewEquationNumber();
    } else {
        OOFEM_ERROR("giveNewPrescribedEquationNumber:: Unknown DofIDItem");
    }

    return 0;
}

int PFEM :: giveNumberOfDomainEquations(int d, const UnknownNumberingScheme &numberingScheme) {            //
    // returns number of equations of current problem
    // this method is implemented here, because some method may add some
    // conditions in to system and this may results into increased number of
    // equations.
    //

    if ( !equationNumberingCompleted ) {
        EngngModel :: forceEquationNumbering();
    }

    return numberingScheme.giveRequiredNumberOfDomainEquation();
}


void
PFEM :: resetEquationNumberings()
{
    vns.reset();
    prescribedVns.reset();
    pns.reset();
    avns.reset();
}

void
PFEM :: deactivateTooCloseParticles()
{
    Domain *d = this->giveDomain(1);
    // deactivating particles
    if ( particleRemovalRatio > 1.e-6 ) { // >0
        for ( int i = 1; i <= d->giveNumberOfElements(); i++ ) {
            TR1_2D_PFEM *element = dynamic_cast< TR1_2D_PFEM * >( d->giveElement(i) );

            PFEMParticle *particle1 = dynamic_cast< PFEMParticle * >( element->giveNode(1) );
            PFEMParticle *particle2 = dynamic_cast< PFEMParticle * >( element->giveNode(2) );
            PFEMParticle *particle3 = dynamic_cast< PFEMParticle * >( element->giveNode(3) );

            double l12 = particle1->giveCoordinates()->distance( particle2->giveCoordinates() );
            double l23 = particle2->giveCoordinates()->distance( particle3->giveCoordinates() );
            double l31 = particle3->giveCoordinates()->distance( particle1->giveCoordinates() );

            double maxLength = max( l12, max(l23, l31) );
            double minLength = min( l12, min(l23, l31) );

            if ( minLength / maxLength < particleRemovalRatio ) {
                if ( fabs(l12 - minLength) < 1.e-6 ) {             // ==
                    if ( particle1->isActive() && particle2->isActive() ) {
                        bool isSupported = false;

                        for ( Dof *jDof: *particle1 ) {//int j = 1; j <= particle1->giveNumberOfDofs();
                            DofIDItem type  =  jDof->giveDofID();
                            if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                                if ( jDof->giveBcId() ) {
                                    isSupported = true;
                                }
                            }
                        }
                        if ( isSupported == false ) {
                            particle1->deactivate();
                        } else {
                            particle2->deactivate();
                        }
                    } else if ( particle1->isActive() == false || particle2->isActive() == false ) {
                        ;                        // it was deactivated by other element
                    } else {
                        OOFEM_ERROR("Both particles deactivated");
                    }
                }

                if ( fabs(l23 - minLength) < 1.e-6 ) {             // ==
                    if ( particle2->isActive() && particle3->isActive() ) {
                        bool isSupported = false;

                        for ( Dof *jDof  : *particle2 ) {
                            DofIDItem type  =  jDof->giveDofID();
                            if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                                if ( jDof->giveBcId() ) {
                                    isSupported = true;
                                }
                            }
                        }
                        if ( isSupported == false ) {
                            particle2->deactivate();
                        } else {
                            particle3->deactivate();
                        }
                    } else if ( particle2->isActive() == false || particle3->isActive() == false ) {
                        ;                        // it was deactivated by other element
                    } else {
                        OOFEM_ERROR("Both particles deactivated");
                    }
                }

                if ( fabs(l31 - minLength) < 1.e-6 ) {             // ==
                    if ( particle3->isActive() && particle1->isActive() ) {
                        bool isSupported = false;
                        for ( Dof *jDof  : *particle3 ) {
                            DofIDItem type  =  jDof->giveDofID();
                            if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                                if ( jDof->giveBcId() ) {
                                    isSupported = true;
                                }
                            }
                        }
                        if ( isSupported == false ) {
                            particle3->deactivate();
                        } else {
                            particle1->deactivate();
                        }
                    } else if ( particle3->isActive() == false || particle1->isActive() == false ) {
                        ;                        // it was deactivated by other element
                    } else {
                        OOFEM_ERROR("Both particles deactivated");
                    }
                }
            }
        }
    }
}
} // end namespace oofem
