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


#include "mpmsemodels.h"
#include "classfactory.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "dof.h"
#include "dictionary.h"
#include "verbose.h"
#include "mathfem.h"
#include "assemblercallback.h"
#include "unknownnumberingscheme.h"
#include "dofdistributedprimaryfield.h"
#include "primaryfield.h"
#include "maskedprimaryfield.h"
#include "nrsolver.h"
#include "activebc.h"
#include "boundarycondition.h"
#include "boundaryload.h"
#include "outputmanager.h"
#include "mpm.h"
#include "mpmassemblers.h"

namespace oofem {
    REGISTER_EngngModel(StationaryMPMSProblem);   
    REGISTER_EngngModel(NonStationaryMPMSProblem); 
    REGISTER_EngngModel_Alt(NonStationaryMPMSProblem, mpmproblem); // for compatibility with old input files

void
NonStationaryMPMSProblem :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    EngngModel :: initializeFrom(ir);

    int val = SMT_Skyline;


    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    this->sparseMtrxType = ( SparseMtrxType ) val;

    if ( ir->hasField("initt") ) {
        IR_GIVE_FIELD(ir, initT, "initt");
    }

    if ( ir->hasField("deltat") ) {
        IR_GIVE_FIELD(ir, deltaT, "deltat");
    } else if ( ir->hasField("deltatfunction") ) {
        IR_GIVE_FIELD(ir, dtFunction, "deltatfunction");
    } else if ( ir->hasField("prescribedtimes") ) {
        IR_GIVE_FIELD(ir, prescribedTimes, "prescribedtimes");
    } else {
        throw ValueInputException(ir, "none", "Time step not defined");
    }

    IR_GIVE_FIELD(ir, alpha, "alpha");
    problemType = "up"; // compatibility mode @TODO Remove default value
    IR_GIVE_OPTIONAL_FIELD(ir, problemType, "ptype");
    if (!((problemType == "up")||(problemType == "tm")||(problemType == "symbolic"))) {
      throw ValueInputException(ir, "none", "Problem type not recognized");
    }
    OOFEM_LOG_RELEVANT("MPM: %s formulation\n", problemType.c_str());

    if (problemType == "symbolic") {
        IR_GIVE_FIELD(ir, lhsIntegrals, "lhsterms");
        IR_GIVE_FIELD(ir, lhsdotIntegrals, "lhsdotterms");
        IR_GIVE_FIELD(ir, rhsIntegrals, "rhsterms");
    }
    
    this->keepTangent = ir->hasField("keeptangent");
    field = std::make_unique<DofDistributedPrimaryField>(this, 1, FT_TransportProblemUnknowns, 2, this->alpha);

    // read field export flag
    exportFields.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, exportFields, "exportfields");
    if ( exportFields.giveSize() ) {
        FieldManager *fm = this->giveContext()->giveFieldManager();
        for ( int i = 1; i <= exportFields.giveSize(); i++ ) {
            if ( exportFields.at(i) == FT_Displacements ) {
                FieldPtr _displacementField( new MaskedPrimaryField ( ( FieldType ) exportFields.at(i), this->field.get(), {D_u, D_v, D_w} ) );
                fm->registerField( _displacementField, ( FieldType ) exportFields.at(i) );
            } else if ( exportFields.at(i) == FT_Pressure ) {
                FieldPtr _pressureField( new MaskedPrimaryField ( ( FieldType ) exportFields.at(i), this->field.get(), {P_f} ) );
                fm->registerField( _pressureField, ( FieldType ) exportFields.at(i) );
            }
        }
    }
    IR_GIVE_OPTIONAL_FIELD(ir, this->debug, "debug");
}

double NonStationaryMPMSProblem :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    return this->field->giveUnknownValue(dof, mode, tStep);
}


double
NonStationaryMPMSProblem :: giveDeltaT(int n)
{
    if ( this->dtFunction ) {
        return this->giveDomain(1)->giveFunction(this->dtFunction)->evaluateAtTime(n);
    } else if ( this->prescribedTimes.giveSize() > 0 ) {
        return this->giveDiscreteTime(n) - this->giveDiscreteTime(n - 1);
    } else {
        return this->deltaT;
    }
}

double
NonStationaryMPMSProblem :: giveDiscreteTime(int iStep)
{
    if ( iStep > 0 && iStep <= this->prescribedTimes.giveSize() ) {
        return ( this->prescribedTimes.at(iStep) );
    } else if ( iStep == 0 ) {
        return initT;
    }

    OOFEM_ERROR("invalid iStep");
    return 0.0;
}


TimeStep *NonStationaryMPMSProblem :: giveNextStep()
{
    return this->timeStepController->giveNextStep();
}

TimeStep *NonStationaryMPMSProblem :: giveSolutionStepWhenIcApply(bool force)
{
    if ( master && (!force)) {
        return master->giveSolutionStepWhenIcApply();
    } else {
        if ( !stepWhenIcApply ) {
            double dt = this->giveDeltaT(1);
            stepWhenIcApply = std::make_unique<TimeStep>(giveNumberOfTimeStepWhenIcApply(), this, 0, this->initT, dt, 0);
            // The initial step goes from [-dt, 0], so the intrinsic time is at: -deltaT  + alpha*dt
            stepWhenIcApply->setIntrinsicTime(-dt + alpha * dt);
        }

        return stepWhenIcApply.get();
    }
}

void NonStationaryMPMSProblem :: solveYourselfAt(TimeStep *tStep)
{
    OOFEM_LOG_INFO( "\nSolving [step number %5d, time %e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
    
    Domain *d = this->giveDomain(1);
    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );

    if ( tStep->isTheFirstStep() ) {
        this->applyIC();
    }

    field->advanceSolution(tStep);
    field->initialize(VM_Total, tStep, solution, EModelDefaultEquationNumbering());

    if ( !effectiveMatrix ) {
        effectiveMatrix = classFactory.createSparseMtrx(sparseMtrxType);
        effectiveMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
    }

    OOFEM_LOG_INFO("Assembling external forces\n");
    FloatArray externalForces(neq);
    externalForces.zero();
    if (this->problemType == "symbolic") {
        // loop over rhs integrals
        for (auto i: rhsIntegrals) {
            Integral* integral = this->integralList[i-1].get();
            integral->assemble_rhs (externalForces, EModelDefaultEquationNumbering(), tStep); 
        }
        // experimental: allow traditional BCs to contribute to rhs as well by treating them as integrals
        this->assembleVectorFromBC(externalForces, tStep, ExternalForceAssembler(), VM_Total,
                                    EModelDefaultEquationNumbering(), this->giveDomain(1) );
    } else {
        this->assembleVector( externalForces, tStep, ExternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), d );
    }
    this->updateSharedDofManagers(externalForces, EModelDefaultEquationNumbering(), LoadExchangeTag);

    // set-up numerical method
    this->giveNumericalMethod( this->giveCurrentMetaStep() );
    OOFEM_LOG_INFO("Solving for %d unknowns...\n", neq);

    internalForces.resize(neq);

    FloatArray incrementOfSolution;
    double loadLevel;
    int currentIterations = 0;
    this->updateInternalRHS(this->internalForces, tStep, this->giveDomain(1), NULL); /// @todo Hack to ensure that internal RHS is evaluated before the tangent. This is not ideal, causing this to be evaluated twice for a linearproblem. We have to find a better way to handle this.
    ConvergedReason status = this->nMethod->solve(*this->effectiveMatrix,
                                                  externalForces,
                                                  nullptr, // ignore
                                                  this->solution,
                                                  incrementOfSolution,
                                                  this->internalForces,
                                                  this->eNorm,
                                                  loadLevel, // ignore
                                                  SparseNonLinearSystemNM :: rlm_total, // ignore
                                                  currentIterations, // ignore
                                                  tStep);
    if (status != CR_CONVERGED) {
      OOFEM_WARNING("No success in solving problem at step %d", tStep->giveNumber());
    }

    tStep->numberOfIterations = currentIterations;
    tStep->convergedReason = status;
}


void
NonStationaryMPMSProblem :: updateSolution(FloatArray &solutionVector, TimeStep *tStep, Domain *d)
{
    ///@todo NRSolver should report when the solution changes instead of doing it this way.
    this->field->update(VM_Total, tStep, solutionVector, EModelDefaultEquationNumbering());
    ///@todo Need to reset the boundary conditions properly since some "update" is doing strange
    /// things such as applying the (wrong) boundary conditions. This call will be removed when that code can be removed.
    this->field->applyBoundaryCondition(tStep);
}


void
NonStationaryMPMSProblem :: updateInternalRHS(FloatArray &answer, TimeStep *tStep, Domain *d, FloatArray *eNorm)
{
    if ( eNorm ) {
        int maxdofids = this->giveDomain(1)->giveMaxDofID();
#ifdef __MPI_PARALLEL_MODE
        if ( this->isParallel() ) {
            int val;
            MPI_Allreduce(& maxdofids, & val, 1, MPI_INT, MPI_MAX, this->comm);
            maxdofids = val;
        }
#endif
        eNorm->resize(maxdofids);
        eNorm->zero();
    }
    // F_eff = F(T^(k)) + C * dT/dt^(k)
    answer.zero();
    if (this->problemType == "up") {
      this->assembleVector(answer, tStep, UPResidualAssembler(this->alpha, tStep->giveTimeIncrement()), VM_Total, EModelDefaultEquationNumbering(), d, eNorm);
    } else if (this->problemType == "tm") {
      this->assembleVector(answer, tStep, TMResidualAssembler(this->alpha, tStep->giveTimeIncrement()), VM_Total, EModelDefaultEquationNumbering(), d, eNorm);
    } else if (this->problemType == "symbolic") {
        for (auto i: lhsIntegrals) {
            Integral* integral = this->integralList[i-1].get();
            integral->assemble_rhs (answer, EModelDefaultEquationNumbering(), tStep); 
        }
        for (auto i: lhsdotIntegrals) {
            Integral* integral = this->integralList[i-1].get();
            integral->assemble_rhs (answer, EModelDefaultEquationNumbering(), tStep, eNorm); 
        }
    } else {
      OOFEM_ERROR ("unsupported problemType");
    }
    this->updateSharedDofManagers(answer, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);
}


void
NonStationaryMPMSProblem :: updateMatrix(SparseMtrx &mat, TimeStep *tStep, Domain *d)
{
    // K_eff = (a*K + C/dt)
    if ( !this->keepTangent || !this->hasTangent ) {
        mat.zero();
        if (this->problemType == "up") {
          UPLhsAssembler jacobianAssembler(this->alpha, tStep->giveTimeIncrement());
          //Assembling left hand side 
          this->assemble( *effectiveMatrix, tStep, jacobianAssembler,
                          EModelDefaultEquationNumbering(), d );
        } else if (this->problemType == "tm") {
            TMLhsAssembler jacobianAssembler(this->alpha, tStep->giveTimeIncrement());
            //Assembling left hand side 
            this->assemble( *effectiveMatrix, tStep, jacobianAssembler,
                            EModelDefaultEquationNumbering(), d );
        } else if (this->problemType == "symbolic") {
            for (auto i: lhsIntegrals) {
                Integral* integral = this->integralList[i-1].get();
                integral->assemble_lhs (mat, EModelDefaultEquationNumbering(), tStep, this->alpha); 
            }
            for (auto i: lhsdotIntegrals) {
                Integral* integral = this->integralList[i-1].get();
                integral->assemble_lhs (mat, EModelDefaultEquationNumbering(), tStep, 1.0/tStep->giveTimeIncrement()); 
            }
        } else {
          OOFEM_ERROR ("unsupported problemType");
        }
        
        this->hasTangent = true;
    }
}


void
NonStationaryMPMSProblem :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
{
    ///@todo NRSolver should report when the solution changes instead of doing it this way.
    this->field->update(VM_Total, tStep, solution, EModelDefaultEquationNumbering());
    ///@todo Need to reset the boundary conditions properly since some "update" is doing strange
    /// things such as applying the (wrong) boundary conditions. This call will be removed when that code can be removed.
    this->field->applyBoundaryCondition(tStep);

    if ( cmpn == InternalRhs ) {
        // F_eff = F(T^(k)) + C * dT/dt^(k)
        this->internalForces.zero();
        this->updateInternalRHS(this->internalForces, tStep, d, &this->eNorm);
    } else if ( cmpn == NonLinearLhs ) {
        this->updateMatrix(*this->effectiveMatrix, tStep, d);
    } else {
        OOFEM_ERROR("Unknown component");
    }
}


void
NonStationaryMPMSProblem :: applyIC()
{
    Domain *domain = this->giveDomain(1);
    OOFEM_LOG_INFO("Applying initial conditions\n");

    this->field->applyDefaultInitialCondition();

    // set initial field IP values (needed by some nonlinear materials)
    TimeStep *s = this->giveSolutionStepWhenIcApply();
    for ( auto &elem : domain->giveElements() ) {
        Element *element = elem.get();
        element->updateInternalState(s);
        element->updateYourself(s);
    }
}


bool
NonStationaryMPMSProblem :: requiresEquationRenumbering(TimeStep *tStep)
{
    ///@todo This method should be set as the default behavior instead of relying on a user specified flag. Then this function should be removed.
    if ( tStep->isTheFirstStep() ) {
        return true;
    }
    // Check if Dirichlet b.c.s has changed.
    Domain *d = this->giveDomain(1);
    for ( auto &gbc : d->giveBcs() ) {
        ActiveBoundaryCondition *active_bc = dynamic_cast< ActiveBoundaryCondition * >(gbc.get());
        BoundaryCondition *bc = dynamic_cast< BoundaryCondition * >(gbc.get());
        // We only need to consider Dirichlet b.c.s
        if ( bc || ( active_bc && ( active_bc->requiresActiveDofs() || active_bc->giveNumberOfInternalDofManagers() ) ) ) {
            // Check of the dirichlet b.c. has changed in the last step (if so we need to renumber)
            if ( gbc->isImposed(tStep) != gbc->isImposed(tStep->givePreviousStep()) ) {
                return true;
            }
        }
    }
    return false;
}

int
NonStationaryMPMSProblem :: forceEquationNumbering()
{
    this->effectiveMatrix = nullptr;
    return EngngModel :: forceEquationNumbering();
}


void
NonStationaryMPMSProblem :: printOutputAt(FILE *file, TimeStep *tStep)
{
    if ( !this->giveDomain(1)->giveOutputManager()->testTimeStepOutput(tStep) ) {
        return; // Do not print even Solution step header
    }

    EngngModel :: printOutputAt(file, tStep);

    if (this->debug) {
        /* custom export */
        // plot profile of the pressure fields vs x
        Domain *d = this->giveDomain(1);
        double t = tStep->giveTargetTime();
        char filename[256];
        sprintf(filename, "time_%lf.csv", t);
        FILE *file1 = fopen(filename, "w");
        fprintf(file1, "t=%lf[s]\nx, pw, pa, T, u\n", t);
        for ( int i =1; i<= 141; i++) {
            DofManager *dm = d->giveDofManager(i);
            double x = dm->giveCoordinate(1);
            FloatArray pw,pa,u,T,Sw;
            dm->giveUnknownVector(pw, {P_f}, VM_Total, tStep);
            dm->giveUnknownVector(pa, {P_f2}, VM_Total, tStep);
            dm->giveUnknownVector(u, {D_u}, VM_Total, tStep);
            dm->giveUnknownVector(T, {T_f}, VM_Total, tStep);

            // dm->giveUnknownVector(Sw, {S_w}, VM_Total, tStep);
            fprintf(file1, "%le, %le, %le, %le, %le\n", x, pw[0], pa[0], T[0], u[0]);
        }
        fclose(file1);
    }

}


void
NonStationaryMPMSProblem :: updateYourself(TimeStep *tStep)
{
    EngngModel :: updateYourself(tStep);
}

void
NonStationaryMPMSProblem :: saveContext(DataStream &stream, ContextMode mode)
{
    EngngModel :: saveContext(stream, mode);
    field->saveContext(stream);
}


void
NonStationaryMPMSProblem :: restoreContext(DataStream &stream, ContextMode mode)
{
    EngngModel :: restoreContext(stream, mode);
    field->restoreContext(stream);
}

int
NonStationaryMPMSProblem :: checkConsistency()
{
  return EngngModel :: checkConsistency();
}


void
NonStationaryMPMSProblem :: updateDomainLinks()
{
    EngngModel :: updateDomainLinks();
    this->giveNumericalMethod( this->giveCurrentMetaStep() )->setDomain( this->giveDomain(1) );
}


FieldPtr NonStationaryMPMSProblem::giveField(FieldType key, TimeStep *tStep)
{
    /* Note: the current implementation uses MaskedPrimaryField, that is automatically updated with the model progress, 
        so the returned field always refers to active solution step. 
    */

    if ( tStep != this->giveCurrentStep()) {
        OOFEM_ERROR("Unable to return field representation for non-current time step");
    }
    if ( key == FT_Displacements ) {
      return std::make_shared<MaskedPrimaryField>( key, this->field.get(), IntArray{D_u, D_v, D_w} );
    } else if ( key == FT_Pressure ) {
        return std::make_shared<MaskedPrimaryField>( key, this->field.get(), IntArray{P_f} );
    } else if ( key == FT_Temperature ) {
        return std::make_shared<MaskedPrimaryField>( key, this->field.get(), IntArray{T_f} );
    } else {
        return FieldPtr();
    }
}
    
} // end namespace oofem
