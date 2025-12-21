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

#ifndef mpmsemodels_h
#define mpmsemodels_h

#include "engngm.h"
#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"
#include "function.h"
#include "dofdistributedprimaryfield.h"
#include "unknownnumberingscheme.h"
#include "integral.h"
#include "nrsolver.h"
#include "connectivitytable.h"

#define _IFT_StationaryMPMSProblem_Name "mpmsymbolicstationaryproblem"


namespace oofem {

    class StationaryMPMSProblem : public EngngModel
    {
    protected:

        SparseMtrxType sparseMtrxType = SMT_Skyline;
        LinSystSolverType solverType = ST_Direct;
        /// This field stores solution vector. For fixed size of problem, the PrimaryField is used, for growing/decreasing size, DofDistributedPrimaryField applies.
        std :: unique_ptr< PrimaryField > unknownsField;

        std :: unique_ptr< SparseMtrx > effectiveMatrix;
        FloatArray solutionVector;
        FloatArray residualVector;
        FloatArray eNorm;

        /// Numerical method used to solve the problem
        std::unique_ptr<SparseNonLinearSystemNM> nMethod;

        bool keepTangent = false;

        // list of integrals contributing to lhs and rhs
        IntArray lhsIntegrals;
        IntArray rhsIntegrals;

    public:
        StationaryMPMSProblem(int i, EngngModel * _master) : EngngModel(i, _master), nMethod(nullptr) { ndomains = 1;}

        void initializeFrom(InputRecord &ir) override {
            EngngModel::initializeFrom(ir);
            IR_GIVE_FIELD(ir, lhsIntegrals, "lhsterms");
            IR_GIVE_FIELD(ir, rhsIntegrals, "rhsterms");
            if ( !unknownsField ) { 
              //unknownsField = std::make_unique<DofDistributedPrimaryField>(this, 1, FT_Unknown, 0);
              unknownsField = std::make_unique<PrimaryField>(this, 1, FT_Unknown, 0);
            }
        }

        void postInitialize() override {
            this->giveDomain(1)->giveConnectivityTable()->buildSharedBoundaryEntities(this->giveDomain(1));
            for (const auto& i: integralList) {
                i->initialize();
            } 
            EngngModel::postInitialize();
        }

        void solveYourselfAt(TimeStep *tStep) override {
            unknownsField->advanceSolution(tStep);

            this->forceEquationNumbering();
            int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );

            OOFEM_LOG_INFO("MPM Symbolic solver\n");
            OOFEM_LOG_INFO("Number of equations %d\n", this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering()) );

            if ( tStep->giveNumber() == 1 ) {
                // allocate space for solution vector
                FloatArray *solutionVector = unknownsField->giveSolutionVector(tStep);
                solutionVector->resize(neq);
                solutionVector->zero();

                effectiveMatrix = classFactory.createSparseMtrx(sparseMtrxType);
                if ( !effectiveMatrix ) {
                    OOFEM_ERROR("LHS sparse matrix creation failed");
                }

                effectiveMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
                if ( this->keepTangent ) {
                    this->effectiveMatrix->zero();
                    // loop over lhs integrals
                    for (auto i: lhsIntegrals) {
                        Integral* integral = this->integralList[i-1].get();
                        integral->assemble_lhs (*effectiveMatrix, EModelDefaultEquationNumbering(), tStep); 
                    }
                }
            }
            // assemble rhs
            FloatArray rhs(this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ));
            // loop over rhs integrals
            for (auto i: rhsIntegrals) {
                Integral* integral = this->integralList[i-1].get();
                integral->assemble_rhs (rhs, EModelDefaultEquationNumbering(), tStep); 
            }
            this->updateSharedDofManagers(rhs, EModelDefaultEquationNumbering(), LoadExchangeTag);

            residualVector.resize(neq);

            // set-up numerical method
            this->giveNumericalMethod( this->giveCurrentMetaStep() );
            FloatArray incrementOfSolution;
            double loadLevel;
            int currentIterations;
            this->nMethod->solve(*this->effectiveMatrix,
                                rhs,
                                NULL,
                                *unknownsField->giveSolutionVector(tStep),
                                incrementOfSolution,
                                residualVector,
                                this->eNorm,
                                loadLevel, // Only relevant for incrementalBCLoadVector
                                SparseNonLinearSystemNM :: rlm_total,
                                currentIterations,
                                tStep);
        }

        void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d) override
        {
            if ( cmpn == InternalRhs ) {
                this->residualVector.zero();
                int maxdofids = d->giveMaxDofID();
                eNorm.resize(maxdofids);
                eNorm.zero();
                for (auto i: lhsIntegrals) {
                    Integral* integral = this->integralList[i-1].get();
                    integral->assemble_rhs (residualVector, EModelDefaultEquationNumbering(), tStep, &eNorm); 
                }
                this->updateSharedDofManagers(this->residualVector, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);
                return;
            } else if ( cmpn == NonLinearLhs ) {
                if ( !this->keepTangent ) {
                    // Optimization for linear problems, we can keep the old matrix (which could save its factorization)
                    this->effectiveMatrix->zero();
                    // loop over lhs integrals
                    for (auto i: lhsIntegrals) {
                        Integral* integral = this->integralList[i-1].get();
                        integral->assemble_lhs (*effectiveMatrix, EModelDefaultEquationNumbering(), tStep); 
                    }
                }
                return;
            } else {
                OOFEM_ERROR("Unknown component");
            }
        }

        void updateSolution(FloatArray &solutionVector, TimeStep *tStep, Domain *d) override
        {}


        TimeStep* giveNextStep() override
        {
            int istep = this->giveNumberOfFirstStep();
            StateCounterType counter = 1;

            if ( currentStep ) {
                istep = currentStep->giveNumber() + 1;
                counter = currentStep->giveSolutionStateCounter() + 1;
            }

            previousStep = std :: move(currentStep);
            currentStep = std::make_unique<TimeStep>(istep, this, 1, ( double ) istep, 0., counter);
            return currentStep.get();
        }

        NumericalMethod *giveNumericalMethod(MetaStep *mStep) override
        {
           if ( !nMethod ) {
                nMethod = std::make_unique<NRSolver>(this->giveDomain(1), this);
            }
            return nMethod.get();
        }

        double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof) override
        // returns unknown quantity like displacement, velocity of equation eq
        // This function translates this request to numerical method language
        {
        #ifdef DEBUG
            int eq = dof->__giveEquationNumber();
            if ( eq == 0 ) {
                OOFEM_ERROR("invalid equation number");
            }
        #endif

            if (mode == VM_TotalIntrinsic) mode = VM_Total;
            return unknownsField->giveUnknownValue(dof, mode, tStep);
        }
        // identification
        const char *giveInputRecordName() const { return _IFT_StationaryMPMSProblem_Name; }
        const char *giveClassName() const override { return "StationaryMPMSProblem"; }
        fMode giveFormulation() override { return TL; }
    };
        
} // end namespace oofem
#endif // mpmsemodels_h
