#include "engngm.h"
#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"
#include "function.h"
#include "dofdistributedprimaryfield.h"
#include "unknownnumberingscheme.h"
#include "integral.h"
#include "fei2dquadlin.h"
#include "termlibrary.h"
#include "crosssection.h"
#include "material.h"
#include "masterdof.h"
#include "gaussintegrationrule.h"

///@name Input fields for testproblem
//@{
#define _IFT_TestProblem_Name "test"
//@}

namespace oofem {
    class MyTerm : public Term {
        protected:
        public:
        MyTerm (const Variable &testField, const Variable& unknownField, MaterialMode m)  : Term(testField, unknownField, m) {};

        /**
         * @brief Evaluates the linearization of $B^T\sigma(u)$, i.e. $B^TDBu$
         * 
         * @param answer 
         * @param e 
         * @param coords 
         */
        void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override {
            FloatMatrix D, B, DB;
            e.giveCrossSection()->giveMaterial(gp)->giveCharacteristicMatrix(D, StiffnessMatrix, gp, tstep);
            this->grad(B, this->field, this->field.interpolation, e, gp->giveNaturalCoordinates(), gp->giveMaterialMode());
            DB.beProductOf(D, B);
            //answer.plusProductSymmUpper(B, DB, 1.0);
            answer.beTProductOf(B,DB);
        }
        /**
         * @brief Evaluates Internal forces vector, i.e. $b^T\sigma(u)$
         * 
         * @param cell 
         * @param coords 
         */
        void evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override {
            FloatArray u, eps, sig;
            FloatMatrix B;
            cell.getUnknownVector(u, this->field, VM_TotalIntrinsic, tstep);
            this->grad(B, this->field, this->field.interpolation, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode());
            eps.beProductOf(B, u);
            cell.giveCrossSection()->giveMaterial(gp)->giveCharacteristicVector(sig, eps, InternalForcesVector, gp, tstep);
            answer.beTProductOf(B, sig);
        }
        void getDimensions(Element& cell) const override {}
        void initializeCell(Element& cell) const override {
            // allocate necessary DOFs
            IntArray enodes, einteranlnodes, dofIDs;
            // process term field
            dofIDs = this->field.getDofManDofIDs();
            this->field.interpolation.giveCellDofMans(enodes, einteranlnodes, &cell);
            for (auto i: enodes) {
                DofManager* dman =  cell.giveDofManager(i);
                for (auto d: dofIDs) {
                    if (!dman->hasDofID((DofIDItem) d)) {
                        // create a DOF
                        MasterDof* dof = new MasterDof(dman, (DofIDItem)d);
                        dman->appendDof(dof);
                    }
                }
            }
            for (auto i: einteranlnodes) {
                DofManager* dman =  cell.giveInternalDofManager(i);
                for (auto d: dofIDs) {
                    if (!dman->hasDofID((DofIDItem) d)) {
                        // create a DOF
                        MasterDof* dof = new MasterDof(dman, (DofIDItem)d);
                        dman->appendDof(dof);
                    }
                }
            }
            // process testField
            dofIDs = this->testField.getDofManDofIDs();
            this->testField.interpolation.giveCellDofMans(enodes, einteranlnodes, &cell);
            for (auto i: enodes) {
                DofManager* dman =  cell.giveDofManager(i);
                for (auto d: dofIDs) {
                    if (!dman->hasDofID((DofIDItem)d)) {
                        // create a DOF
                        MasterDof* dof = new MasterDof(dman, (DofIDItem)d);
                        dman->appendDof(dof);
                    }
                }
            }
            for (auto i: einteranlnodes) {
                DofManager* dman =  cell.giveInternalDofManager(i);
                for (auto d: dofIDs) {
                    if (!dman->hasDofID((DofIDItem)d)) {
                        // create a DOF
                        MasterDof* dof = new MasterDof(dman, (DofIDItem)d);
                        dman->appendDof(dof);
                    }
                }
            }
            // set up the integration rule on cell
            // get required number of IPs
            int myorder = this->field.interpolation.giveInterpolationOrder() *  this->testField.interpolation.giveInterpolationOrder(); 
            GaussIntegrationRule ir(0, &cell);
            int nip = ir.getRequiredNumberOfIntegrationPoints(cell.giveIntegrationDomain(), myorder);
            // create nd insert it toelement if not exist yet.
            std::vector< std :: unique_ptr< IntegrationRule > > &irvec = cell.giveIntegrationRulesArray();
            bool found = false;
            int size = irvec.size();
            for (int i = 0; i< size; i++) {
                if (irvec[i].get()->giveNumberOfIntegrationPoints() == nip) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                // need to insert right one
                irvec.resize( size +1);
                irvec [ size] = std::make_unique<GaussIntegrationRule>(size, &cell);
                irvec [ size ]->SetUpPointsOnSquare(nip, this->mode);
            }
        }
        IntegrationRule* giveElementIntegrationRule(Element* e) const override {
            int myorder = this->field.interpolation.giveInterpolationOrder() *  this->testField.interpolation.giveInterpolationOrder(); 
            GaussIntegrationRule ir(0, e);
            int nip = ir.getRequiredNumberOfIntegrationPoints(e->giveIntegrationDomain(), myorder);
            std::vector< std :: unique_ptr< IntegrationRule > > &irvec = e->giveIntegrationRulesArray();
            int size = irvec.size();
            for (int i = 0; i< size; i++) {
                if (irvec[i].get()->giveNumberOfIntegrationPoints() == nip) {
                    return irvec[i].get();
                }
            }
            return NULL;
        }
        protected:
        /**
         * @brief Evaluates B matrix; i.e. $LN$ where $L$ is operator matrix and $N$ is interpolation matrix of unknowns
         * 
         * @param answer B matrix
         * @param v 
         * @param interpol 
         * @param cell 
         * @param coords 
         */
        void grad(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords, const MaterialMode mmode) const  {
            FloatMatrix dndx;
            int nnodes = interpol.giveNumberOfNodes();
            int ndofs = v.size;
            // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
            interpol.evaldNdx(dndx, coords, FEIElementGeometryWrapper(&cell));

            if ((mmode == _3dUP) || (mmode == _3dUPV)) {
                // 3D mode only now
                answer.resize(6, nnodes*ndofs);
                for (int i = 0; i< nnodes; i++) {
                    answer(0, i*ndofs+0) = dndx(i, 0);
                    answer(1, i*ndofs+1) = dndx(i, 1);
                    answer(2, i*ndofs+2) = dndx(i, 2);

                    answer(3, i*ndofs+1) = dndx(i, 2);
                    answer(3, i*ndofs+2) = dndx(i, 1);

                    answer(4, i*ndofs+0) = dndx(i, 2);
                    answer(4, i*ndofs+2) = dndx(i, 0);

                    answer(5, i*ndofs+0) = dndx(i, 1);
                    answer(5, i*ndofs+1) = dndx(i, 0);
                }   
            } else if ((mmode == _2dUP) || (mmode == _2dUPV)) {
                answer.resize(6, nnodes*ndofs);
                for (int i = 0; i< nnodes; i++) {
                    answer(0, i*ndofs+0) = dndx(i, 0);
                    answer(1, i*ndofs+1) = dndx(i, 1);

                    answer(5, i*ndofs+0) = dndx(i, 1);
                    answer(5, i*ndofs+1) = dndx(i, 0);
                }
            }
        }
        
    };



    class TestProblem : public EngngModel
    {
    protected:
        SparseMtrxType sparseMtrxType = SMT_Skyline;
        std :: unique_ptr< DofDistributedPrimaryField > field;

        std :: unique_ptr< SparseMtrx > effectiveMatrix;

        FloatArray solution;
        FloatArray internalForces;
        FloatArray eNorm;

        /// Numerical method used to solve the problem
        std :: unique_ptr< SparseLinearSystemNM > nMethod;

    public:
        TestProblem(int i, EngngModel * _master) : EngngModel(i, _master) { ndomains = 1;}


        void solveYourselfAt(TimeStep *tStep) override {
            Domain *domain = this->giveDomain(1);
            Set myset (1, domain);
            FEI2dQuadLin interpol(1,2);
            Variable u = Variable(interpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 2, NULL, {1,2});
	        MyTerm mt(u,u, _2dUP);
            myset.setElementList({1});
            this->integralList.push_back(std::make_unique<Integral>(domain, myset, mt));
            Integral *i = this->integralList[0].get();
            i->initialize();
           
            this->forceEquationNumbering();
            OOFEM_LOG_DEBUG("Number of equations %d\n", this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering()) );

            if ( !effectiveMatrix ) {
                effectiveMatrix = classFactory.createSparseMtrx(sparseMtrxType);
                effectiveMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
            }
            i->assemble_dw (*effectiveMatrix, EModelDefaultEquationNumbering(), tStep); 
            effectiveMatrix->printYourself();
        }

        TimeStep *giveNextStep() override
        {
            if ( !currentStep ) {
                // first step -> generate initial step
                //currentStep = std::make_unique<TimeStep>(*giveSolutionStepWhenIcApply());
                currentStep = std::make_unique<TimeStep>(giveNumberOfTimeStepWhenIcApply(), this, 1, 0., 1., 0);
            }
            previousStep = std :: move(currentStep);
            currentStep = std::make_unique<TimeStep>(*previousStep, 1.);

            return currentStep.get();
        }

        // identification
        const char *giveInputRecordName() const { return _IFT_TestProblem_Name; }
        const char *giveClassName() const override { return "TestProblem"; }
        fMode giveFormulation() override { return TL; }
    };
    


    class Q1Element : public MPElement  {
        protected:
            const static FEInterpolation & gInterpol;
        public:
            Q1Element (int n, Domain* d): MPElement(n,d) {
                numberOfDofMans  = 4; 
            }
            const char *giveInputRecordName() const override {return "Q1";}
            const char *giveClassName() const override {return "Q1";}
            const FEInterpolation& getGeometryInterpolation() const override {return this->gInterpol;}
  
            Element_Geometry_Type giveGeometryType() const override {
               return EGT_quad_1;
            }
            // MPElement requirements
            void getDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {}
            void getInternalDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {}
            int getNumberOfSurfaceDOFs() const override {return 0;}
            int getNumberOfEdgeDOFs() const override {return 0;}
            void getSurfaceLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {}
            void getEdgeLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {}

    };
    const FEInterpolation & Q1Element::gInterpol = FEI2dQuadLin(1,2);
    #define _IFT_Q1Element_Name "q1"

} // end namespace oofem
