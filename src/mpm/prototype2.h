
#ifndef prototype2_h
#define prototype2_h

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
#include "fei2dlinelin.h"
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

    class LinearInterpolation : public FEInterpolation {
    private:
        FEI2dQuadLin  fei2dQuadLin;
        FEI2dLineLin  fei2dLineLin;
    public:
        LinearInterpolation () : FEInterpolation(1), fei2dQuadLin(1,2), fei2dLineLin(1,2) {}

        virtual integrationDomain giveIntegrationDomain(const Element_Geometry_Type egt) const override {
            if (egt == EGT_quad_1) return _Square;
            else if (egt == EGT_line_1) return _Line;
            else return _UnknownIntegrationDomain;
        }
        virtual const Element_Geometry_Type giveGeometryType() const override {
            return EGT_unknown;
        }
        void giveCellDofMans(IntArray& nodes, IntArray& internalDofMans, Element* elem) const override {
            return this->getCellInterpolation(elem->giveGeometryType())->giveCellDofMans(nodes,internalDofMans, elem );
        }
        void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->evalN(answer, lcoords, cellgeo);
        }
        double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->evaldNdx(answer, lcoords, cellgeo);
        }
        void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->local2global(answer, lcoords, cellgeo);
        }
        virtual int global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->global2local(answer, gcoords, cellgeo);
        }
        double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->giveTransformationJacobian(lcoords, cellgeo);
        }
        void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->giveJacobianMatrixAt(jacobianMatrix, lcoords, cellgeo );
        }
        std::unique_ptr<IntegrationRule> giveIntegrationRule(int order, Element_Geometry_Type egt) const override {
            return this->getCellInterpolation(egt)->giveIntegrationRule(order, egt);
        }
        virtual void boundaryEdgeEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->boundaryEdgeEvalN (answer, boundary, lcoords, cellgeo);
        }
        double boundaryEdgeEvalNormal(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->boundaryEdgeEvalNormal(answer, boundary, lcoords, cellgeo);
        }
        double boundaryEdgeGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->boundaryEdgeGiveTransformationJacobian(boundary, lcoords, cellgeo);
        }
        void boundaryEdgeLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->boundaryEdgeLocal2Global(answer, boundary, lcoords, cellgeo);
        }
        integrationDomain giveBoundaryEdgeIntegrationDomain(int boundary, Element_Geometry_Type egt) const override {
            return this->getCellInterpolation(egt)->giveBoundaryEdgeIntegrationDomain(boundary, egt);
        }
        std::unique_ptr<IntegrationRule> giveBoundaryEdgeIntegrationRule(int order, int boundary, Element_Geometry_Type egt) const override {
            return this->getCellInterpolation(egt)->giveBoundaryEdgeIntegrationRule(order, boundary, egt);
        }
        IntArray boundaryEdgeGiveNodes(int boundary, Element_Geometry_Type gt) const override {
            return this->getCellInterpolation(gt)->boundaryEdgeGiveNodes(boundary, gt);
        }
        
        void boundarySurfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->boundarySurfaceEvalN(answer, isurf, lcoords, cellgeo);
        }
        void boundarySurfaceEvaldNdx(FloatMatrix &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->boundarySurfaceEvaldNdx(answer, isurf, lcoords, cellgeo);
        }
        double boundarySurfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->boundarySurfaceEvalNormal(answer, isurf, lcoords, cellgeo);
        }

        void boundarySurfaceLocal2global(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->boundarySurfaceLocal2global(answer, isurf, lcoords, cellgeo);
        }
        double boundarySurfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->boundarySurfaceGiveTransformationJacobian(isurf, lcoords, cellgeo);
        }
        integrationDomain giveBoundarySurfaceIntegrationDomain(int boundary, const Element_Geometry_Type egt) const override {
            return this->getCellInterpolation(egt)->giveBoundarySurfaceIntegrationDomain(boundary, egt);
        }
        std::unique_ptr<IntegrationRule> giveBoundarySurfaceIntegrationRule(int order, int boundary, const Element_Geometry_Type egt) const override {
            return this->getCellInterpolation(egt)->giveBoundarySurfaceIntegrationRule(order, boundary, egt) ;
        }
        IntArray boundarySurfaceGiveNodes(int boundary, const Element_Geometry_Type gt) const override {
            return this->getCellInterpolation(gt)->boundarySurfaceGiveNodes(boundary, gt);
        }
        
        IntArray boundaryGiveNodes(int boundary, const Element_Geometry_Type gt) const override {
            return this->getCellInterpolation(gt)->boundaryGiveNodes(boundary, gt);
        }
        void boundaryEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->boundaryEvalN(answer, boundary, lcoords, cellgeo);
        }
        double boundaryEvalNormal(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->boundaryEvalNormal(answer, boundary, lcoords, cellgeo);
        }
        double boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->boundaryGiveTransformationJacobian(boundary, lcoords, cellgeo);
        }
        void boundaryLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
            return this->getCellInterpolation(cellgeo.giveGeometryType())->boundaryLocal2Global(answer, boundary, lcoords, cellgeo);
        }
        integrationDomain giveBoundaryIntegrationDomain(int boundary, Element_Geometry_Type egt) const override {
            return this->getCellInterpolation(egt)->giveBoundaryIntegrationDomain(boundary,egt);
        }
        std::unique_ptr<IntegrationRule> giveBoundaryIntegrationRule(int order, int boundary, Element_Geometry_Type egt) const override {
            return this->getCellInterpolation(egt)->giveBoundaryIntegrationRule(order, boundary, egt);
        }
        int giveKnotSpanBasisFuncMask(const IntArray &knotSpan, IntArray &mask) const override{ return 0; }
        int giveNumberOfKnotSpanBasisFunctions(const IntArray &knotSpan) const override { return 0; }
        bool hasSubPatchFormulation() const  override { return false; }
        const FloatArray *giveKnotVector() const override { return nullptr; }
        int giveNumberOfKnotSpans(int dim) const override { return 0; }
        const FloatArray *giveKnotValues(int dim) const  override{ return nullptr; }
        const IntArray *giveKnotMultiplicity(int dim) const  override{ return nullptr; }
        int giveNsd(const Element_Geometry_Type egt) const override {
            if (egt == EGT_quad_1) return 2;
            else if (egt == EGT_line_1) return 1;
            else return 0;
        }
        int giveNumberOfEdges(Element_Geometry_Type gt) const override
        { 
            return this->getCellInterpolation(gt)->giveNumberOfEdges(gt);
        }
        int giveNumberOfNodes(Element_Geometry_Type gt) const override
        { 
           return this->getCellInterpolation(gt)-> giveNumberOfNodes(gt);
        }

        const FEInterpolation* getCellInterpolation (Element_Geometry_Type egt) const {
            if (egt==EGT_quad_1) return &fei2dQuadLin;
            else if (egt == EGT_line_1) return &fei2dLineLin;
            else return NULL;
        }
    };


    class MPMSymbolicTerm : public Term {
        public:
        MPMSymbolicTerm (const Variable &testField, const Variable& unknownField, MaterialMode m)  : Term(testField, unknownField, m) {};
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
                //irvec [ size ]->SetUpPointsOnSquare(nip, this->mode);
                irvec[size]->setUpIntegrationPoints(cell.giveIntegrationDomain(), nip, this->mode);
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
    };

    class BTSigmaTerm2 : public MPMSymbolicTerm {
        protected:
        public:
        BTSigmaTerm2 (const Variable &testField, const Variable& unknownField, MaterialMode m)  : MPMSymbolicTerm(testField, unknownField, m) {};

        /**
         * @brief Evaluates the linearization of $B^T\sigma(u)$, i.e. $B^TDBu$
         * 
         * @param answer 
         * @param e 
         * @param coords 
         */
        void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override {
            FloatMatrix D, B, DB;
            e.giveCrossSection()->giveMaterial(gp)->giveCharacteristicMatrix(D, TangentStiffness, gp, tstep);
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
            cell.giveCrossSection()->giveMaterial(gp)->giveCharacteristicVector(sig, eps, Stress, gp, tstep);
            answer.beTProductOf(B, sig);
        }
        void getDimensions(Element& cell) const override {}

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
            int nnodes = interpol.giveNumberOfNodes(cell.giveGeometryType());
            int ndofs = v.size;
            // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
            interpol.evaldNdx(dndx, coords, FEIElementGeometryWrapper(&cell));

            if ((mmode == _3dUP) || (mmode == _3dUPV) || (mode==_3dMat)) {
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
            } else if ((mmode == _PlaneStress)) {
                answer.resize(3, nnodes*ndofs);
                for (int i = 0; i< nnodes; i++) {
                    answer(0, i*ndofs+0) = dndx(i, 0);
                    answer(1, i*ndofs+1) = dndx(i, 1);

                    answer(2, i*ndofs+0) = dndx(i, 1);
                    answer(2, i*ndofs+1) = dndx(i, 0);
                }
            }
        }
        
    };


    class NTfTerm : public MPMSymbolicTerm {
        protected:
        public:
        NTfTerm (const Variable &testField, const Variable& unknownField, MaterialMode m)  : MPMSymbolicTerm(testField, unknownField, m) {};

        void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override {}

        /**
         * @brief Evaluates the residual contribution (rhs)
         * 
         * @param answer 
         * @param e 
         * @param coords 
         */
        void evaluate (FloatArray& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override {
            FloatMatrix N;
            FloatArray nvec, flux={1.,0.};
            const FloatArray& lc = gp->giveNaturalCoordinates();
            
            this->testField.interpolation.evalN(nvec, lc, FEIElementGeometryWrapper(&e));
            N.beNMatrixOf(nvec, testField.size);
            answer.beTProductOf(N, flux);
        }
        
        void getDimensions(Element& cell) const override {}
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
	        BTSigmaTerm2 mt(u,u, _2dUP);
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
            i->assemble_lhs (*effectiveMatrix, EModelDefaultEquationNumbering(), tStep); 
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

    class L1Element : public MPElement  {
        protected:
            const static FEInterpolation & gInterpol;
        public:
            L1Element (int n, Domain* d): MPElement(n,d) {
                numberOfDofMans  = 2; 
            }
            const char *giveInputRecordName() const override {return "L1";}
            const char *giveClassName() const override {return "L1";}
            const FEInterpolation& getGeometryInterpolation() const override {
                return this->gInterpol;
            }
  
            Element_Geometry_Type giveGeometryType() const override {
               return EGT_line_1;
            }
            // MPElement requirements
            void getDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {}
            void getInternalDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {}
            int getNumberOfSurfaceDOFs() const override {return 0;}
            int getNumberOfEdgeDOFs() const override {return 0;}
            void getSurfaceLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {}
            void getEdgeLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {}

    };
    const FEInterpolation & L1Element::gInterpol = FEI2dLineLin(1,2);
    #define _IFT_L1Element_Name "l1"

} // end namespace oofem
#endif // prototype2_h