
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
#include "fei2dquadconst.h"
#include "fei2dlineconst.h"
#include "fei2dquadlin.h"
#include "fei2dlinelin.h"
#include "fei2dquadquad.h"
#include "fei2dlinequad.h"
#include "termlibrary.h"
#include "crosssection.h"
#include "material.h"
#include "masterdof.h"
#include "gaussintegrationrule.h"
#include "dofmanager.h"
#include "connectivitytable.h"
#include "matresponsemode.h"

///@name Input fields for testproblem
//@{
#define _IFT_TestProblem_Name "test"
//@}

namespace oofem {

    /**
     * @brief CellTypeUnifiedInterpolation represent unification of interpolations of the same type (linear, quadratic) 
     * defined for specific elements set of cell types to work on set of different cell types.
     * 
     */
    class CellTypeUnifiedInterpolation : public FEInterpolation {
    public:
        CellTypeUnifiedInterpolation (int o) : FEInterpolation(o) {}
        virtual const FEInterpolation* getCellInterpolation (Element_Geometry_Type egt) const = 0;

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
        std::unique_ptr<IntegrationRule> giveIntegrationRule(int _order, Element_Geometry_Type egt) const override {
            return this->getCellInterpolation(egt)->giveIntegrationRule(_order, egt);
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
        std::unique_ptr<IntegrationRule> giveBoundaryEdgeIntegrationRule(int _order, int boundary, Element_Geometry_Type egt) const override {
            return this->getCellInterpolation(egt)->giveBoundaryEdgeIntegrationRule(_order, boundary, egt);
        }
        IntArray boundaryEdgeGiveNodes(int boundary, Element_Geometry_Type gt, bool includeHierarchical=false) const override {
            return this->getCellInterpolation(gt)->boundaryEdgeGiveNodes(boundary, gt, includeHierarchical);
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
        std::unique_ptr<IntegrationRule> giveBoundarySurfaceIntegrationRule(int _order, int boundary, const Element_Geometry_Type egt) const override {
            return this->getCellInterpolation(egt)->giveBoundarySurfaceIntegrationRule(_order, boundary, egt) ;
        }
        IntArray boundarySurfaceGiveNodes(int boundary, const Element_Geometry_Type gt, bool includeHierarchical=false) const override {
            return this->getCellInterpolation(gt)->boundarySurfaceGiveNodes(boundary, gt, includeHierarchical);
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
        std::unique_ptr<IntegrationRule> giveBoundaryIntegrationRule(int _order, int boundary, Element_Geometry_Type egt) const override {
            return this->getCellInterpolation(egt)->giveBoundaryIntegrationRule(_order, boundary, egt);
        }
        int giveKnotSpanBasisFuncMask(const IntArray &knotSpan, IntArray &mask) const override{ return 0; }
        int giveNumberOfKnotSpanBasisFunctions(const IntArray &knotSpan) const override { return 0; }
        bool hasSubPatchFormulation() const  override { return false; }
        const FloatArray *giveKnotVector() const override { return nullptr; }
        int giveNumberOfKnotSpans(int dim) const override { return 0; }
        const FloatArray *giveKnotValues(int dim) const  override{ return nullptr; }
        const IntArray *giveKnotMultiplicity(int dim) const  override{ return nullptr; }
        
        int giveNumberOfEdges(Element_Geometry_Type gt) const override
        { 
            return this->getCellInterpolation(gt)->giveNumberOfEdges(gt);
        }
        int giveNumberOfNodes(Element_Geometry_Type gt) const override
        { 
           return this->getCellInterpolation(gt)-> giveNumberOfNodes(gt);
        } 
        /**
         * Return dofmanager numbering offsets for given cell type and interpolation order
         * 
         * @param egt cell geometry type
         * @param order interpolation order
         * @param dofManOffset dofman numbering offset for new dofmanagers to be introduced for given interpolation order
         * @param internalDofManOffset internal dofman numbering offset for new dofmanagers to be introduced for given interpolation order
         */
        virtual void giveCellDofManNumberingOffsets(Element_Geometry_Type egt, int order, int &dofManOffset, int&internalDofManOffset) const {
            if (egt == EGT_line_1) {
                if ((order == 0) || (order == 1)) { // constant
                    dofManOffset = 0;
                    internalDofManOffset = 0;
                    return;
                } else if (order == 2) {
                    dofManOffset = 2;
                    internalDofManOffset = 0;
                    return;
                }
            } else if (egt==EGT_quad_1) {
                if ((order == 0) || (order == 1)) { // constant
                    dofManOffset = 0;
                    internalDofManOffset = 0;
                    return;
                } else if (order == 2) {
                    dofManOffset = 4;
                    internalDofManOffset = 0;
                    return;
                }
            }
            OOFEM_ERROR("Unsupported element geometry type (%d) and interpolation order (%d)", egt, order);
        }

    };

    class ConstantInterpolation : public CellTypeUnifiedInterpolation {
    private:
        FEI2dQuadConst  fei2dQuadConst;
        FEI2dLineConst  fei2dLineConst;   
    public:
        ConstantInterpolation () : CellTypeUnifiedInterpolation(0), fei2dQuadConst(1,2), fei2dLineConst(1,2) {}

        virtual integrationDomain giveIntegrationDomain(const Element_Geometry_Type egt) const override {
            if (egt == EGT_quad_1) return _Square;
            else if (egt == EGT_line_1) return _Line;
            else return _UnknownIntegrationDomain;
        }
        virtual const Element_Geometry_Type giveGeometryType() const override {
            return EGT_unknown;
        }
        virtual const Element_Geometry_Type giveBoundaryGeometryType(int boundary) const override {
            return EGT_unknown;
        }
        
        int giveNsd(const Element_Geometry_Type egt) const override {
            if (egt == EGT_quad_1) return 2;
            else if (egt == EGT_line_1) return 1;
            else return 0;
        }
        
        const FEInterpolation* getCellInterpolation (Element_Geometry_Type egt) const override {
            if (egt==EGT_quad_1) return &fei2dQuadConst;
            else if (egt == EGT_line_1) return &fei2dLineConst;
            else return NULL;
        }
        void initializeCell(Element* e) const override {
            // initialize cell to be compatible with interpolation
            // typically 1) new nodes are allocated on shared edges and surfaces with neighboring elements 
            // and (2)new internal dofmans are allocated
            Domain *d = e->giveDomain();
            Element_Geometry_Type egt=e->giveGeometryType();
            if ((egt == EGT_quad_2) || (egt == EGT_quad_1)) { 
                if (e->giveNumberOfInternalDofManagers() < 1) {
                    std::unique_ptr<DofManager> dm = std::make_unique<DofManager>(1, d);
                    e->setInternalDofManager(1, std::move(dm));
                }
            } else {
                OOFEM_ERROR ("Unsupported element geometry type (%d)", egt);
            }
        }
    };

    class LinearInterpolation : public CellTypeUnifiedInterpolation {
    private:
        FEI2dQuadLin  fei2dQuadLin;
        FEI2dLineLin  fei2dLineLin;
    public:
        LinearInterpolation () : CellTypeUnifiedInterpolation(1), fei2dQuadLin(1,2), fei2dLineLin(1,2) {}

        virtual integrationDomain giveIntegrationDomain(const Element_Geometry_Type egt) const override {
            if (egt == EGT_quad_1) return _Square;
            else if (egt == EGT_line_1) return _Line;
            else return _UnknownIntegrationDomain;
        }
        virtual const Element_Geometry_Type giveGeometryType() const override {
            return EGT_unknown;
        }
        virtual const Element_Geometry_Type giveBoundaryGeometryType(int boundary) const override {
            return EGT_unknown;
        }
        
        int giveNsd(const Element_Geometry_Type egt) const override {
            if (egt == EGT_quad_1) return 2;
            else if (egt == EGT_line_1) return 1;
            else return 0;
        }
        
        const FEInterpolation* getCellInterpolation (Element_Geometry_Type egt) const override {
            if (egt==EGT_quad_1) return &fei2dQuadLin;
            else if (egt == EGT_line_1) return &fei2dLineLin;
            else return NULL;
        }
        void initializeCell(Element* e) const override {}
    };

    class QuadraticInterpolation : public CellTypeUnifiedInterpolation {
    private:
        FEI2dQuadQuad  fei2dQuadQuad;
        FEI2dLineQuad  fei2dLineQuad;
    public:
        QuadraticInterpolation () : CellTypeUnifiedInterpolation(2), fei2dQuadQuad(1,2), fei2dLineQuad(1,2) {}

        virtual integrationDomain giveIntegrationDomain(const Element_Geometry_Type egt) const override {
            if ((egt == EGT_quad_1) || (egt==EGT_quad_2)) return _Square;
            else if ((egt == EGT_line_1) || (egt == EGT_line_2)) return _Line;
            else return _UnknownIntegrationDomain;
        }
        virtual const Element_Geometry_Type giveGeometryType() const override {
            return EGT_unknown;
        }
        virtual const Element_Geometry_Type giveBoundaryGeometryType(int boundary) const override {
            return EGT_unknown;
        }
        
        int giveNsd(const Element_Geometry_Type egt) const override {
            if ((egt == EGT_quad_1)||(egt==EGT_quad_2)) return 2;
            else if ((egt == EGT_line_1)||(egt==EGT_line_2)) return 1;
            else return 0;
        }
        
        const FEInterpolation* getCellInterpolation (Element_Geometry_Type egt) const override {
            if ((egt==EGT_quad_1)||(egt==EGT_quad_2)) return &fei2dQuadQuad;
            else if ((egt == EGT_line_1)||(egt==EGT_line_2)) return &fei2dLineQuad;
            else return NULL;
        }
        void initializeCell(Element* e) const override {
            // initialize cell to be compatible with interpolation
            // typically 1) new nodes are allocated on shared edges and surfaces with neighboring elements 
            // and (2)new internal dofmans are allocated
            Domain *d = e->giveDomain();
            Element_Geometry_Type egt=e->giveGeometryType();
            if ((egt == EGT_quad_2) || (egt == EGT_line_2)) {
                return;
            } else if (egt == EGT_line_1) {
                this->allocateDofMans(e); 
                if ((e->giveNumberOfDofManagers() > 2) && (e->giveDofManagerNumber(3) ==0)) {
                    int ndm = d->giveNumberOfDofManagers()+1;
                    std::unique_ptr<DofManager> dm = std::make_unique<DofManager>(ndm, d);
                    d->resizeDofManagers(ndm);
                    d->setDofManager(ndm, std::move(dm));
                    e->setDofManager(3, ndm);
                }
            } else if (egt == EGT_quad_1) {
                this->allocateDofMans(e);
                // loop over edges shared with other element5s
                const IntArray* edges = e->giveSharedEdgeIDs();
                for (int i = 1; i<= 4; i++) {
                    if (e->giveDofManagerNumber(4+i) == 0) {
                        int ndm = d->giveNumberOfDofManagers()+1;
                        std::unique_ptr<DofManager> dm= std::make_unique<DofManager>(ndm, d);
                        d->resizeDofManagers(ndm);
                        d->setDofManager(ndm, std::move(dm));
                        e->setDofManager(4+i, ndm);
                    
                        // loop over neighbors
                        SharedBoundaryEntity* sbe = d->giveConnectivityTable()->giveBoundaryEntity(edges->at(i));
                        for (auto ier: sbe->elements) {
                            // find corresponding edge in neighbor element and set DofManager
                            int dofManOffset, internalDofManOffset;
                            giveCellDofManNumberingOffsets(d->giveElement(ier.elementID)->giveGeometryType(), this->giveInterpolationOrder(), dofManOffset, internalDofManOffset);
                            this->allocateDofMans(d->giveElement(ier.elementID));
                            d->giveElement(ier.elementID)->setDofManager(dofManOffset+ier.boundaryID, ndm);
                        }
                    }
                }
            } else {
                OOFEM_ERROR ("Unsupported element geometry type (%d)", egt);
            }
        }
        protected:
        void allocateDofMans(Element* e) const {
            Element_Geometry_Type egt=e->giveGeometryType();
            if (egt == EGT_line_1) {
                if (e->giveNumberOfDofManagers() < 3) {
                    IntArray enodes = {e->giveDofManagerNumber(1), e->giveDofManagerNumber(2), 0};
                    e->setNumberOfDofManagers(3);
                    e->setDofManagers(enodes);
                }
            } if (egt == EGT_quad_1) {
                const IntArray &nodes = e->giveDofManArray();
                int nnodes = nodes.giveSize();
                if (nnodes < 8) {
                    // basic cell with linear geometry 

                    IntArray enodes(8);
                    for (int i=1; i<=nnodes; i++) {
                        enodes.at(i) = nodes.at(i);
                    }
                    e->setNumberOfDofManagers(8);
                    e->setDofManagers(enodes);
                }
            }
        }
    };


    class BTSigmaTerm2 : public MPMSymbolicTerm {
        protected:
            MatResponseMode lhsmatmode = MatResponseMode::TangentStiffness;
            MatResponseMode rhsmatmode = MatResponseMode::Stress;

        public:
        BTSigmaTerm2() : MPMSymbolicTerm() {}
        BTSigmaTerm2 (const Variable *testField, const Variable* unknownField, MaterialMode m)  : MPMSymbolicTerm(testField, unknownField, m) {};

        /**
         * @brief Evaluates the linearization of $B^T\sigma(u)$, i.e. $B^TDBu$
         * 
         * @param answer 
         * @param e 
         * @param coords 
         */
        void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override {
            FloatMatrix D, B, DB;
            e.giveCrossSection()->giveMaterial(gp)->giveCharacteristicMatrix(D, this->lhsmatmode, gp, tstep);
            this->grad(B, this->field, this->field->interpolation, e, gp->giveNaturalCoordinates(), gp->giveMaterialMode());
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
            this->grad(B, this->field, this->field->interpolation, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode());
            eps.beProductOf(B, u);
            cell.giveCrossSection()->giveMaterial(gp)->giveCharacteristicVector(sig, eps, this->rhsmatmode, gp, tstep);
            answer.beTProductOf(B, sig);
        }
        void getDimensions(Element& cell) const override {}
        void initializeFrom(InputRecord &ir, EngngModel* problem) override {
            MPMSymbolicTerm::initializeFrom(ir, problem);
            int value;
            if (ir.hasField("lhsmatmode")) {
                IR_GIVE_FIELD(ir, value, "lhsmatmode");
                lhsmatmode = (MatResponseMode) value;
            }
            if (ir.hasField("rhsmatmode")) {
                IR_GIVE_FIELD(ir, value, "rhsmatmode");
                rhsmatmode = (MatResponseMode) value;
            }
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
        void grad(FloatMatrix& answer, const Variable *v, const FEInterpolation* interpol, const MPElement& cell, const FloatArray& coords, const MaterialMode mmode) const  {
            FloatMatrix dn, dndx, jacobianMatrix, inv;
            int nnodes = interpol->giveNumberOfNodes(cell.giveGeometryType());
            int ndofs = v->size;
            // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj

            interpol->evaldNdx(dndx, coords, FEIElementGeometryWrapper(&cell));   // won't work for cells with lower geometry intyerpolation than unknown interpolation
            // use cell geometry interpolation to evaluate jacobian 
            // @TODO would be better if this is handled by interpolation class internally by passing geometry interpolation somhow (part of FEIElementGeometryWrapper?)
            //interpol->evaldNdxi(dn, coords, FEIElementGeometryWrapper(&cell));
            //const FEInterpolation *gi = cell.getGeometryInterpolation();
            //gi->giveJacobianMatrixAt(jacobianMatrix, coords, FEIElementGeometryWrapper(&cell));
            //inv.beInverseOf(jacobianMatrix);
            //dndx.beProductTOf(dn, inv);


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
            } else if (mmode == _PlaneStrain) {
                answer.resize(4, nnodes*ndofs);
                for (int i = 0; i< nnodes; i++) {
                    answer(0, i*ndofs+0) = dndx(i, 0);
                    answer(1, i*ndofs+1) = dndx(i, 1);

                    answer(3, i*ndofs+0) = dndx(i, 1);
                    answer(3, i*ndofs+1) = dndx(i, 0);
                }
            } else {
                OOFEM_ERROR("Unsupported material mode %d", mmode);
            }
        }
        
    };
    #define _IFT_BTSigmaTerm2_Name "BTSigmaTerm"

    class NTfTerm : public MPMSymbolicTerm {
        protected:
        FloatArray flux;
        public:
        NTfTerm () : MPMSymbolicTerm() {}
        NTfTerm (const Variable *testField, const Variable* unknownField, MaterialMode m)  : MPMSymbolicTerm(testField, unknownField, m) {};
        NTfTerm (const Variable *testField, const Variable* unknownField, MaterialMode m, const FloatArray* f)  : MPMSymbolicTerm(testField, unknownField, m) {flux = *f;};

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
            FloatArray nvec;
            const FloatArray& lc = gp->giveNaturalCoordinates();
            
            this->testField->interpolation->evalN(nvec, lc, FEIElementGeometryWrapper(&e));
            N.beNMatrixOf(nvec, testField->size);
            answer.beTProductOf(N, flux);
        }
        
        void getDimensions(Element& cell) const override {}
        void initializeFrom(InputRecord &ir, EngngModel* problem) override {
            MPMSymbolicTerm::initializeFrom(ir, problem);
            IR_GIVE_FIELD(ir, flux, "flux");
        }
    };

    #define _IFT_NTfTerm_Name "NTfTerm"

    class TestProblem : public EngngModel
    {
    protected:
        SparseMtrxType sparseMtrxType = SMT_Skyline;
        LinSystSolverType solverType = ST_Direct;
        std :: unique_ptr< DofDistributedPrimaryField > field;

        std :: unique_ptr< SparseMtrx > effectiveMatrix;

        FloatArray solution;
        FloatArray internalForces;
        FloatArray eNorm;

        /// Numerical method used to solve the problem
        std :: unique_ptr< SparseLinearSystemNM > nMethod;

        // list of integrals contributing to lhs and rhs
        IntArray lhsIntegrals;
        IntArray rhsIntegrals;

    public:
        TestProblem(int i, EngngModel * _master) : EngngModel(i, _master) { ndomains = 1;}

        void initializeFrom(InputRecord &ir) override {
            EngngModel::initializeFrom(ir);
            IR_GIVE_FIELD(ir, lhsIntegrals, "lhsterms");
            IR_GIVE_FIELD(ir, rhsIntegrals, "rhsterms");
        }

        void postInitialize() override {
            
            this->giveDomain(1)->giveConnectivityTable()->buildSharedBoundaryEntities(this->giveDomain(1));
            for (const auto& i: integralList) {
                i->initialize();
            } 
            EngngModel::postInitialize();
            // initialize integrals
            //this->giveDomain(1)->postInitialize();
        }

        void solveYourselfAt(TimeStep *tStep) override {
            FloatArray rhs;
            /*
            Domain *domain = this->giveDomain(1);
            Set myset (1, domain);
            FEI2dQuadLin interpol(1,2);
            Variable u = Variable(&interpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 2, NULL, {1,2});
	        BTSigmaTerm2 mt(&u,&u, _2dUP);
            myset.setElementList({1});
            this->integralList.push_back(std::make_unique<Integral>(domain, &myset, &mt));
            Integral *i = this->integralList[0].get();
            i->initialize();
            */

            this->forceEquationNumbering();
            OOFEM_LOG_INFO("MPM Symbolic solver\n");
            OOFEM_LOG_INFO("Number of equations %d\n", this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering()) );

#if 0
            // print element connectivity and nodes after interpolation initialization
            for (auto &n: this->giveDomain(1)->dofManagerList) {
                printf("Node %3d, code numbers: %4d %4d\n", n->giveNumber(), n->giveDofWithID(D_u)->__giveEquationNumber(), n->giveDofWithID(D_v)->__giveEquationNumber());
            }

            for (auto& e: this->giveDomain(1)->elementList) {
                printf("Element %4d, nodes ", e->giveNumber());
                IntArray nodes = e->giveDofManArray();
                for (auto n: nodes) {
                    printf("%4d ",n);
                }
                printf("\n");
            }
#endif

            if ( !effectiveMatrix ) {
                effectiveMatrix = classFactory.createSparseMtrx(sparseMtrxType);
                effectiveMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
            }
            // loop over lhs integrals
            for (auto i: lhsIntegrals) {
                Integral* integral = this->integralList[i-1].get();
                integral->assemble_lhs (*effectiveMatrix, EModelDefaultEquationNumbering(), tStep); 
            }
            //effectiveMatrix->printYourself();
            // assemble rhs
            rhs.resize(this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ));
            // loop over rhs integrals
            for (auto i: rhsIntegrals) {
                Integral* integral = this->integralList[i-1].get();
                integral->assemble_rhs (rhs, EModelDefaultEquationNumbering(), tStep); 
            }
            //rhs.printYourself();
            solution.resize(rhs.giveSize());
            // solve the system
            nMethod->solve(*effectiveMatrix, rhs, solution);
            //solution.printYourself();
            
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
        NumericalMethod *giveNumericalMethod(MetaStep *mStep) override
        {
            if ( !nMethod ) {
                if ( isParallel() ) {
                    if ( ( solverType == ST_Petsc ) || ( solverType == ST_Feti ) ) {
                        nMethod = classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this);
                    }
                } else {
                    nMethod = classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this);
                }
                if ( !nMethod ) {
                    OOFEM_ERROR("linear solver creation failed for lstype %d", solverType);
                }
            }
            return nMethod.get();
        }

        double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof) override
        // returns unknown quantity like displacement, velocity of equation eq
        // This function translates this request to numerical method language
        {
            int eq = dof->__giveEquationNumber();
        #ifdef DEBUG
            if ( eq == 0 ) {
                OOFEM_ERROR("invalid equation number");
            }
        #endif

            if ( tStep != this->giveCurrentStep() ) {
                OOFEM_ERROR("unknown time step encountered");
            }

            switch ( mode ) {
            case VM_Total:
            case VM_Incremental:
                if ( solution.isNotEmpty() ) {
                    return solution.at(eq);
                } else {
                    return 0.;
                }

            default:
                OOFEM_ERROR("Unknown is of undefined type for this problem");
            }

            return 0.;
        }
        // identification
        const char *giveInputRecordName() const { return _IFT_TestProblem_Name; }
        const char *giveClassName() const override { return "TestProblem"; }
        fMode giveFormulation() override { return TL; }
    };
    


    class Q1Element : public MPElement  {
        protected:
            const static FEInterpolation & gInterpol;
            std::vector<std::shared_ptr<DofManager>> internalDofManagers; // stored in the domain not internally
        public:
            Q1Element (int n, Domain* d): MPElement(n,d), internalDofManagers() {
                numberOfDofMans  = 4; 
            }
            const char *giveInputRecordName() const override {return "Q1";}
            const char *giveClassName() const override {return "Q1";}
            const FEInterpolation* getGeometryInterpolation() const override {return &this->gInterpol;}
  
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
            DofManager *giveInternalDofManager(int i) const override {
                    return this->internalDofManagers.at(i-1).get();
            }
            void setInternalDofManager(int num, std::unique_ptr<DofManager> dm) override{
                if (num > (int)this->internalDofManagers.size()) {
                    this->internalDofManagers.resize(num);
                }
                this->internalDofManagers[num-1] = std::move(dm);
            }
            int giveNumberOfInternalDofManagers() const override  { return this->internalDofManagers.size(); }

            IntArray giveBoundaryEdgeNodes(int boundary, bool includeHierarchical=false) const override {
                IntArray answer = MPElement::giveBoundaryEdgeNodes(boundary, false);
                int nnodes = this->dofManArray.giveSize();
                if (includeHierarchical && nnodes > 4) {  
                    if (nnodes <= 8) { // qudratic
                        int nbaseenodes = answer.giveSize();
                        int offset = 3+boundary;
                        int nhiernodes = 1;
                        answer.resizeWithValues(nbaseenodes+nhiernodes);
                        for (int i=1; i<=nhiernodes; i++) {
                            answer.at(nbaseenodes+i) = i+offset;
                        } 
                    } else {
                        OOFEM_ERROR("Unsupported hierarchical node scheme (nnodes=%d)", nnodes);
                    }
                }
                return answer;
            }
    };
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
            const FEInterpolation* getGeometryInterpolation() const override {
                return &this->gInterpol;
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
    #define _IFT_L1Element_Name "l1"

} // end namespace oofem
#endif // prototype2_h
