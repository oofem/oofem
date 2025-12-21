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


#include "mpm.h"
#include "termlibrary4.h"
#include "element.h"
#include "gausspoint.h"
#include "feinterpol.h"
#include "intarray.h"
#include "connectivitytable.h"
#include "classfactory.h"
#include "gaussintegrationrule.h"
#include "elementinternaldofman.h"
#include "masterdof.h"
#include "engngm.h"

#include "fei2dlinelin.h"
#include "fei2dtrlin.h"
#include "fei3dhexalin.h"
#include "fei3dquadlin.h"

#include "mathfem.h"

#include "material.h"
#include "matstatus.h"

#include "boundaryload.h"
#include "bodyload.h"
#include "zznodalrecoverymodel.h"

#include <vector>

namespace oofem {


/**
 * @brief Base class for Discontinuous Galerkin scalar advection elements 
 * 
 */
class SADGElement : public MPElement {
        
    protected:
        virtual int  giveNumberOfSDofs() const = 0;
        virtual const Variable* getScalarVariable() const = 0;

    protected:
        std :: vector<ElementDofManager*> internalDofManagers;

    public:
    SADGElement(int n, Domain* d): MPElement(n,d) { }

    void getDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {
          answer={num};
    }

    void getInternalDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {
        answer={};
    }
    // Note: performance can be probably improved once it will be possible 
    // to directly assemble multiple term contributions to the system matrix.
    // template metaprogramming?
    void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep) override {
        IntegrationRule* ir = this->giveDefaultIntegrationRulePtr();

        if (type == MassMatrix) { // Mass matrix
            int sdofs = this->giveNumberOfSDofs();
            answer.resize(sdofs,sdofs);
            answer.zero();
            this->integrateTerm_dw (answer, NTN (getScalarVariable(),getScalarVariable()), ir, tStep) ;
        } else if (type == StiffnessMatrix) {
            int sdofs = this->giveNumberOfSDofs();
            answer.resize(sdofs,sdofs);
            answer.zero();
            this->integrateTerm_dw (answer, dnTaN (getScalarVariable(),getScalarVariable(), this->giveDomain()->giveEngngModel()->giveField(FieldType::FT_Velocity, tStep)), ir, tStep) ;
            answer.times(-1.0);
        } else if (type == InternalFluxVector) {
            answer.clear();
        } else {
	        OOFEM_ERROR("Unknown characteristic matrix type");
	    }
    }

    void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) override {
        OOFEM_ERROR("Unknown characteristic vector type");
    }

    void giveCharacteristicMatrixFromBC(FloatMatrix &answer, CharType type, TimeStep *tStep, GeneralBoundaryCondition *bc, int boundaryID) override {
        answer.clear();
    }


    virtual void giveCharacteristicVectorFromBC(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep, GeneralBoundaryCondition *bc, int boundaryID) override {
        answer.clear();
    }

    void computeBoundarySurfaceLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global = true) override {
        OOFEM_ERROR("Unsupported boundary condition type");
    }

    void computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global=true) override {
        OOFEM_ERROR("Unsupported boundary condition type");
    }

    void computeLoadVector(FloatArray &answer, BodyLoad *load, CharType type, ValueModeType mode, TimeStep *tStep) override {
        OOFEM_ERROR("Unsupported load type");
    }


    int computeFluxLBToLRotationMatrix(FloatMatrix &answer, int iSurf, const FloatArray& lc, const Variable::VariableQuantity q, char btype) override {
        answer.clear(); 
        return 0;
    }
};

class SADGBoundaryElement : public SADGElement {
    protected:
        /**
         *  Node numbers of dofs on shared edges/surfaces
         *  They are supposed to be stored for each neighbor boundary in matching order
         *  So for example for 2D triangle element:
         *  - the codeNumbers for shared edge represented by Boundary element the codeNumbers are {2,3}
         *  - and matching codeNumbers for the neighbor element are {5,4}
         * The numbering should be such that the normal vector points from the first element to the second element
         * so the normal vector should point outwards.
         *
         * 
         *   3|\  \4``|6
         *    | \  \  |
         *    | A\  \B|
         *    |   \  \|
         *    1----2   5
         */
    public:
    SADGBoundaryElement(int n, Domain* d): SADGElement(n,d) {}
    /*
        this->neighbors = neighbors;
        this->boundaryIDs = boundaryIDs;

        // attempt to determine matching internal node numbers
        // assuming that i-th internal DofNMan is matching the i-th node  
        for (int i=1; i<=neighbors.giveSize(); i++) {
            Element* neighbor = d->giveElement(neighbors.at(i));
            IntArray bnodes; 
            if (i==1) {
                bnodes = neighbor->giveInterpolation()->boundaryGiveNodes(boundaryIDs.at(i), neighbor->giveGeometryType()); // local numbers
                internalNodeNumbers.push_back(bnodes);
            } else { // need to find matching node numbers on opposite edge/surface 
                IntArray neighborNodes = neighbor->giveDofManArray();
                // loop over nodes of first element and find matching node
                for (int j=1; j<=internalNodeNumbers[0].giveSize(); j++) {
                    int node = d->giveElement(neighbors.at(1))->giveDofManagerNumber(internalNodeNumbers[0].at(j));
                    int indx;
                    if ((indx = neighborNodes.findFirstIndexOf(node))) {
                        bnodes.followedBy(indx);
                        continue;
                    }
                    OOFEM_ERROR("Matching node not found");
                }
                internalNodeNumbers.push_back(bnodes);
            }
        }
    }
    */
    virtual int giveNumberOfSharedElements() const =0;

    void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep) override {
        if (type == InternalFluxVector) {
            
            answer.resize(this->giveNumberOfDofManagers(), this->giveNumberOfDofManagers()); // scalar advection
            answer.zero();

            // code numbers
            IntArray rows, cols;
            // set up integration rule
            IntegrationRule* ir = this->giveDefaultIntegrationRulePtr();
            for (int j=0;j<ir->giveNumberOfIntegrationPoints();j++) {
                FloatArray lc = ir->getIntegrationPoint(j)->giveNaturalCoordinates();
                FloatArray gc;
                this->giveInterpolation()->local2global(gc, lc, FEIElementGeometryWrapper(this));
                FloatArray normal, v, N;
                this->giveInterpolation()->boundaryEvalNormal(normal, 1, lc, FEIElementGeometryWrapper(this));
                this->giveInterpolation()->evalN(N, lc, FEIElementGeometryWrapper(this));
                // get velocity vector first
                // @BP: todo -> connect external field 
                //this->domain->giveEngngModel()->giveField(FT_Velocity, tStep)->evaluateAt(v, gc, ValueModeType::VM_Total, tStep);
                int nsd = this->giveInterpolation()->giveNsd(this->giveGeometryType());
                v.resize(nsd);
                //v.at(2) = 1.0;
                this->giveDomain()->giveEngngModel()->giveField(FT_Velocity, tStep)->evaluateAt(v, gc, ValueModeType::VM_Total, tStep);
                // v.at(1) = sqrt(0.5); v.at(2) = sqrt(0.5); // dummy velocity
                // evaluate N^T (a\cdot n) N
                FloatMatrix contrib;
                contrib.beDyadicProductOf(N,N);
                contrib.times(normal.dotProduct(v));
                contrib.times(ir->getIntegrationPoint(j)->giveWeight()*this->giveInterpolation()->giveTransformationJacobian(lc, FEIElementGeometryWrapper(this)));
                int nnodes = N.giveSize(); // number of nodes on boundary
                rows.resize(nnodes);
                cols.resize(nnodes);
                // determine upwind boundary
                if (normal.dotProduct(v) > 0) {
                    // upwind element is the first neighbor element
                    // set up code numbers for e1 boundary
                    for (int k=1;k<=nnodes;k++) {
                        rows.at(k)=k;
                        cols.at(k)=k;
                    }
                    answer.assemble(contrib, rows, cols);
                    // set up code numbers for e2 boundary, if exists
                    if (this->giveNumberOfSharedElements() > 1) {
                        for (int k=1;k<=nnodes;k++) {
                            rows.at(k)=k+nnodes;
                            cols.at(k)=k;
                        }
                        contrib.times(-1.0); // e2 normal = -e1 normal
                        answer.assemble(contrib, rows, cols);
                    }
                } else {
                    // upwind element is the second neighbor element (if exist)
                    if (this->giveNumberOfSharedElements() > 1) {
                        // e1 contribution
                        for (int k=1;k<=nnodes;k++) {
                            rows.at(k)=k;
                            cols.at(k)=k+nnodes;
                        }
                        answer.assemble(contrib, rows, cols);
                        // e2 contribution
                        for (int k=1;k<=nnodes;k++) {
                            rows.at(k)=k+nnodes;
                            cols.at(k)=k+nnodes;
                        }
                        contrib.times(-1.0); // e2 normal = -e1 normal
                        answer.assemble(contrib, rows, cols);
                    }
                }
            }
        } else if ((type == MassMatrix) || (type == StiffnessMatrix)) {
            answer.clear(); //resize(this->giveNumberOfDofManagers(), this->giveNumberOfDofManagers()); // scalar advection
        } else {
            OOFEM_ERROR("Unknown characteristic matrix type");
        }
    }

     void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) override {
        OOFEM_ERROR("Unknown characteristic vector type");
    }

    void giveCharacteristicMatrixFromBC(FloatMatrix &answer, CharType type, TimeStep *tStep, GeneralBoundaryCondition *bc, int boundaryID) override {
        answer.clear();
    }
};

/**
 * @brief 1D Line linear Element for Discontinuous Galerkin scalar advection
 * 
 */
class SADGBLine1 : public SADGBoundaryElement {
    protected:
        //FEI3dTetLin pInterpol;
        //FEI3dTetQuad uInterpol;
        const static FEI2dLineLin  interpol;
        const static Variable scalarVariable;

      
    public:
    SADGBLine1(int n, Domain* d): 
        SADGBoundaryElement(n,d) { 
            this->numberOfGaussPoints = 2;
    }
    
    void initializeFrom(InputRecord &ir, int priority) override {
        SADGBoundaryElement::initializeFrom(ir, priority);
        this->numberOfDofMans = this->dofManArray.giveSize();
        if (!((numberOfDofMans == 2) || (numberOfDofMans == 4))) {
            OOFEM_ERROR("Invalid number of dofs");
        }
        
        this->computeGaussPoints();
    }


    void giveDofManDofIDMask(int inode, IntArray &answer) const override { 
            int dofid = this->scalarVariable.getDofManDofIDs().at(1);
            answer = {dofid};
    }
    int giveNumberOfDofs() override { return numberOfDofMans; }
    const char *giveInputRecordName() const override {return "sadgbline1";}
    const char *giveClassName() const override { return "SADGBLine1"; }

    
    const FEInterpolation* getGeometryInterpolation() const override {return &this->interpol;}
  
    Element_Geometry_Type giveGeometryType() const override {
        return EGT_line_1;
    }
    int getNumberOfSurfaceDOFs() const override {return 0;}
    int getNumberOfEdgeDOFs() const override {return 0;}
    void getSurfaceLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {
        answer={};
    }
    void getEdgeLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {}
    Interface *giveInterface(InterfaceType it) override {
        return NULL;
    }
    int  giveNumberOfSDofs() const override {return 4;} 
    int giveNumberOfSharedElements() const override {return this->numberOfDofMans/2;}   



private:
        virtual const Variable* getScalarVariable() const override {return &scalarVariable;}
        void computeGaussPoints() override {
            if ( integrationRulesArray.size() == 0 ) {
                integrationRulesArray.resize( 1 );
                integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this);
                integrationRulesArray [ 0 ]->SetUpPointsOnLine(numberOfGaussPoints, _Unknown);
            }
        }
};

const FEI2dLineLin SADGBLine1::interpol = FEI2dLineLin(1,2);
const Variable SADGBLine1::scalarVariable(&SADGBLine1::interpol, Variable::VariableQuantity::VolumeFraction, Variable::VariableType::scalar, 1, NULL, {DofIDItem::C_1});

#define _IFT_SADGBLine1_Name "sadgbline1"
REGISTER_Element(SADGBLine1)


/**
 * @brief 3D Equal order linear Brick TS Element
 * 
 */
class SADGTriangle1 : public SADGElement {
    protected:
        //FEI3dTetLin pInterpol;
        //FEI3dTetQuad uInterpol;
        const static FEI2dTrLin scalarInterpol;
        const static Variable scalarVariable;
      
    public:
    SADGTriangle1(int n, Domain* d): 
        SADGElement(n,d)
    {
        numberOfDofMans  = 3;
        numberOfGaussPoints = 4;
        this->computeGaussPoints();
        /*
        // set up internal dof managers
        // regular nodes just define element geometry
        // internal dof managers are used to store dofs (Discontinuous Galerkin)
        for (int i=0; i<3; i++) {
            ElementDofManager* edm = new ElementDofManager(i, d, this);
            edm->appendDof(new MasterDof(edm, (DofIDItem)scalarVariable.getDofManDofIDs().at(0)));
            internalDofManagers.push_back(edm);
        }
        */
    }

    //int giveNumberOfInternalDofManagers() const override {return 3;}


    void giveDofManDofIDMask(int inode, IntArray &answer) const override { 
            int dofid = this->scalarVariable.getDofManDofIDs().at(1);
            answer = {dofid};
    }
    int giveNumberOfDofs() override { return 3; }
    const char *giveInputRecordName() const override {return "sadgtriangle1";}
    const char *giveClassName() const override { return "SADGTriangle1"; }

    
    const FEInterpolation* getGeometryInterpolation() const override {return &this->scalarInterpol;}
  
    Element_Geometry_Type giveGeometryType() const override {
        return EGT_triangle_1;
    }
    int getNumberOfSurfaceDOFs() const override {return 0;}
    int getNumberOfEdgeDOFs() const override {return 2;}
    void getSurfaceLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {
        answer={};
    }
    void getEdgeLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {}
    Interface *giveInterface(InterfaceType it) override {
        return NULL;
    }

private:
        virtual int  giveNumberOfSDofs() const override {return 3;} 
        virtual const Variable* getScalarVariable() const override {return &scalarVariable;}
        void computeGaussPoints() override {
            if ( integrationRulesArray.size() == 0 ) {
                integrationRulesArray.resize( 1 );
                integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this);
                integrationRulesArray [ 0 ]->SetUpPointsOnTriangle(numberOfGaussPoints, _Unknown);
            }
        }
};

const FEI2dTrLin  SADGTriangle1::scalarInterpol = FEI2dTrLin(1,2);
const Variable SADGTriangle1::scalarVariable(&SADGTriangle1::scalarInterpol, Variable::VariableQuantity::VolumeFraction, Variable::VariableType::scalar, 1, NULL, {DofIDItem::C_1});

#define _IFT_SADGTriangle1_Name "sadgtria1"
REGISTER_Element(SADGTriangle1)

/**
 * @brief 3D SADG linear brick Element
 * 
 */
class SADGBrick1 : public SADGElement {
    protected:
        //FEI3dTetLin pInterpol;
        //FEI3dTetQuad uInterpol;
        const static FEI3dHexaLin scalarInterpol;
        const static Variable scalarVariable;
      
    public:
    SADGBrick1(int n, Domain* d): 
        SADGElement(n,d)
    {
        numberOfDofMans  = 8;
        numberOfGaussPoints = 8;
        this->computeGaussPoints();
        /*
        // set up internal dof managers
        // regular nodes just define element geometry
        // internal dof managers are used to store dofs (Discontinuous Galerkin)
        for (int i=0; i<3; i++) {
            ElementDofManager* edm = new ElementDofManager(i, d, this);
            edm->appendDof(new MasterDof(edm, (DofIDItem)scalarVariable.getDofManDofIDs().at(0)));
            internalDofManagers.push_back(edm);
        }
        */
    }

    //int giveNumberOfInternalDofManagers() const override {return 3;}

    void giveDofManDofIDMask(int inode, IntArray &answer) const override { 
            int dofid = this->scalarVariable.getDofManDofIDs().at(1);
            answer = {dofid};
    }
    int giveNumberOfDofs() override { return 8; }
    const char *giveInputRecordName() const override {return "sadgbrick1";}
    const char *giveClassName() const override { return "SADGBrick1"; }

    
    const FEInterpolation* getGeometryInterpolation() const override {return &this->scalarInterpol;}
  
    Element_Geometry_Type giveGeometryType() const override {
        return EGT_hexa_1;
    }
    int getNumberOfSurfaceDOFs() const override {return 4;}
    int getNumberOfEdgeDOFs() const override {return 2;}
    void getSurfaceLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {
        answer={};
    }
    void getEdgeLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {}
    Interface *giveInterface(InterfaceType it) override {
        return NULL;
    }

private:
        virtual int  giveNumberOfSDofs() const override {return 8;} 
        virtual const Variable* getScalarVariable() const override {return &scalarVariable;}
        void computeGaussPoints() override {
            if ( integrationRulesArray.size() == 0 ) {
                integrationRulesArray.resize( 1 );
                integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this);
                integrationRulesArray [ 0 ]->SetUpPointsOnCube(numberOfGaussPoints, _Unknown);
            }
        }
};

const FEI3dHexaLin  SADGBrick1::scalarInterpol = FEI3dHexaLin();
const Variable SADGBrick1::scalarVariable(&SADGBrick1::scalarInterpol, Variable::VariableQuantity::VolumeFraction, Variable::VariableType::scalar, 1, NULL, {DofIDItem::C_1});

#define _IFT_SADGBrick1_Name "sadgbrick1"
REGISTER_Element(SADGBrick1)

class SADGBQuad1 : public SADGBoundaryElement {
    protected:
        //FEI3dTetLin pInterpol;
        //FEI3dTetQuad uInterpol;
        const static FEI3dQuadLin interpol;
        const static Variable scalarVariable;

      
    public:
    SADGBQuad1(int n, Domain* d): 
        SADGBoundaryElement(n,d) {
            this->numberOfGaussPoints = 8;
        }
    
    void initializeFrom(InputRecord &ir, int priority) override {
        SADGBoundaryElement::initializeFrom(ir, priority);
        this->numberOfDofMans = this->dofManArray.giveSize();
        if (!((numberOfDofMans == 4) || (numberOfDofMans == 8))) {
            OOFEM_ERROR("Invalid number of dofs");
        }
        this->computeGaussPoints();
    }


    void giveDofManDofIDMask(int inode, IntArray &answer) const override { 
            int dofid = this->scalarVariable.getDofManDofIDs().at(1);
            answer = {dofid};
    }
    int giveNumberOfDofs() override { return numberOfDofMans; }
    const char *giveInputRecordName() const override {return "sadgbquad1";}
    const char *giveClassName() const override { return "SADGBQuad1"; }

    
    const FEInterpolation* getGeometryInterpolation() const override {return &this->interpol;}
  
    Element_Geometry_Type giveGeometryType() const override {
        return EGT_quad_1;
    }
    int getNumberOfSurfaceDOFs() const override {return 0;}
    int getNumberOfEdgeDOFs() const override {return 0;}
    void getSurfaceLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {
        answer={};
    }
    void getEdgeLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {}
    Interface *giveInterface(InterfaceType it) override {
        return NULL;
    }
    int  giveNumberOfSDofs() const override {return 8;} 
    int giveNumberOfSharedElements() const override {return this->numberOfDofMans/4;}   



private:
        virtual const Variable* getScalarVariable() const override {return &scalarVariable;}
        void computeGaussPoints() override {
            if ( integrationRulesArray.size() == 0 ) {
                integrationRulesArray.resize( 1 );
                integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this);
                integrationRulesArray [ 0 ]->SetUpPointsOnSquare(numberOfGaussPoints, _Unknown);
            }
        }
};

const FEI3dQuadLin SADGBQuad1::interpol = FEI3dQuadLin();
const Variable SADGBQuad1::scalarVariable(&SADGBQuad1::interpol, Variable::VariableQuantity::VolumeFraction, Variable::VariableType::scalar, 1, NULL, {DofIDItem::C_1});

#define _IFT_SADGBQuad1_Name "sadgbquad1"
REGISTER_Element(SADGBQuad1)



} // end namespace oofem
