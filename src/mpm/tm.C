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
#include "termlibrary3.h"
#include "termlibrary.h"
#include "element.h"
#include "gausspoint.h"
#include "feinterpol.h"
#include "intarray.h"
#include "classfactory.h"
#include "gaussintegrationrule.h"

#include "fei3dtetlin.h"
#include "fei3dtetquad.h"
#include "fei3dhexalin.h"
#include "fei2dquadlin.h"
#include "mathfem.h"

#include "material.h"
#include "matstatus.h"

#include "boundaryload.h"
#include "bodyload.h"
#include "zznodalrecoverymodel.h"

namespace oofem {


/**
 * @brief Base class for fully coupled, nonlinear thermo mechanical elements
 * 
 */
class TMElement : public MPElement {
        
    private:
        virtual int  giveNumberOfUDofs() const = 0;
        virtual int  giveNumberOfTDofs() const = 0;
        virtual const Variable* getU() const = 0;
        virtual const Variable* getT() const = 0;

    public:
    TMElement(int n, Domain* d): 
        MPElement(n,d) { }

    // Note: performance can be probably improved once it will be possible 
    // to directly assemble multiple term contributions to the system matrix.
    // template metaprogramming?
    void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep) override {
        IntegrationRule* ir = this->giveDefaultIntegrationRulePtr();

        if (type == MomentumBalance_StiffnessMatrix) {
            int udofs = this->giveNumberOfUDofs();
            answer.resize(udofs,udofs);
            answer.zero();
            this->integrateTerm_dw (answer, TMBTSigTerm (getU(),getU(), getT()), ir, tStep) ;
        } else if (type == MomentumBalance_ThermalCouplingMatrix) {
            int udofs = this->giveNumberOfUDofs();
            int tdofs = this->giveNumberOfTDofs(); 
            answer.resize(udofs,tdofs);
            answer.zero();
            this->integrateTerm_dw (answer, BDalphaPiTerm (getU(),getT(), VM_TotalIntrinsic), ir, tStep) ;
            this->integrateTerm_dw (answer, BTdSigmadT(getU(),getT()), ir, tStep) ;
        } else if (type == EnergyBalance_ConductivityMatrix) {
            int tdofs = this->giveNumberOfTDofs();
            answer.resize(tdofs,tdofs);
            answer.zero();
            this->integrateTerm_dw (answer, TMgNTfTerm(getT(),getT(), Conductivity, Flux), ir, tStep) ;
        } else if (type == EnergyBalance_CapacityMatrix) {
            int tdofs = this->giveNumberOfTDofs();
            answer.resize(tdofs,tdofs);
            answer.zero();
            this->integrateTerm_dw (answer, NTcN(getT(), getT(), Capacity), ir, tStep) ;
        } else {
	        OOFEM_ERROR("Unknown characteristic matrix type");
	    }
    }

    void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) override {
        IntegrationRule* ir = this->giveDefaultIntegrationRulePtr();
        if (type == MomentumBalance_StressResidual) {
            answer.resize(this->giveNumberOfUDofs());
            answer.zero();
            this->integrateTerm_c (answer, TMBTSigTerm(getU(),getU(),getT()), ir, tStep) ;
        } else if (type == EnergyBalance_Residual) {
            answer.resize(this->giveNumberOfTDofs());
            answer.zero();
            this->integrateTerm_c(answer, TMgNTfTerm(getT(),getT(), Conductivity, Flux), ir, tStep) ;
            answer.times(-1.0);
            this->integrateTerm_c (answer, NTcN(getT(), getT(), Capacity), ir, tStep) ;
            // add internal (material generated) heat source r(T) to the residual
            FloatArray temp;
            this->integrateTerm_c(temp, InternalTMFluxSourceTerm(getT(),getU(),getT()), ir, tStep);
            answer.subtract(temp);
        } else if (type == ExternalForcesVector) {
          answer.zero();
        } else {
	        OOFEM_ERROR("Unknown characteristic vector type");
	    }
    }

    void giveCharacteristicMatrixFromBC(FloatMatrix &answer, CharType type, TimeStep *tStep, GeneralBoundaryCondition *bc, int boundaryID) override {
        if (bc->giveType() == ConvectionBC) {
            BoundaryLoad *bbc = dynamic_cast<BoundaryLoad*>(bc);
            if (bbc) {
                FloatMatrix contrib;
                IntArray loc;
                int iorder = getU()->interpolation->giveInterpolationOrder()+bbc->giveApproxOrder();
                std::unique_ptr<IntegrationRule> ir;
                answer.clear();
                if (bbc->giveBCGeoType() == bcGeomType::SurfaceLoadBGT) {
                        ir = this->getGeometryInterpolation()->giveBoundarySurfaceIntegrationRule(iorder, 1, this->giveGeometryType());
                        //this->integrateSurfaceTerm_dw(contrib, NTf_Surface(getT(), BoundaryFluxFunctor(bbc, boundaryID, getT().dofIDs,'s'), boundaryID), ir.get(), boundaryID, tStep);
                        this->integrateSurfaceTerm_dw(answer, NTaTmTe(getT(), getT(), bbc, boundaryID, 's'), ir.get(), boundaryID, tStep);
                } else if (bbc->giveBCGeoType() == bcGeomType::EdgeLoadBGT) {
                        ir = this->getGeometryInterpolation()->giveBoundaryEdgeIntegrationRule(iorder, 1, this->giveGeometryType());
                        //this->integrateSurfaceTerm_dw(contrib, NTf_Surface(getT(), BoundaryFluxFunctor(bbc, boundaryID, getT().dofIDs,'e'), boundaryID), ir.get(), boundaryID, tStep);
                        this->integrateSurfaceTerm_dw(answer, NTaTmTe(getT(), getT(), bbc, boundaryID, 'e'), ir.get(), boundaryID, tStep);
                } else {
                    OOFEM_ERROR("Unsupported boundary condition geometry type");
                }
            }
        } else {
            answer.clear();
        }
    }


    virtual void giveCharacteristicVectorFromBC(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep, GeneralBoundaryCondition *bc, int boundaryID) override {
        
        if ((type == EnergyBalance_ConvectionBCResidual) && (bc->giveType() == ConvectionBC)) {
            BoundaryLoad *bl = dynamic_cast<BoundaryLoad*>(bc);
            IntArray loct, tc;
            FloatArray contrib2;
            answer.clear();            
            int o = getT()->interpolation->giveInterpolationOrder()+bl->giveApproxOrder();
            if (bc->giveBCGeoType() == bcGeomType::SurfaceLoadBGT) {
                std::unique_ptr<IntegrationRule> ir2 = this->getGeometryInterpolation()->giveBoundarySurfaceIntegrationRule(o, boundaryID, this->giveGeometryType());
                this->integrateSurfaceTerm_c(answer, NTaTmTe(getT(), getT(), bl, boundaryID, 's'), ir2.get(), boundaryID, tStep);
            } else {
                std::unique_ptr<IntegrationRule> ir2 = this->getGeometryInterpolation()->giveBoundaryEdgeIntegrationRule(o, boundaryID, this->giveGeometryType());
                this->integrateEdgeTerm_c(answer, NTaTmTe(getT(), getT(), bl, boundaryID, 'e'), ir2.get(), boundaryID, tStep);
            }
        } else {
            answer.clear();
        }

    }

    void computeBoundarySurfaceLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global = true) override {
        answer.resize(this->getNumberOfSurfaceDOFs());
        answer.zero();
        
        if ( type != ExternalForcesVector ) {
            return;
        }

        bcType bct = load->giveType();
        if (bct == TransmissionBC ) {
        
            IntArray locu, loct;
            FloatArray contrib, contrib2;
            getSurfaceLocalCodeNumbers (locu, Variable::VariableQuantity::Displacement) ;
            getSurfaceLocalCodeNumbers (loct, Variable::VariableQuantity::Temperature) ;

            // integrate traction contribution (momentum balance)
            int o = getU()->interpolation->giveInterpolationOrder()+load->giveApproxOrder();
            std::unique_ptr<IntegrationRule> ir = this->getGeometryInterpolation()->giveBoundarySurfaceIntegrationRule(o, boundary, this->giveGeometryType());
            this->integrateSurfaceTerm_c(contrib, NTf_Surface(getU(), BoundaryFluxFunctor(load, boundary, getU()->dofIDs, 's'), boundary), ir.get(), boundary, tStep);
            answer.assemble(contrib, locu);

            // integrate mass (fluid) flux normal to the boundary (mass balance) 
            o = getT()->interpolation->giveInterpolationOrder()+load->giveApproxOrder();
            std::unique_ptr<IntegrationRule> ir2 = this->getGeometryInterpolation()->giveBoundarySurfaceIntegrationRule(o, boundary, this->giveGeometryType());
            this->integrateSurfaceTerm_c(contrib2, NTf_Surface(getT(), BoundaryFluxFunctor(load, boundary, getT()->dofIDs,'s'), boundary), ir2.get(), boundary, tStep);
            contrib2.times(-1.0);
            answer.assemble(contrib2, loct);

        } else if (bct == ConvectionBC) {
            // convection handled in residual evaluation
        } else {
            OOFEM_ERROR("Unsupported boundary condition type");
        }
    }

    void computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global=true) override {
        answer.resize(this->getNumberOfEdgeDOFs());
        answer.zero();
        
        if ( type != ExternalForcesVector ) {
            return;
        }
        
        bcType bct = load->giveType();
        if (bct == TransmissionBC ) {
            IntArray locu, loct;
            FloatArray contrib, contrib2;
            getEdgeLocalCodeNumbers (locu, Variable::VariableQuantity::Displacement) ;
            getEdgeLocalCodeNumbers (loct, Variable::VariableQuantity::Pressure) ;

            // integrate traction contribution (momentum balance)
            int o = getU()->interpolation->giveInterpolationOrder()+load->giveApproxOrder();
            std::unique_ptr<IntegrationRule> ir = this->getGeometryInterpolation()->giveBoundaryEdgeIntegrationRule(o, boundary, this->giveGeometryType());
            this->integrateEdgeTerm_c(contrib, NTf_Edge(getU(), BoundaryFluxFunctor(load, boundary, getU()->dofIDs,'e'), boundary), ir.get(), boundary, tStep);
            answer.assemble(contrib, locu);

            // integrate mass (fluid) flux normal to the boundary (mass balance) 
            o = getT()->interpolation->giveInterpolationOrder()+load->giveApproxOrder();
            std::unique_ptr<IntegrationRule> ir2 = this->getGeometryInterpolation()->giveBoundaryEdgeIntegrationRule(o, boundary, this->giveGeometryType());
            this->integrateEdgeTerm_c(contrib2, NTf_Edge(getT(), BoundaryFluxFunctor(load, boundary, getT()->dofIDs,'e'), boundary), ir2.get(), boundary, tStep);
            contrib2.times(-1.0);
            answer.assemble(contrib2, loct);
        } else if (bct == ConvectionBC) {
            // convection handled in residual evaluation
        } else {
            OOFEM_ERROR("Unsupported boundary condition type");
        }
    }

    void computeLoadVector(FloatArray &answer, BodyLoad *load, CharType type, ValueModeType mode, TimeStep *tStep) override {
        answer.resize(this->giveNumberOfDofs());
        answer.zero();
        if (type == ExternalForcesVector) {
            if (load->giveBCValType() == ForceLoadBVT) {
                FloatArray contribu, contribt;
                IntArray locu, loct;
                getLocalCodeNumbers (locu, Variable::VariableQuantity::Displacement) ;
                getLocalCodeNumbers (loct, Variable::VariableQuantity::Temperature) ;
                this->integrateTerm_c(contribu, NTf_Body(getU(), BodyFluxFunctor(load, getU()->dofIDs)), this->giveDefaultIntegrationRulePtr(), tStep);
                this->integrateTerm_c(contribt, NTf_Body(getT(), BodyFluxFunctor(load, getT()->dofIDs)), this->giveDefaultIntegrationRulePtr(), tStep);
                answer.assemble(contribu, locu);
                answer.assemble(contribt, loct);
            } else {
                OOFEM_ERROR("Unsupported body load type");
            }
        } else {
            OOFEM_ERROR("Unsupported load type");
        }
    }


    int computeFluxLBToLRotationMatrix(FloatMatrix &answer, int iSurf, const FloatArray& lc, const Variable::VariableQuantity q, char btype) override {
        if (q == Variable::VariableQuantity::Displacement) {
            // better to integrate this into FEInterpolation class 
            FloatArray nn, h1(3), h2(3);
            answer.resize(3,3);
            if (btype == 's') {
                this->getGeometryInterpolation()->boundarySurfaceEvalNormal(nn, iSurf, lc, FEIElementGeometryWrapper(this));
            } else {
                OOFEM_ERROR ("Unsupported boundary entity");
            }
            nn.normalize();
            for ( int i = 1; i <= 3; i++ ) {
                answer.at(i, 3) = nn.at(i);
            }   

            // determine lcs of surface
            // local x axis in xy plane
            double test = fabs(fabs( nn.at(3) ) - 1.0);
            if ( test < 1.e-5 ) {
                h1.at(1) = answer.at(1, 1) = 1.0;
                h1.at(2) = answer.at(2, 1) = 0.0;
            } else {
                h1.at(1) = answer.at(1, 1) = answer.at(2, 3);
                h1.at(2) = answer.at(2, 1) = -answer.at(1, 3);
            }

            h1.at(3) = answer.at(3, 1) = 0.0;
            // local y axis perpendicular to local x,z axes
            h2.beVectorProductOf(nn, h1);
            for ( int i = 1; i <= 3; i++ ) {
                answer.at(i, 2) = h2.at(i);
            }

            return 1;
        } else {
            answer.clear(); 
            return 0;
        }
    }
  //virtual void getLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q ) const = 0; 
    //virtual void giveDofManDofIDMask(int inode, IntArray &answer) const =0;

};



/**
 * @brief 3D Equal order linear Brick TS Element
 * 
 */
class TMBrick11 : public TMElement, public ZZNodalRecoveryModelInterface {
    protected:
        //FEI3dTetLin pInterpol;
        //FEI3dTetQuad uInterpol;
        const static FEI3dHexaLin  tInterpol;
        const static FEI3dHexaLin  uInterpol;
        const static Variable t;
        const static Variable u;
      
    public:
    TMBrick11(int n, Domain* d): 
        TMElement(n,d), ZZNodalRecoveryModelInterface(this)
    {
        numberOfDofMans  = 8;
        numberOfGaussPoints = 8;
        this->computeGaussPoints();
    }

  void getDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {
        /* dof ordering: u1 v1 w1 p1  u2 v2 w2 p2  u3 v3 w3 p3  u4 v4 w4   u5 v5 w5  u6 v6 w6*/
        if (q == Variable::VariableQuantity::Displacement) {
          //answer={1,2,3, 5,6,7, 9,10,11, 13,14,15, 17,18,19, 21,22,23, 25,26,27, 29,30,31 };
          int o = (num-1)*4+1;
          answer={o, o+1, o+2};
        } else if (q == Variable::VariableQuantity::Temperature) {
          //answer = {4, 8, 12, 16, 20, 24, 28, 32};
          answer={num*4};
        }
    }
    void getInternalDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {
        answer={};
    }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override { 
            answer = {1,2,3,10};
    }
    int giveNumberOfDofs() override { return 32; }
    const char *giveInputRecordName() const override {return "tmbrick11";}
    const char *giveClassName() const override { return "TMBrick11"; }

    
    const FEInterpolation* getGeometryInterpolation() const override {return &this->tInterpol;}
  
    Element_Geometry_Type giveGeometryType() const override {
        return EGT_hexa_1;
    }
    int getNumberOfSurfaceDOFs() const override {return 16;}
    int getNumberOfEdgeDOFs() const override {return 0;}
    void getSurfaceLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {
        if (q == Variable::VariableQuantity::Displacement) {
        answer={1,2,3, 5,6,7, 9,10,11, 13,14,15};
        } else {
        answer ={4, 8, 12, 16};
        }
    }
    void getEdgeLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {}
    Interface *giveInterface(InterfaceType it) override {
        if (it == ZZNodalRecoveryModelInterfaceType) {
            return this;
        } else {
            return NULL;
        }
    }



private:
        virtual int  giveNumberOfUDofs() const override {return 24;} 
        virtual int  giveNumberOfTDofs() const override {return 8;}
        virtual const Variable* getU() const override {return &u;}
        virtual const Variable* getT() const override {return &t;}
        void computeGaussPoints() override {
            if ( integrationRulesArray.size() == 0 ) {
                integrationRulesArray.resize( 1 );
                integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this);
                integrationRulesArray [ 0 ]->SetUpPointsOnCube(numberOfGaussPoints, _3dMat);
            }
        }
};

const FEI3dHexaLin TMBrick11::uInterpol;
const FEI3dHexaLin TMBrick11::tInterpol;
const Variable TMBrick11::t(&TMBrick11::tInterpol, Variable::VariableQuantity::Temperature, Variable::VariableType::scalar, 1, NULL, {10});
const Variable TMBrick11::u(&TMBrick11::uInterpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 3, NULL, {1,2,3});

#define _IFT_TMBrick11_Name "tmbrick11"
REGISTER_Element(TMBrick11)


/**
 * @brief 3D Equal order linear Tetrahedra TS Element
 * 
 */
class TMTetra11 : public TMElement, public ZZNodalRecoveryModelInterface {
    protected:
        //FEI3dTetLin pInterpol;
        //FEI3dTetQuad uInterpol;
        const static FEI3dTetLin tInterpol;
        const static FEI3dTetLin uInterpol;
        const static Variable t;
        const static Variable u;
      
    public:
    TMTetra11(int n, Domain* d): 
        TMElement(n,d), ZZNodalRecoveryModelInterface(this)
    {
        numberOfDofMans  = 4;
        numberOfGaussPoints = 4;
        this->computeGaussPoints();
    }

  void getDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {
        /* dof ordering: u1 v1 w1 p1  u2 v2 w2 p2  u3 v3 w3 p3  u4 v4 w4   u5 v5 w5  u6 v6 w6*/
        if (q == Variable::VariableQuantity::Displacement) {
          //answer={1,2,3, 5,6,7, 9,10,11, 13,14,15};
          int o = (num-1)*4+1;
          answer={o, o+1, o+2};
        } else if (q == Variable::VariableQuantity::Temperature) {
          //answer = {4, 8, 12, 16};
          answer={num*4};
        }
    }
    void getInternalDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {
        answer={};
    }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override { 
            answer = {1,2,3,10};
    }
    int giveNumberOfDofs() override { return 16; }
    const char *giveInputRecordName() const override {return "tmtetra11";}
    const char *giveClassName() const override { return "TMTetra11"; }

    
    const FEInterpolation* getGeometryInterpolation() const override {return &this->tInterpol;}
  
    Element_Geometry_Type giveGeometryType() const override {
        return EGT_tetra_1;
    }
    int getNumberOfSurfaceDOFs() const override {return 12;}
    int getNumberOfEdgeDOFs() const override {return 0;}
    void getSurfaceLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {
        if (q == Variable::VariableQuantity::Displacement) {
        answer={1,2,3, 5,6,7, 9,10,11};
        } else {
        answer ={4, 8, 12};
        }
    }
    void getEdgeLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {}
    Interface *giveInterface(InterfaceType it) override {
        if (it == ZZNodalRecoveryModelInterfaceType) {
            return this;
        } else {
            return NULL;
        }
    }



private:
        virtual int  giveNumberOfUDofs() const override {return 12;} 
        virtual int  giveNumberOfTDofs() const override {return 4;}
        virtual const Variable* getU() const override {return &u;}
        virtual const Variable* getT() const override {return &t;}
        void computeGaussPoints() override {
            if ( integrationRulesArray.size() == 0 ) {
                integrationRulesArray.resize( 1 );
                integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this);
                integrationRulesArray [ 0 ]->SetUpPointsOnTetrahedra(numberOfGaussPoints, _3dMat);
            }
        }
};

const FEI3dTetLin TMTetra11::uInterpol;
const FEI3dTetLin TMTetra11::tInterpol;
const Variable TMTetra11::t(&TMTetra11::tInterpol, Variable::VariableQuantity::Temperature, Variable::VariableType::scalar, 1, NULL, {10});
const Variable TMTetra11::u(&TMTetra11::uInterpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 3, NULL, {1,2,3});

#define _IFT_TMTetra11_Name "tmtetra11"
REGISTER_Element(TMTetra11)


#define _IFT_TMSimpleMaterial_Name "tmm"
#define _IFT_TMSimpleMaterial_E "e"
#define _IFT_TMSimpleMaterial_nu "nu"
#define _IFT_TMSimpleMaterial_lambda "lambda"
#define _IFT_TMSimpleMaterial_alpha "alpha"
#define _IFT_TMSimpleMaterial_c "c"

class TMMaterialStatus : public MaterialStatus
{
protected:
    /// Equilibrated strain vector in reduced form
    FloatArray strainVector;
    /// Equilibrated stress vector in reduced form
    FloatArray stressVector;
    /// Temporary stress vector in reduced form (increments are used mainly in nonlinear analysis)
    FloatArray tempStressVector;
    /// Temporary strain vector in reduced form (to find balanced state)
    FloatArray tempStrainVector;
    /// Temporary flux 
    FloatArray tempFluxVector;
    /// Equilibrated flux
    FloatArray fluxVector;
public:
    /// Constructor. Creates new StructuralMaterialStatus with IntegrationPoint g.
    TMMaterialStatus (GaussPoint * g) : MaterialStatus(g), strainVector(), stressVector(),
    tempStressVector(), tempStrainVector() 
    {}

/// Returns the const pointer to receiver's strain vector.
    const FloatArray &giveStrainVector() const { return strainVector; }
    /// Returns the const pointer to receiver's stress vector.
    const FloatArray &giveStressVector() const { return stressVector; }
    /// Returns the const pointer to receiver's temporary strain vector.
    const FloatArray &giveTempStrainVector() const { return tempStrainVector; }
    /// Returns the const pointer to receiver's temporary stress vector.
    const FloatArray &giveTempStressVector() const { return tempStressVector; }
    /// Returns the const pointer to receiver's temporary flux vector.
    const FloatArray &giveTempFluxVector() const { return tempFluxVector; }
    /// Returns the const pointer to receiver's flux vector.
    const FloatArray &giveFluxVector() const { return fluxVector; }
    /// Assigns tempStressVector to given vector v.
    void letTempStressVectorBe(const FloatArray &v) { tempStressVector = v; }
    /// Assigns tempStrainVector to given vector v
    void letTempStrainVectorBe(const FloatArray &v) { tempStrainVector = v; }
    /// Assigns tempFluxVector to given vector v
    void letTempFluxVectorBe(const FloatArray &v) { tempFluxVector = v; }

    void printOutputAt(FILE *file, TimeStep *tStep) const override {
        MaterialStatus :: printOutputAt(file, tStep);
        fprintf(file, "  strains ");
        for ( auto &var : strainVector ) {
            fprintf( file, " %+.4e", var );
        }
      
        fprintf(file, "\n              stresses");  
        for ( auto &var : stressVector ) {
            fprintf( file, " %+.4e", var );
        }

        fprintf(file, "\n              fluxes");
        for ( auto &var : fluxVector ) {
            fprintf( file, " %+.4e", var );
        }
        fprintf(file, "\n");
    }

    void initTempStatus() override {
        MaterialStatus :: initTempStatus();
        tempStressVector = stressVector;
        tempStrainVector = strainVector;
        tempFluxVector = fluxVector;
    }
    void updateYourself(TimeStep *tStep) override {
        MaterialStatus :: updateYourself(tStep);
        stressVector = tempStressVector;
        strainVector = tempStrainVector;
        fluxVector = tempFluxVector;
    }
    const char *giveClassName() const override {return "TMMaterialStatus";}

};

class TMSimpleMaterial : public Material {
    protected:
        double e, nu; // elastic isotropic constants
        double lambda; // isotropic conductivity
        double alpha; // thermal expansion coefficient
        double c; // thermal capacity
    public:
  TMSimpleMaterial (int n, Domain* d) : Material (n,d) {e=1.0; nu=0.15; lambda=1.0; alpha=1.0; c=0.1;}

    void giveCharacteristicMatrix(FloatMatrix &answer, MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const override {
        MaterialMode mmode = gp->giveMaterialMode();
        if (type == TangentStiffness) {
            double ee;

            ee = e / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) =  1. - nu;
            answer.at(1, 2) =  nu;
            answer.at(1, 3) =  nu;
            answer.at(2, 1) =  nu;
            answer.at(2, 2) =  1. - nu;
            answer.at(2, 3) =  nu;
            answer.at(3, 1) =  nu;
            answer.at(3, 2) =  nu;
            answer.at(3, 3) =  1. - nu;

            answer.at(4, 4) =  ( 1. - 2. * nu ) * 0.5;
            answer.at(5, 5) =  ( 1. - 2. * nu ) * 0.5;
            answer.at(6, 6) =  ( 1. - 2. * nu ) * 0.5;

            answer.times(ee);
        } else if (type == DSigmaDT) {
            answer.resize(6,1);
            answer.zero();
        } else if (type == Conductivity) {
            if (mmode == _3dMat) {
                answer.resize(3,3);
                answer.beUnitMatrix();
                answer.times(this->lambda);
            }
        } else {
            OOFEM_ERROR("Unknown characteristic matrix type");
        }
    }
    /**
     * @param flux Generalized strain vector, flux.at(1-6) containing total strain vector, flux(7-9) temperature gradient, flux(10) temperature
     */
    void giveCharacteristicVector(FloatArray &answer, FloatArray& flux, MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const override {
        TMMaterialStatus *status = static_cast< TMMaterialStatus * >( this->giveStatus(gp) );
        if (type == Stress) {
            FloatMatrix d;
            FloatArray eps(6);
            for (int i=0; i<6; i++) {
                eps(i) = flux(i);
            }
            double t = flux(9);
            eps(0)-= t*alpha;
            eps(1)-= t*alpha;
            eps(2)-= t*alpha;
 
            this->giveCharacteristicMatrix(d, TangentStiffness, gp, tStep);

            answer.beProductOf(d, eps);
            // update gp status
            status->letTempStrainVectorBe(flux);
            status->letTempStressVectorBe(answer);

        } else if (type == Flux) {
            FloatMatrix k;
            FloatArray grad(3);
            this->giveCharacteristicMatrix(k, Conductivity, gp, tStep);
            grad(0) = -flux(6);
            grad(1) = -flux(7);
            grad(2) = -flux(8);
            answer.beProductOf(k, grad);
            status->letTempFluxVectorBe(answer);
        } else if (type == IntSource) {
            answer.resize(1);
            answer.zero();
        } else {
            OOFEM_ERROR("Unknown characteristic vector type");
        }
    }


    double giveCharacteristicValue(MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const override {
        if (type == BiotConstant) {
            return alpha;
        } else if (type == Capacity) {
            return c;
        } else {
            return 0.0;
        }
    };

    void initializeFrom(InputRecord &ir) override {
        Material :: initializeFrom(ir);

        IR_GIVE_OPTIONAL_FIELD(ir, e, _IFT_TMSimpleMaterial_E);
        IR_GIVE_OPTIONAL_FIELD(ir, nu, _IFT_TMSimpleMaterial_nu);
        IR_GIVE_OPTIONAL_FIELD(ir, lambda, _IFT_TMSimpleMaterial_lambda);
        IR_GIVE_OPTIONAL_FIELD(ir, alpha, _IFT_TMSimpleMaterial_alpha);
        IR_GIVE_OPTIONAL_FIELD(ir, c, _IFT_TMSimpleMaterial_c);

    };
    //   void giveInputRecord(DynamicInputRecord &input) override {};
    void giveInputRecord(DynamicInputRecord &input) override {};
    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<TMMaterialStatus>(gp); }

    const char *giveClassName() const override {return "TMSimpleMaterial";}
    const char *giveInputRecordName() const override {return _IFT_TMSimpleMaterial_Name;}
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override {
        TMMaterialStatus *status = static_cast< TMMaterialStatus * >( this->giveStatus(gp) );
        if ( type == IST_StrainTensor ) {
            answer = status->giveStrainVector(); 
            return 1;
        }
        if ( type == IST_StressTensor ) {
            answer = status->giveStressVector();
            return 1;
        } else {
        return Material::giveIPValue(answer, gp, type, tStep);
        }
    }

};
REGISTER_Material(TMSimpleMaterial)

} // end namespace oofem
