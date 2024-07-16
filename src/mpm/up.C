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


#include "mpm.h"
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
#include "zznodalrecoverymodel.h"

namespace oofem {


/**
 * @brief Base class for (3D) UP elements
 * 
 */
class UPElement : public MPElement {
        
    public:
    UPElement(int n, Domain* d): 
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
            this->integrateTerm_dw (answer, BTSigTerm(getU(),getU()), ir, tStep) ;
        } else if (type == MomentumBalance_PressureCouplingMatrix) {
            answer.resize(this->giveNumberOfUDofs(),this->giveNumberOfPDofs());
            answer.zero();
            this->integrateTerm_dw (answer, BTamNTerm(getU(),getP()), ir, tStep) ;
        } else if (type == MassBalance_PermeabilityMatrix) {
            int pdofs = this->giveNumberOfPDofs();
            answer.resize(pdofs,pdofs);
            answer.zero();
            this->integrateTerm_dw (answer, gNTfTerm(getP(), getP(), Permeability, FluidMassBalancePressureContribution), ir, tStep) ;
        } else if (type == MassBalance_CompresibilityMatrix) {
            int pdofs = this->giveNumberOfPDofs();
            answer.resize(pdofs,pdofs);
            answer.zero();
            this->integrateTerm_dw (answer, NTcN(getP(), getP(), CompressibilityCoefficient), ir, tStep) ;
        } else if (type == MassBalance_StressCouplingMatrix) {
            answer.resize(this->giveNumberOfPDofs(),this->giveNumberOfUDofs());
            answer.zero();
            this->integrateTerm_dw (answer, NTamTBTerm(getP(), getU()), ir, tStep) ;
        } else {
	        OOFEM_ERROR("Unknown characteristic matrix type");
	    }
    }

    void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) override {
        IntegrationRule* ir = this->giveDefaultIntegrationRulePtr();
        if (type == MomentumBalance_StressResidual) {
            answer.resize(this->giveNumberOfUDofs());
            answer.zero();
            this->integrateTerm_c (answer, BTSigTerm(getU(),getU()), ir, tStep) ;
        } else if (type == MomentumBalance_PressureResidual) {
            answer.resize(this->giveNumberOfUDofs());
            answer.zero();
            this->integrateTerm_c(answer, BTamNTerm(getU(),getP()), ir, tStep) ;
        } else if (type == MassBalance_StressRateResidual) {
            answer.resize(this->giveNumberOfPDofs());
            answer.zero();
            this->integrateTerm_c (answer, NTamTBTerm(getP(), getU()), ir, tStep) ;
        } else if (type == MassBalance_PressureResidual) {
            answer.resize(this->giveNumberOfPDofs());
            answer.zero();
            this->integrateTerm_c (answer, gNTfTerm(getP(), getP(), Permeability, FluidMassBalancePressureContribution), ir, tStep) ;
        } else if (type == MassBalance_PressureRateResidual) {
            answer.resize(this->giveNumberOfPDofs());
            answer.zero();
            this->integrateTerm_c (answer, NTcN(getP(), getP(), CompressibilityCoefficient), ir, tStep) ;   
        } else if (type == ExternalForcesVector) {
          answer.zero();
        } else {
	        OOFEM_ERROR("Unknown characteristic vector type");
	    }
    }

    void computeBoundarySurfaceLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global = true) override {
        answer.resize(giveNumberOfDofs());
        answer.zero();
        if ( type != ExternalForcesVector ) {
            return;
        }

        IntArray locu, locp;
        FloatArray contrib, contrib2;
        getSurfaceLocalCodeNumbers (locu, Variable::VariableQuantity::Displacement) ;
        getSurfaceLocalCodeNumbers (locp, Variable::VariableQuantity::Pressure) ;

        // integrate traction contribution (momentum balance)
        int o = getU().interpolation.giveInterpolationOrder()+load->giveApproxOrder();
        std::unique_ptr<IntegrationRule> ir = this->getGeometryInterpolation().giveBoundarySurfaceIntegrationRule(o, boundary, this->giveGeometryType());
        this->integrateSurfaceTerm_c(contrib, NTf_Surface(getU(), BoundaryFluxFunctor(load, boundary, getU().dofIDs, 's'), boundary), ir.get(), boundary, tStep);

        answer.resize(this->getNumberOfSurfaceDOFs());
        answer.zero();
        answer.assemble(contrib, locu);

        // integrate mass (fluid) flux normal to the boundary (mass balance) 
        o = getP().interpolation.giveInterpolationOrder()+load->giveApproxOrder();
        std::unique_ptr<IntegrationRule> ir2 = this->getGeometryInterpolation().giveBoundarySurfaceIntegrationRule(o, boundary, this->giveGeometryType());
        this->integrateSurfaceTerm_c(contrib2, NTf_Surface(getP(), BoundaryFluxFunctor(load, boundary, getP().dofIDs,'s'), boundary), ir2.get(), boundary, tStep);
        answer.assemble(contrib2, locp);
    }

    void computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global=true) override {
        answer.resize(giveNumberOfDofs());
        answer.zero();
        if ( type != ExternalForcesVector ) {
            return;
        }

        IntArray locu, locp;
        FloatArray contrib, contrib2;
        getEdgeLocalCodeNumbers (locu, Variable::VariableQuantity::Displacement) ;
        getEdgeLocalCodeNumbers (locp, Variable::VariableQuantity::Pressure) ;

        // integrate traction contribution (momentum balance)
        int o = getU().interpolation.giveInterpolationOrder()+load->giveApproxOrder();
        std::unique_ptr<IntegrationRule> ir = this->getGeometryInterpolation().giveBoundaryEdgeIntegrationRule(o, boundary, this->giveGeometryType());
        this->integrateEdgeTerm_c(contrib, NTf_Edge(getU(), BoundaryFluxFunctor(load, boundary, getU().dofIDs,'e'), boundary), ir.get(), boundary, tStep);

        answer.resize(this->getNumberOfEdgeDOFs());
        answer.zero();
        answer.assemble(contrib, locu);

        // integrate mass (fluid) flux normal to the boundary (mass balance) 
        o = getP().interpolation.giveInterpolationOrder()+load->giveApproxOrder();
        std::unique_ptr<IntegrationRule> ir2 = this->getGeometryInterpolation().giveBoundaryEdgeIntegrationRule(o, boundary, this->giveGeometryType());
        this->integrateEdgeTerm_c(contrib2, NTf_Edge(getP(), BoundaryFluxFunctor(load, boundary, getP().dofIDs,'e'), boundary), ir2.get(), boundary, tStep);
        answer.assemble(contrib2, locp);
    }


    int computeFluxLBToLRotationMatrix(FloatMatrix &answer, int iSurf, const FloatArray& lc, const Variable::VariableQuantity q, char btype) override {
        if (q == Variable::VariableQuantity::Displacement) {
            // better to integrate this into FEInterpolation class 
            FloatArray nn, h1(3), h2(3);
            answer.resize(3,3);
            if (btype == 's') {
                this->getGeometryInterpolation().boundarySurfaceEvalNormal(nn, iSurf, lc, FEIElementGeometryWrapper(this));
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
    private:
        virtual int  giveNumberOfUDofs() const = 0;
        virtual int  giveNumberOfPDofs() const = 0;
        virtual const Variable& getU() const = 0;
        virtual const Variable& getP() const = 0;
};

/**
 * @brief 3D Tetrahedra element with quadratic interpolation for displacements, linear interpolation for pressure
 * 
 */
class UPTetra21 : public UPElement {
    protected:
        //FEI3dTetLin pInterpol;
        //FEI3dTetQuad uInterpol;
        const static FEInterpolation & pInterpol;
        const static FEInterpolation & uInterpol;
        const static Variable& p;
        const static Variable& u;

      
    public:
    UPTetra21(int n, Domain* d): 
        UPElement(n,d) 
    {
        numberOfDofMans  = 10;
        numberOfGaussPoints = 4;
        this->computeGaussPoints();
    }

  void getDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {
        /* dof ordering: u1 v1 w1 p1  u2 v2 w2 p2  u3 v3 w3 p3  u4 v4 w4   u5 v5 w5  u6 v6 w6*/
        if (q == Variable::VariableQuantity::Displacement) {
          //answer={1,2,3, 5,6,7, 9,10,11, 13,14,15, 17,18,19, 20,21,22, 23,24,25, 26,27,28, 29,30,31, 32,33,34 };
          int o = (num-1)*4+1-(num>4)*(num-5);
          answer = {o, o+1, o+2};
        } else if (q == Variable::VariableQuantity::Pressure) {
          if (num<=4) {
            //answer = {4, 8, 12, 16};
            answer={num*4};
          } else {
            answer={};
          }
        }
    }
    void getInternalDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {
        answer={};
    }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override { 
        if (inode >0 && inode <5) {
            answer = {1,2,3,11};
        } else {
            answer= {1,2,3};
        }
    }
    int giveNumberOfDofs() override { return 34; }
    const char *giveInputRecordName() const override {return "uptetra21";}
    const FEInterpolation& getGeometryInterpolation() const override {return this->uInterpol;}
  
    Element_Geometry_Type giveGeometryType() const override {
        return EGT_tetra_2;
    }
    int getNumberOfSurfaceDOFs() const override {return 21;}
    int getNumberOfEdgeDOFs() const override {return 0;}
    void getSurfaceLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {
        if (q == Variable::VariableQuantity::Displacement) {
        answer={1,2,3, 5,6,7, 9,10,11, 13,14,15, 16,17,18, 19,20,21};
        } else {
        answer ={4, 8, 12};
        }
    }
  void getEdgeLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {}

    private:
        virtual int  giveNumberOfUDofs() const override {return 30;} 
        virtual int  giveNumberOfPDofs() const override {return 4;}
        virtual const Variable& getU() const override {return u;}
        virtual const Variable& getP() const override {return p;}
        void computeGaussPoints() override {
            if ( integrationRulesArray.size() == 0 ) {
                integrationRulesArray.resize( 1 );
                integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this);
                integrationRulesArray [ 0 ]->SetUpPointsOnTetrahedra(numberOfGaussPoints, _3dUP);
            }
        }
};

const FEInterpolation & UPTetra21::uInterpol = FEI3dTetQuad();
const FEInterpolation & UPTetra21::pInterpol = FEI3dTetLin();
const Variable& UPTetra21::p = Variable(UPTetra21::pInterpol, Variable::VariableQuantity::Pressure, Variable::VariableType::scalar, 1, NULL, {11});
const Variable& UPTetra21::u = Variable(UPTetra21::uInterpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 3, NULL, {1,2,3});

#define _IFT_UPTetra21_Name "uptetra21"
REGISTER_Element(UPTetra21)

/**
 * @brief 3D Equal order linear Brick UP Element
 * 
 */
class UPBrick11 : public UPElement, public ZZNodalRecoveryModelInterface {
    protected:
        //FEI3dTetLin pInterpol;
        //FEI3dTetQuad uInterpol;
        const static FEInterpolation & pInterpol;
        const static FEInterpolation & uInterpol;
        const static Variable& p;
        const static Variable& u;
      
    public:
    UPBrick11(int n, Domain* d): 
        UPElement(n,d), ZZNodalRecoveryModelInterface(this)
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
        } else if (q == Variable::VariableQuantity::Pressure) {
          //answer = {4, 8, 12, 16, 20, 24, 28, 32};
          answer={num*4};
        }
    }
    void getInternalDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {
        answer={};
    }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override { 
            answer = {1,2,3,11};
    }
    int giveNumberOfDofs() override { return 32; }
    const char *giveInputRecordName() const override {return "upbrick11";}
    const char *giveClassName() const override { return "UPBrick11"; }

    
    const FEInterpolation& getGeometryInterpolation() const override {return this->pInterpol;}
  
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
        virtual int  giveNumberOfPDofs() const override {return 8;}
        virtual const Variable& getU() const override {return u;}
        virtual const Variable& getP() const override {return p;}
        void computeGaussPoints() override {
            if ( integrationRulesArray.size() == 0 ) {
                integrationRulesArray.resize( 1 );
                integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this);
                integrationRulesArray [ 0 ]->SetUpPointsOnCube(numberOfGaussPoints, _3dUP);
            }
        }
};

const FEInterpolation & UPBrick11::uInterpol = FEI3dHexaLin();
const FEInterpolation & UPBrick11::pInterpol = FEI3dHexaLin();
const Variable& UPBrick11::p = Variable(UPBrick11::pInterpol, Variable::VariableQuantity::Pressure, Variable::VariableType::scalar, 1, NULL, {11});
const Variable& UPBrick11::u = Variable(UPBrick11::uInterpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 3, NULL, {1,2,3});

#define _IFT_UPBrick11_Name "upbrick11"
REGISTER_Element(UPBrick11)

/**
 * @brief 2D Equal order linear Quad UP Element
 * 
 */
class UPQuad11 : public UPElement {
    protected:
        //FEI3dTetLin pInterpol;
        //FEI3dTetQuad uInterpol;
        const static FEInterpolation & pInterpol;
        const static FEInterpolation & uInterpol;
        const static Variable& p;
        const static Variable& u;
       
    public:
    UPQuad11(int n, Domain* d): 
        UPElement(n,d)
    {
        numberOfDofMans  = 4;
        numberOfGaussPoints = 4;
        this->computeGaussPoints();
    }

  void getDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {
        /* dof ordering: u1 v1 w1 p1  u2 v2 w2 p2  u3 v3 w3 p3  u4 v4 w4 p4*/
        if (q == Variable::VariableQuantity::Displacement) {
          //answer={1,2,3, 5,6,7, 9,10,11, 13,14,15 };
          int o = (num-1)*3+1;
          answer={o, o+1};
        } else if (q == Variable::VariableQuantity::Pressure) {
          //answer = {4, 8, 12, 16};
          answer={num*3};
        }
    }
    void getInternalDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {
        answer={};
    }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override { 
            answer = {1,2,11};
    }
    int giveNumberOfDofs() override { return 12; }
    const char *giveInputRecordName() const override {return "upquad11";}
    
    const FEInterpolation& getGeometryInterpolation() const override {return this->pInterpol;}
  
    Element_Geometry_Type giveGeometryType() const override {
        return EGT_quad_1;
    }
    int getNumberOfSurfaceDOFs() const override {return 0;}
    void getSurfaceLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {
        answer={};
    }

    int getNumberOfEdgeDOFs() const override  {return 6;}
    void getEdgeLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override  {
        if (q == Variable::VariableQuantity::Displacement) {
            answer={1,2, 4,5};
        } else {
            answer ={3, 6};
        }
    }
  


private:
        virtual int  giveNumberOfUDofs() const override {return 8;} 
        virtual int  giveNumberOfPDofs() const override {return 4;}
        virtual const Variable& getU() const override {return u;}
        virtual const Variable& getP() const override {return p;}
        void computeGaussPoints() override {
            if ( integrationRulesArray.size() == 0 ) {
                integrationRulesArray.resize( 1 );
                integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this);
                integrationRulesArray [ 0 ]->SetUpPointsOnSquare(numberOfGaussPoints, _2dUP);
            }
        }
};

const FEInterpolation & UPQuad11::uInterpol = FEI2dQuadLin(1,2);
const FEInterpolation & UPQuad11::pInterpol = FEI2dQuadLin(1,2);
const Variable& UPQuad11::p = Variable(UPQuad11::pInterpol, Variable::VariableQuantity::Pressure, Variable::VariableType::scalar, 1, NULL, {11});
const Variable& UPQuad11::u = Variable(UPQuad11::uInterpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 2, NULL, {1,2});

#define _IFT_UPQuad11_Name "upquad11"
REGISTER_Element(UPQuad11)



#define _IFT_UPSimpleMaterial_Name "upm"
#define _IFT_UPSimpleMaterial_E "e"
#define _IFT_UPSimpleMaterial_nu "nu"
#define _IFT_UPSimpleMaterial_k "k"
#define _IFT_UPSimpleMaterial_alpha "alpha"
#define _IFT_UPSimpleMaterial_c "c"

class UPMaterialStatus : public MaterialStatus
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
public:
    /// Constructor. Creates new StructuralMaterialStatus with IntegrationPoint g.
    UPMaterialStatus (GaussPoint * g) : MaterialStatus(g), strainVector(), stressVector(),
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
    /// Assigns tempStressVector to given vector v.
    void letTempStressVectorBe(const FloatArray &v) { tempStressVector = v; }
    /// Assigns tempStrainVector to given vector v
    void letTempStrainVectorBe(const FloatArray &v) { tempStrainVector = v; }

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
        fprintf(file, "\n");
    }

    void initTempStatus() override {
        MaterialStatus :: initTempStatus();
        tempStressVector = stressVector;
        tempStrainVector = strainVector;
    }
    void updateYourself(TimeStep *tStep) override {
        MaterialStatus :: updateYourself(tStep);
        stressVector = tempStressVector;
        strainVector = tempStrainVector;
    }
    const char *giveClassName() const override {return "UPMaterialStatus";}

};

class UPSimpleMaterial : public Material {
    protected:
        double e, nu; // elastic isotropic constants
        double k; // isotropic permeability
        double alpha; // Biot constant = 1-K_t/K_s (Kt bulk moduli of the porous medium, Ks bulk moduli of solid phase)
        double c; // 1/Q, where Q is combined compressibility of the fluid and solid phases (1/Q=n/Kt+(b-n)/Ks, where n is porosity)
        double muw; // dynamic viscosity of water
    public:
  UPSimpleMaterial (int n, Domain* d) : Material (n,d) {e=1.0; nu=0.15; k=1.0; alpha=1.0; c=0.1; muw = 1.0;}

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
        } else if (type == Permeability) {
            if (mmode == _3dUP) {
                answer.resize(3,3);
                answer.beUnitMatrix();
                answer.times(this->k/this->muw);
            } else if (mmode == _2dUP) {
                answer.resize(2,2);
                answer.beUnitMatrix();
                answer.times(this->k);
            }
        }
    }

    void giveCharacteristicVector(FloatArray &answer, FloatArray& flux, MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const override {
        if (type == Stress) {
            FloatMatrix d;
            UPMaterialStatus *status = static_cast< UPMaterialStatus * >( this->giveStatus(gp) );

            this->giveCharacteristicMatrix(d, TangentStiffness, gp, tStep);
            answer.beProductOf(d, flux);
            // update gp status
            status->letTempStrainVectorBe(flux);
            status->letTempStressVectorBe(answer);

        }else if (type == FluidMassBalancePressureContribution) {
            FloatMatrix k;
            this->giveCharacteristicMatrix(k, Permeability, gp, tStep);
            answer.beProductOf(k, flux);
        }
    }


    double giveCharacteristicValue(MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const override {
        if (type == BiotConstant) {
            return alpha;
        } else if (type == CompressibilityCoefficient) {
            return c;
        } else {
            return 0.0;
        }
    };

    void initializeFrom(InputRecord &ir) override {
        Material :: initializeFrom(ir);

        IR_GIVE_OPTIONAL_FIELD(ir, e, _IFT_UPSimpleMaterial_E);
        IR_GIVE_OPTIONAL_FIELD(ir, nu, _IFT_UPSimpleMaterial_nu);
        IR_GIVE_OPTIONAL_FIELD(ir, k, _IFT_UPSimpleMaterial_k);
        IR_GIVE_OPTIONAL_FIELD(ir, alpha, _IFT_UPSimpleMaterial_alpha);
        IR_GIVE_OPTIONAL_FIELD(ir, c, _IFT_UPSimpleMaterial_c);

    };
    //   void giveInputRecord(DynamicInputRecord &input) override {};
    void giveInputRecord(DynamicInputRecord &input) override {};
    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new UPMaterialStatus(gp); }

    const char *giveClassName() const override {return "UPSimpleMaterial";}
    const char *giveInputRecordName() const override {return _IFT_UPSimpleMaterial_Name;}
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override {
        UPMaterialStatus *status = static_cast< UPMaterialStatus * >( this->giveStatus(gp) );
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
REGISTER_Material(UPSimpleMaterial)

} // end namespace oofem
