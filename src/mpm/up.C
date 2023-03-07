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
#include "mathfem.h"

#include "material.h"

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
            this->integrateTerm_dw (answer, gNTfTerm(getP(), getP()), ir, tStep) ;
        } else if (type == MassBalance_CompresibilityMatrix) {
            int pdofs = this->giveNumberOfPDofs();
            answer.resize(pdofs,pdofs);
            answer.zero();
            this->integrateTerm_dw (answer, NTcN(getP(), getP()), ir, tStep) ;
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
            this->integrateTerm_c (answer, gNTfTerm(getP(), getP()), ir, tStep) ;
        } else if (type == MassBalance_PressureRateResidual) {
            answer.resize(this->giveNumberOfPDofs());
            answer.zero();
            this->integrateTerm_c (answer, NTcN(getP(), getP()), ir, tStep) ;   
        } else {
	        OOFEM_ERROR("Unknown characteristic vector type");
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

        GaussIntegrationRule ir;
        
    public:
    UPTetra21(int n, Domain* d): 
        UPElement(n,d), 
        ir(1, this)
    {
        numberOfDofMans  = 10;
        numberOfGaussPoints = 4;
        ir.SetUpPointsOnTetrahedra(numberOfGaussPoints, _Unknown);
    }

    void getLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q ) const  override {
        /* dof ordering: u1 v1 w1 p1  u2 v2 w2 p2  u3 v3 w3 p3  u4 v4 w4   u5 v5 w5  u6 v6 w6*/
        if (q == Variable::VariableQuantity::Displacement) {
            answer={1,2,3, 5,6,7, 9,10,11, 13,14,15, 17,18,19, 20,21,22, 23,24,25, 26,27,28, 29,30,31, 32,33,34 };
        } else if (q == Variable::VariableQuantity::Pressure) {
            answer = {4, 8, 12, 16};
        }
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
    
    double computeVolumeAround(GaussPoint *gp) override {
        double determinant = fabs( this->pInterpol.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
        double weight = gp->giveWeight();
        return determinant * weight ;
    }
    IntegrationRule *giveDefaultIntegrationRulePtr() override {
        return &ir;
    }
    Element_Geometry_Type giveGeometryType() const override {
        return EGT_tetra_2;
    }
    private:
        virtual int  giveNumberOfUDofs() const override {return 30;} 
        virtual int  giveNumberOfPDofs() const override {return 4;}
        virtual const Variable& getU() const override {return u;}
        virtual const Variable& getP() const override {return p;}
};

const FEInterpolation & UPTetra21::uInterpol = FEI3dTetQuad();
const FEInterpolation & UPTetra21::pInterpol = FEI3dTetLin();
const Variable& UPTetra21::p = Variable(UPTetra21::pInterpol, Variable::VariableQuantity::Pressure, Variable::VariableType::scalar, 3, NULL, {11});
const Variable& UPTetra21::u = Variable(UPTetra21::uInterpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 3, NULL, {1,2,3});

#define _IFT_UPTetra21_Name "uptetra21"
REGISTER_Element(UPTetra21)

/**
 * @brief 3D Equal order linear Brick UP Element
 * 
 */
class UPBrick11 : public UPElement {
    protected:
        //FEI3dTetLin pInterpol;
        //FEI3dTetQuad uInterpol;
        const static FEInterpolation & pInterpol;
        const static FEInterpolation & uInterpol;
        const static Variable& p;
        const static Variable& u;

        GaussIntegrationRule ir;
        
    public:
    UPBrick11(int n, Domain* d): 
        UPElement(n,d), 
        ir(1, this)
    {
        numberOfDofMans  = 8;
        numberOfGaussPoints = 4;
        ir.SetUpPointsOnCube(numberOfGaussPoints, _Unknown);
    }

    void getLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q ) const  override {
        /* dof ordering: u1 v1 w1 p1  u2 v2 w2 p2  u3 v3 w3 p3  u4 v4 w4   u5 v5 w5  u6 v6 w6*/
        if (q == Variable::VariableQuantity::Displacement) {
            answer={1,2,3, 5,6,7, 9,10,11, 13,14,15, 17,18,19, 21,22,23, 25,26,27, 29,30,31 };
        } else if (q == Variable::VariableQuantity::Pressure) {
            answer = {4, 8, 12, 16, 20, 24, 28, 32};
        }
    }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override { 
            answer = {1,2,3,11};
    }
    int giveNumberOfDofs() override { return 32; }
    const char *giveInputRecordName() const override {return "upbrick11";}
    
    double computeVolumeAround(GaussPoint *gp) override {
        double determinant = fabs( this->pInterpol.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
        double weight = gp->giveWeight();
        return determinant * weight ;
    }
    IntegrationRule *giveDefaultIntegrationRulePtr() override {
        return &ir;
    }
    Element_Geometry_Type giveGeometryType() const override {
        return EGT_hexa_1;
    }
    private:
        virtual int  giveNumberOfUDofs() const override {return 24;} 
        virtual int  giveNumberOfPDofs() const override {return 8;}
        virtual const Variable& getU() const override {return u;}
        virtual const Variable& getP() const override {return p;}
};

const FEInterpolation & UPBrick11::uInterpol = FEI3dHexaLin();
const FEInterpolation & UPBrick11::pInterpol = FEI3dHexaLin();
const Variable& UPBrick11::p = Variable(UPBrick11::pInterpol, Variable::VariableQuantity::Pressure, Variable::VariableType::scalar, 3, NULL, {11});
const Variable& UPBrick11::u = Variable(UPBrick11::uInterpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 3, NULL, {1,2,3});

#define _IFT_UPBrick11_Name "upbrick11"
REGISTER_Element(UPBrick11)



#define _IFT_UPMaterial_Name "upm"
#define _IFT_UPMaterial_E "e"
#define _IFT_UPMaterial_nu "nu"
#define _IFT_UPMaterial_k "k"
#define _IFT_UPMaterial_alpha "alpha"
#define _IFT_UPMaterial_c "c"

class UPMaterial : public Material {
    protected:
        double e, nu; // elastic isotropic constants
        double k; // isotropic permeability
        double alpha; // Biot constant = 1-K_t/K_s (Kt bulk moduli of the porous medium, Ks bulk moduli of solid phase)
        double c; // 1/Q, where Q is combined compressibility of the fluid and solid phases (1/Q=n/Kt+(b-n)/Ks, where n is porosity) 
    public:
    UPMaterial (int n, Domain* d) : Material (n,d) {e=1.0; nu=0.15; k=1.0; alpha=1.0; c=0.1;}

    void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, GaussPoint* gp, TimeStep *tStep) override {
        if (type == StiffnessMatrix) {
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
        } else if (type == PermeabilityMatrix) {
            answer.resize(3,3);
            answer.beUnitMatrix();
            answer.times(this->k);
        }
    }

    void giveCharacteristicVector(FloatArray &answer, FloatArray& flux, CharType type, GaussPoint* gp, TimeStep *tStep) override {
        if (type == InternalForcesVector) {
            FloatMatrix d;
            this->giveCharacteristicMatrix(d, StiffnessMatrix, gp, tStep);
            answer.beProductOf(d, flux);
        }else if (type == FluidMassBalancePressureContribution) {
            FloatMatrix k;
            this->giveCharacteristicMatrix(k, PermeabilityMatrix, gp, tStep);
            answer.beProductOf(k, flux);
        }
    }


    double giveCharacteristicValue(CharType type, GaussPoint* gp, TimeStep *tStep) override {
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

        IR_GIVE_OPTIONAL_FIELD(ir, e, _IFT_UPMaterial_E);
        IR_GIVE_OPTIONAL_FIELD(ir, nu, _IFT_UPMaterial_nu);
        IR_GIVE_OPTIONAL_FIELD(ir, k, _IFT_UPMaterial_k);
        IR_GIVE_OPTIONAL_FIELD(ir, alpha, _IFT_UPMaterial_alpha);
        IR_GIVE_OPTIONAL_FIELD(ir, c, _IFT_UPMaterial_c);

    };
    //   void giveInputRecord(DynamicInputRecord &input) override {};

    const char *giveClassName() const override {return "UPMaterial";}
    const char *giveInputRecordName() const override {return _IFT_UPMaterial_Name;}

};
REGISTER_Material(UPMaterial)

} // end namespace oofem
