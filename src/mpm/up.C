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
#include "mathfem.h"

#include "material.h"

namespace oofem {

#define _IFT_UPElement_Name "up"


class UPElement : public MPElement {
    protected:
        FEI3dTetLin pInterpol;
        FEI3dTetQuad uInterpol;
        Variable p;
        Variable u;
        // Terms
        BTSigTerm tm;
        gNTfTerm th;
        BTamNTerm tq;
        NTamTBTerm tqt;
        NTcN ts;

        GaussIntegrationRule ir;
        
    public:
    UPElement(int n, Domain* d): 
        MPElement(n,d), 
        pInterpol(), uInterpol(),
        p(this->pInterpol, Variable::VariableQuantity::Pressure, Variable::VariableType::scalar, 3, NULL, {11}), 
        u(this->uInterpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 3, NULL, {1,2,3}),
        tm(u,u), th(p,p), tq(u,p), tqt(p,u), ts(p,p), 
        ir(1, this)
    {
        numberOfDofMans  = 10;
        numberOfGaussPoints = 4;
        ir.SetUpPointsOnTetrahedra(numberOfGaussPoints, _Unknown);
    }

    // Note: performance can be probably improved once it will be possible 
    // to directly assemble multiple term contributions to the system matrix.
    // template metaprogramming?
    void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep) override {

        if (type == MomentumBalance_StiffnessMatrix) {
            answer.resize(30,30);
            answer.zero();
            this->integrateTerm_dw (answer, this->tm, &this->ir, tStep) ;
        } else if (type == MomentumBalance_PressureCouplingMatrix) {
            answer.resize(30,4);
            answer.zero();
            this->integrateTerm_dw (answer, this->tq, &this->ir, tStep) ;
        } else if (type == MassBalance_PermeabilityMatrix) {
            answer.resize(4,4);
            answer.zero();
            this->integrateTerm_dw (answer, this->th, &this->ir, tStep) ;
        } else if (type == MassBalance_CompresibilityMatrix) {
            answer.resize(4,4);
            answer.zero();
            this->integrateTerm_dw (answer, this->ts, &this->ir, tStep) ;
        } else if (type == MassBalance_StressCouplingMatrix) {
            answer.resize(4,30);
            answer.zero();
            this->integrateTerm_dw (answer, this->tqt, &this->ir, tStep) ;
        } else {
	  OOFEM_ERROR("Unknown characteristic matrix type");
	}
    }

    void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) override {
        if (type == MomentumBalance_StressResidual) {
            answer.resize(30);
            answer.zero();
            this->integrateTerm_c (answer, this->tm, &this->ir, tStep) ;
        } else if (type == MomentumBalance_PressureResidual) {
            answer.resize(30);
            answer.zero();
            this->integrateTerm_c(answer, this->tq, &this->ir, tStep) ;
        } else if (type == MassBalance_StressRateResidual) {
            answer.resize(4);
            answer.zero();
            this->integrateTerm_c (answer, this->tqt, &this->ir, tStep) ;
        } else if (type == MassBalance_PressureResidual) {
            answer.resize(4);
            answer.zero();
            this->integrateTerm_c (answer, this->th, &this->ir, tStep) ;
        } else if (type == MassBalance_PressureRateResidual) {
            answer.resize(4);
            answer.zero();
            this->integrateTerm_c (answer, this->ts, &this->ir, tStep) ;   
        } else {
	  OOFEM_ERROR("Unknown characteristic vector type");
	}
    }

    void getLocalCodeNumbers (IntArray& answer, Variable::VariableQuantity q ) override {
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
    const char *giveInputRecordName() const override {return "up";}
    
    double computeVolumeAround(GaussPoint *gp) override {
        double determinant = fabs( this->pInterpol.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
        double weight = gp->giveWeight();
        return determinant * weight ;
    }
};

REGISTER_Element(UPElement)




#define _IFT_UPMaterial_Name "upm"
class UPMaterial : public Material {
    protected:
        double e, nu; // elastic isotropic constants
        double k; // isotropic permeability
        double b; // Biot constant
        double c; // Compressibility coefficient
    public:
    UPMaterial (int n, Domain* d) : Material (n,d) {e=1.0; nu=0.15; k=1.0; b=1.0; c=0.1;}

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
            return b;
        } else if (type == CompressibilityCoefficient) {
            return c;
        } else {
            return 0.0;
        }
    };
    void initializeFrom(InputRecord &ir) override {};
    void giveInputRecord(DynamicInputRecord &input) override {};

    const char *giveClassName() const override {return "UPMaterial";}
    const char *giveInputRecordName() const override {return _IFT_UPMaterial_Name;}

};
REGISTER_Material(UPMaterial)

} // end namespace oofem
