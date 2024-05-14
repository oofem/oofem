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
#include "element.h"
#include "gausspoint.h"
#include "feinterpol.h"
#include "intarray.h"
#include "classfactory.h"

// for demo
#include "fei2dtrlin.h"
#include "gaussintegrationrule.h"


namespace oofem {


class PoissonTerm : public Term {
    protected:
        double c;
    public:
    PoissonTerm (const Variable& unknownField, const Variable &testField, double c) : Term(unknownField, testField) {
        this->c = c;
    }

    void grad(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords) const {
        interpol.evaldNdx(answer, coords, FEIElementGeometryWrapper(&cell));
    }


    void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep *tstep) const override {
        const FEInterpolation & si = field.interpolation;
        const FEInterpolation &ti = testField.interpolation;
        FloatMatrix bs, bt;
        this->grad(bs, this->field,si,e,gp->giveNaturalCoordinates());
        this->grad(bt, this->testField,ti,e,gp->giveNaturalCoordinates());

        FloatMatrix gc, c(2,2);
        c.at(1,1) = this->c;
        c.at(2,2) = this->c;

        gc.beProductOf(bs, c);
        answer.beProductTOf(gc, bt);
    }

    void evaluate (FloatArray&, MPElement& cell, GaussPoint*gp, TimeStep* tstep) const override {}
    void getDimensions(Element& cell) const override {}
    void initializeCell(Element& cell) const override {}
};



#define _IFT_PoissonElement_Name "pe"


class PoissonElement : public MPElement {
     protected:
        FEI2dTrLin interpol;
        Variable t;
        Variable dt;
        PoissonTerm p;
        GaussIntegrationRule ir;
    public:
    PoissonElement(int n, Domain* d): 
        MPElement(n,d), 
        interpol(1,2), 
        t(interpol, Variable::VariableQuantity::Temperature, Variable::VariableType::scalar, 3), 
        dt(interpol, Variable::VariableQuantity::Temperature, Variable::VariableType::scalar, 3, &t),
        p(t,dt,1.0),
        ir(1, this)
    {
        numberOfDofMans  = 3;
        numberOfGaussPoints = 1;
        ir.SetUpPointsOnTriangle(numberOfGaussPoints, _2dHeat);
    }

    // Note: performance can be probably improved once it will be possible 
    // to directly assemble multiple term contributions to the system matrix.
    // template metaprogramming?
    void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep) override {
        if (type == ConductivityMatrix) {
            FloatMatrix term;
            answer.resize(3,3);
            this->integrateTerm_dw (term, this->p, &this->ir, tStep) ;
            this->assembleTermContribution(answer, term, this->p);
            answer.printYourself("Conductivity");

        }
    }

    void getDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int n) const override {
        answer = {n};
    }
    void getInternalDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {
        answer={};
    }

    int getNumberOfSurfaceDOFs() const override {return 0;}
    int getNumberOfEdgeDOFs() const override {return 0;}
    void getSurfaceLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {
        answer.clear();
    }
    void getEdgeLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {
        answer.clear();
    }
    
    const char *giveInputRecordName() const override {return "pe";}
    
    const FEInterpolation& getGeometryInterpolation() const override {return  this->interpol;}
};

REGISTER_Element(PoissonElement)


} // end namespace oofem
