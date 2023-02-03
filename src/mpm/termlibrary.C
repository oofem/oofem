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
#include "material.h"

namespace oofem {

BTSigTerm::BTSigTerm (Variable& unknownField, Variable &testField) : Term(unknownField, testField) {}


// assuming symmmetric form 
void BTSigTerm::evaluate_dw (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep)  {
    FloatMatrix D, B, DB;
    e.giveMaterial()->giveCharacteristicMatrix(D, StiffnessMatrix, gp, tstep);
    this->grad(B, this->field, this->field.interpolation, e, gp->giveNaturalCoordinates());
    DB.beProductOf(D, B);
    answer.plusProductSymmUpper(B, DB, 1.0);
}

void BTSigTerm::evaluate_c (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep)  {
    FloatArray u, eps, sig;
    FloatMatrix B;
    cell.getUnknownVector(u, this->field, tstep);
    this->grad(B, this->field, this->field.interpolation, cell, gp->giveNaturalCoordinates());
    eps.beProductOf(B, u);
    cell.giveMaterial()->giveCharacteristicVector(sig, eps, InertiaForcesVector, gp, tstep);
    answer.beTProductOf(B, eps);
}

void BTSigTerm::getDimensions_dw(Element& cell)  {
    //int nnodes = interpol.giveNumberOfNodes();
    //int ndofs = v.size;
    //return nnodes*ndofs;
}
void BTSigTerm::initializeCell(Element& cell)  {}

void BTSigTerm::grad(FloatMatrix& answer, Variable &v, FEInterpolation& interpol, Element& cell, const FloatArray& coords) {
    FloatMatrix dndx;
    int nnodes = interpol.giveNumberOfNodes();
    int ndofs = v.size;
    // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
    interpol.evaldNdx(dndx, coords, FEIElementGeometryWrapper(&cell));


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
}

//wTgNTfTerm class


gNTfTerm::gNTfTerm (Variable& unknownField, Variable &testField) : Term(unknownField, testField) {}


// assuming symmmetric form 
void gNTfTerm::evaluate_dw (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep)  {
    FloatMatrix D, B, DB;
    e.giveMaterial()->giveCharacteristicMatrix(D, PermeabilityMatrix, gp, tstep); // update
    this->grad(B, this->field, this->field.interpolation, e, gp->giveNaturalCoordinates());
    DB.beProductOf(D, B);
    answer.plusProductSymmUpper(B, DB, 1.0);
}

void gNTfTerm::evaluate_c (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep)  {
    FloatArray p, gradp, fp;
    FloatMatrix B;
    cell.getUnknownVector(p, this->field, tstep);
    this->grad(B, this->field, this->field.interpolation, cell, gp->giveNaturalCoordinates());
    gradp.beProductOf(B, p);
    cell.giveMaterial()->giveCharacteristicVector(fp, gradp, FluidMassBalancePressureContribution, gp, tstep); // update
    answer.beTProductOf(B, fp);
}

void gNTfTerm::getDimensions_dw(Element& cell)  {
    //int nnodes = interpol.giveNumberOfNodes();
    //int ndofs = v.size;
    //return nnodes*ndofs;
}
void gNTfTerm::initializeCell(Element& cell)  {}

void gNTfTerm::grad(FloatMatrix& answer, Variable &v, FEInterpolation& interpol, Element& cell, const FloatArray& coords) {
    FloatMatrix at;
    // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
    interpol.evaldNdx(at, coords, FEIElementGeometryWrapper(&cell));
    answer.beTranspositionOf(at);
}

// BTamN Term

BTamNTerm::BTamNTerm (Variable& unknownField, Variable &testField) : Term(unknownField, testField) {}

void BTamNTerm::evaluate_dw (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep)  {
    FloatMatrix B, mn;
    FloatArray m({1,1,1,0,0,0}), Np;
    this->field.interpolation.evalN(Np, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    m.times(e.giveMaterial()->giveCharacteristicValue(BiotConstant, gp, tstep));
    mn.beDyadicProductOf(m, Np);
    this->grad(B, this->testField, this->testField.interpolation, e, gp->giveNaturalCoordinates());
    answer.beTProductOf(B, mn);
}

void BTamNTerm::evaluate_c (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep)  {
    FloatArray p;
    FloatMatrix Q;
    cell.getUnknownVector(p, this->field, tstep);
    this->evaluate_dw (Q, cell,gp, tstep);
    answer.beProductOf(Q,p);
}

void BTamNTerm::getDimensions_dw(Element& cell)  {
    //int nnodes = interpol.giveNumberOfNodes();
    //int ndofs = v.size;
    //return nnodes*ndofs;
}
void BTamNTerm::initializeCell(Element& cell)  {}

void BTamNTerm::grad(FloatMatrix& answer, Variable &v, FEInterpolation& interpol, Element& cell, const FloatArray& coords) {
 FloatMatrix dndx;
    int nnodes = interpol.giveNumberOfNodes();
    int ndofs = v.size;
    // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
    interpol.evaldNdx(dndx, coords, FEIElementGeometryWrapper(&cell));


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
}


// NTcN Term

NTcN::NTcN (Variable& unknownField, Variable &testField) : Term(unknownField, testField) {}

void NTcN::evaluate_dw (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep)  {
    FloatArray Np;
    this->field.interpolation.evalN(Np, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    answer.beDyadicProductOf(Np, Np);
    answer.times(e.giveMaterial()->giveCharacteristicValue(CompressibilityCoefficient, gp, tstep));
}

void NTcN::evaluate_c (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep)  {
    FloatArray p;
    FloatMatrix S;
    cell.getUnknownVector(p, this->field, tstep);
    this->evaluate_dw (S, cell,gp, tstep);
    answer.beProductOf(S,p);
}

void NTcN::getDimensions_dw(Element& cell)  {
    //int nnodes = interpol.giveNumberOfNodes();
    //int ndofs = v.size;
    //return nnodes*ndofs;
}
void NTcN::initializeCell(Element& cell)  {}



} // end namespace oofem