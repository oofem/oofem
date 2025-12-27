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
#include "termlibrary2.h"
#include "element.h"
#include "material.h"
#include "crosssection.h"

namespace oofem {


TMBTSigTerm::TMBTSigTerm (const Variable *testField, const Variable* unknownField, const Variable* temperatureField) : BTSigTerm(testField, unknownField), temperatureField(temperatureField) {}

void TMBTSigTerm::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    FloatArray eps, sig;
    FloatMatrix B;
    this->computeTMgeneralizedStrain(eps, B, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep);

    
    cell.giveCrossSection()->giveMaterial(gp)->giveCharacteristicVector(sig, eps, MatResponseMode::Stress, gp, tstep);
    answer.beTProductOf(B, sig);
}

void TMBTSigTerm::computeTMgeneralizedStrain (FloatArray& answer, FloatMatrix& B, MPElement& cell, const FloatArray& lcoords, MaterialMode mmode, TimeStep* tstep) const {
    FloatArray u, gradT;
    FloatMatrix dndx ;
    cell.getUnknownVector(u, this->field, VM_TotalIntrinsic, tstep);
    this->grad(B, this->field, this->field->interpolation, cell, lcoords, mmode);
    FloatArray Bu;
    Bu.beProductOf(B, u);

    FloatArray rt, Nt;
    cell.getUnknownVector(rt, temperatureField, VM_TotalIntrinsic, tstep);
    // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
    this->temperatureField->interpolation->evaldNdx(dndx, lcoords, FEIElementGeometryWrapper(&cell));
    // evaluate temperature gradient at given point
    gradT.beTProductOf(dndx, rt);
    // evaluate temperature at given point
    this->temperatureField->interpolation->evalN(Nt, lcoords, FEIElementGeometryWrapper(&cell));
    double t = Nt.dotProduct(rt);
    answer=FloatArray::fromConcatenated({Bu,gradT,Vec1(t)});
}

TMgNTfTerm::TMgNTfTerm (const Variable *testField, const Variable* unknownField, MatResponseMode lhsType, MatResponseMode rhsType) : gNTfTerm(testField, unknownField, lhsType, rhsType) {}
void TMgNTfTerm::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const {
    FloatArray sv(10), Nt, p, gradp, fp;
    FloatMatrix B;
    cell.getUnknownVector(p, this->field, VM_TotalIntrinsic, tstep);
    this->grad(B, this->field, this->field->interpolation, cell, gp->giveNaturalCoordinates());
    this->field->interpolation->evalN(Nt, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&cell));
    double t = Nt.dotProduct(p);
    gradp.beProductOf(B, p);
    sv(6) = gradp(0);
    sv(7) = gradp(1);
    sv(8) = gradp(2);
    sv(9) = t;
    cell.giveCrossSection()->giveMaterial(gp)->giveCharacteristicVector(fp, sv, rhsType, gp, tstep); // update
    answer.beTProductOf(B, fp);
}



BDalphaPiTerm::BDalphaPiTerm (const Variable *testField, const Variable* unknownField, ValueModeType m) : Term(testField, unknownField), m(m) {}


void BDalphaPiTerm::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const {
    FloatMatrix D, alphaPi, B, DaPI, BDaPI;
    FloatArray Nt;
    // alphaPi term
    e.giveCrossSection()->giveMaterial(gp)->giveCharacteristicMatrix(D, Conductivity, gp, tstep);  // 3x3 in 3D
    // expand it 
    alphaPi.resize(6,1);
    alphaPi.zero();
    alphaPi(0,0) = -D(0,0);
    alphaPi(1,0) = -D(1,1);
    alphaPi(2,0) = -D(2,2);
    e.giveCrossSection()->giveMaterial(gp)->giveCharacteristicMatrix(D, TangentStiffness, gp, tstep);  // 3x3 in 3D
    evalB(B, this->testField, this->testField->interpolation, e, gp->giveNaturalCoordinates(), gp->giveMaterialMode());
    this->field->interpolation->evalN(Nt, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));

    DaPI.beProductOf(D, alphaPi);
    BDaPI.beTProductOf(B,DaPI);
    FloatMatrix Ntm=FloatMatrix::fromArray(Nt, true);
    answer.beProductOf(BDaPI,Ntm);

}

void BDalphaPiTerm::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    // this is a partial linearization of BtSigma term with respect to temperature
    // thus the residual (rhs) contribution should come from BtSigma term
    answer.resize(this->testField->interpolation->giveNumberOfNodes(cell.giveGeometryType())*this->testField->size);
    answer.zero();
}

void BDalphaPiTerm::getDimensions(Element& cell) const  {}
void BDalphaPiTerm::initializeCell(Element& cell) const  {}


BTdSigmadT::BTdSigmadT (const Variable *testField, const Variable* unknownField) : Term(testField, unknownField) {}


void BTdSigmadT::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const {
    FloatMatrix D, B, DB;
    FloatArray Nt;
    // aplhaPi term
    e.giveCrossSection()->giveMaterial(gp)->giveCharacteristicMatrix(D, DSigmaDT, gp, tstep);
    this->field->interpolation->evalN(Nt, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    evalB(B, this->testField, this->testField->interpolation, e, gp->giveNaturalCoordinates(), gp->giveMaterialMode());
    FloatMatrix Ntm=FloatMatrix::fromArray(Nt, true);
    DB.beProductOf(D, Ntm);
    //answer.plusProductSymmUpper(B, DB, 1.0);
    answer.beTProductOf(B,DB);
}

void BTdSigmadT::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    // this is a partial linearization of BtSigma term with respect to temperature
    // thus the residual (rhs) contribution should come from BtSigma term
    answer.resize(this->testField->interpolation->giveNumberOfNodes(cell.giveGeometryType())*this->testField->size);
    answer.zero();
}

void BTdSigmadT::getDimensions(Element& cell) const  {}
void BTdSigmadT::initializeCell(Element& cell) const  {}


NTaTmTe:: NTaTmTe (const Variable *testField, const Variable* unknownField, BoundaryLoad* _bl, int bid, char btype) : Term(testField, unknownField), bl(_bl), boundaryID(bid), boundaryType(btype) {};
void NTaTmTe::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const {
    FloatArray Nt;
    if (boundaryType == 's') {
        this->testField->interpolation->boundarySurfaceEvalN(Nt, boundaryID, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    } else {
        this->testField->interpolation->boundaryEdgeEvalN(Nt, boundaryID, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    }
    answer.resize(Nt.giveSize(), Nt.giveSize());
    answer.zero();
    answer.plusDyadUnsym(Nt, Nt, this->bl->giveProperty('a', tstep));
}

void NTaTmTe::evaluate (FloatArray& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const  {
    FloatArray Nt, rt, Te, coords;
    if (boundaryType == 's') {
        this->testField->interpolation->boundarySurfaceEvalN(Nt, boundaryID, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    } else {
        this->testField->interpolation->boundaryEdgeEvalN(Nt, boundaryID, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    }
    // get surface unknown vector
    
    e.getBoundaryUnknownVector(rt, this->field, VM_TotalIntrinsic, this->boundaryID, this->boundaryType, tstep);
    double t = Nt.dotProduct(rt);
    answer= Nt;
    if ( this->bl->giveFormulationType() == Load :: FT_Entity ) {
        coords = gp->giveNaturalCoordinates();
    } else {
        //this->computeSurfIpGlobalCoords(gcoords, gp->giveNaturalCoordinates(), iSurf);
        e.getGeometryInterpolation()->boundarySurfaceLocal2global(coords, this->boundaryID, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    }
    this->bl->computeValues(Te, tstep, coords, this->field->dofIDs, VM_TotalIntrinsic);
    answer *= this->bl->giveProperty('a', tstep)*(t-Te.at(1));
}

void NTaTmTe::getDimensions(Element& cell) const  {}
void NTaTmTe::initializeCell(Element& cell) const  {}

InternalTMFluxSourceTerm::InternalTMFluxSourceTerm (const Variable *testField, const Variable* unknownField, const Variable* temperatureField) : TMBTSigTerm(testField, unknownField, temperatureField) {}

void InternalTMFluxSourceTerm::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    FloatArray eps, n, f;
    FloatMatrix B, N;
    this->computeTMgeneralizedStrain(eps, B, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep);
    cell.giveCrossSection()->giveMaterial(gp)->giveCharacteristicVector(f, eps, MatResponseMode::IntSource, gp, tstep);
    this->testField->interpolation->evalN(n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&cell));    
    N.beNMatrixOf(n, testField->size);
    answer.beTProductOf(N, f);
}


} // end namespace oofem
