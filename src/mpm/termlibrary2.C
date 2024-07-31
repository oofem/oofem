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
#include "termlibrary2.h"
#include "element.h"
#include "material.h"
#include "crosssection.h"

namespace oofem {


void deltaB(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords, const MaterialMode mmode) {
    FloatMatrix dndx;
    int nnodes = interpol.giveNumberOfNodes(cell.giveGeometryType());
    int ndofs = v.size;
    // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
    interpol.evaldNdx(dndx, coords, FEIElementGeometryWrapper(&cell));
    answer.resize(1, nnodes*ndofs);
    answer.zero();


    if (mmode == _3dUPV) {
        // 3D mode only now
        for (int i = 0; i< nnodes; i++) {
            answer(0, i*ndofs+0) = dndx(i, 0);
            answer(0, i*ndofs+1) = dndx(i, 1);
            answer(0, i*ndofs+2) = dndx(i, 2);
        }
    } else if (mmode == _2dUPV) {
        for (int i = 0; i< nnodes; i++) {
            answer(0, i*ndofs+0) = dndx(i, 0);
            answer(0, i*ndofs+1) = dndx(i, 1);
        }
    }
}

void evalB(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords, const MaterialMode mmode) {
    FloatMatrix dndx;
    int nnodes = interpol.giveNumberOfNodes(cell.giveGeometryType());
    int ndofs = v.size;
    // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
    interpol.evaldNdx(dndx, coords, FEIElementGeometryWrapper(&cell));
    answer.resize(6, nnodes*ndofs);
    answer.zero();
    if ((mmode == _3dUPV)||(mmode == _3dMat)) {
        // 3D mode 
        for (int i = 0; i< nnodes; i++) {
            answer(0, i*ndofs+0) = dndx(i, 0); // e_11
            answer(1, i*ndofs+1) = dndx(i, 1); // e_22 
            answer(2, i*ndofs+2) = dndx(i, 2); // e_33 
            answer(3, i*ndofs+1) = dndx(i, 2); // e_23
            answer(3, i*ndofs+2) = dndx(i, 1);
            answer(4, i*ndofs+0) = dndx(i, 2); // e_13
            answer(4, i*ndofs+2) = dndx(i, 0);
            answer(5, i*ndofs+0) = dndx(i, 1); // e_12
            answer(5, i*ndofs+1) = dndx(i, 0);
        }
    } else if (mmode == _2dUPV) {
        for (int i = 0; i< nnodes; i++) {
            answer(0, i*ndofs+0) = dndx(i, 0); // e_11
            answer(1, i*ndofs+1) = dndx(i, 1); // e_22
            answer(5, i*ndofs+0) = dndx(i, 1); // e_12
            answer(5, i*ndofs+1) = dndx(i, 0);
        }
    }
}

double evalVolumeFraction(const Variable&vf, MPElement& e, const FloatArray& coords, TimeStep* tstep)
{
    FloatArray rvf, Nvf;
    e.getUnknownVector(rvf, vf, VM_TotalIntrinsic, tstep);
    vf.interpolation.evalN(Nvf, coords, FEIElementGeometryWrapper(&e));
    return Nvf.dotProduct(rvf);
}

NTBdivTerm::NTBdivTerm (const Variable &testField, const Variable& unknownField, ValueModeType m) : Term(testField, unknownField), m(m) {}


void NTBdivTerm::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const {
    FloatArray Np;
    this->testField.interpolation.evalN(Np, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    FloatMatrix Npm(Np), B;   
    deltaB(B, this->field, this->field.interpolation, e, gp->giveNaturalCoordinates(), gp->giveMaterialMode());
    answer.beTProductOf(Npm,B);
}

void NTBdivTerm::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    FloatMatrix A;
    FloatArray u;
    this->evaluate_lin(A, cell, gp, tstep);
    cell.getUnknownVector(u, this->field, this->m, tstep);
    answer.beProductOf(A, u);
}

void NTBdivTerm::getDimensions(Element& cell) const  {
    //int nnodes = interpol.giveNumberOfNodes();
    //int ndofs = v.size;
    //return nnodes*ndofs;
}
void NTBdivTerm::initializeCell(Element& cell) const  {}


// deltaBTfiNpTerm 

deltaBTfiNpTerm::deltaBTfiNpTerm (const Variable &testField, const Variable& unknownField, const Variable& volumeFraction) : Term(testField, unknownField), volumeFraction(volumeFraction) {}


void deltaBTfiNpTerm::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const {
    FloatArray Nvf, Np;
    FloatMatrix B;   
    deltaB(B, this->testField, this->testField.interpolation, e, gp->giveNaturalCoordinates(), gp->giveMaterialMode());
    double vf = evalVolumeFraction(this->volumeFraction, e, gp->giveNaturalCoordinates(), tstep);
    this->field.interpolation.evalN(Np, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    FloatMatrix Npm(Np);
    answer.beTProductOf(B,Npm);
    answer.times(vf);
}

void deltaBTfiNpTerm::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    FloatMatrix A;
    FloatArray p;
    this->evaluate_lin(A, cell, gp, tstep);
    cell.getUnknownVector(p, this->field, VM_TotalIntrinsic, tstep);
    answer.beProductOf(A, p);
}

void deltaBTfiNpTerm::getDimensions(Element& cell) const  {
    
}
void deltaBTfiNpTerm::initializeCell(Element& cell) const  {}


// NdTdvfNpTerm 

NdTdvfNpTerm::NdTdvfNpTerm (const Variable &testField, const Variable& unknownField, const Variable& volumeFraction) : Term(testField, unknownField), volumeFraction(volumeFraction) {}


void NdTdvfNpTerm::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const {
    FloatArray nvec,dvf,rvf,Np;
    FloatMatrix N,dndx,help;
    this->testField.interpolation.evalN(nvec, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    N.beNMatrixOf(nvec, testField.size);
    // evaluate derivatives of volume fraction
    this->volumeFraction.interpolation.evaldNdx(dndx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    e.getUnknownVector(rvf, this->volumeFraction, VM_TotalIntrinsic, tstep);
    dvf.beTProductOf(dndx, rvf);
    this->field.interpolation.evalN(Np, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    FloatMatrix Npm(Np);
    help.beTProductOf(N,dvf); // Nv * dv
    answer.beProductOf(help, Npm);
}

void NdTdvfNpTerm::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    FloatMatrix A;
    FloatArray p;
    this->evaluate_lin(A, cell, gp, tstep);
    cell.getUnknownVector(p, this->field, VM_TotalIntrinsic, tstep);
    answer.beProductOf(A, p);
}

void NdTdvfNpTerm::getDimensions(Element& cell) const  {
    
}
void NdTdvfNpTerm::initializeCell(Element& cell) const  {}

// BTmuBTerm 

BTmuBTerm::BTmuBTerm (const Variable &testField, const Variable& unknownField) : Term(testField, unknownField) {}


void BTmuBTerm::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const {
    FloatArray nvec,dvf,rvf,Np;
    FloatMatrix B, M(6,6), MB;

    evalB(B, this->field, this->field.interpolation, e,gp->giveNaturalCoordinates(),  gp->giveMaterialMode());
    double m = e.giveCrossSection()->giveMaterial(gp)->giveCharacteristicValue(FluidViscosity, gp, tstep);
    M(0,0) = M(1,1) = M(2,2) = 2*m;
    M(3,3) = M(4,4) = M(5,5) = m;
    MB.beProductOf(M,B);
    answer.beTProductOf(B, MB);
}

void BTmuBTerm::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    FloatMatrix A;
    FloatArray r;
    this->evaluate_lin(A, cell, gp, tstep);
    cell.getUnknownVector(r, this->field, VM_TotalIntrinsic, tstep);
    answer.beTProductOf(A, r);
}

void BTmuBTerm::getDimensions(Element& cell) const  {
    
}
void BTmuBTerm::initializeCell(Element& cell) const  {}

// BTmuVfBTerm 

BTmuVfBTerm::BTmuVfBTerm (const Variable &testField, const Variable& unknownField, const Variable& volumeFraction) : Term(testField, unknownField), volumeFraction(volumeFraction) {}


void BTmuVfBTerm::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const {
    FloatMatrix B, M(6,6), MB;

    evalB(B, this->field, this->field.interpolation, e,gp->giveNaturalCoordinates(),  gp->giveMaterialMode());
    double m = e.giveCrossSection()->giveMaterial(gp)->giveCharacteristicValue(FluidViscosity, gp, tstep);
    double vf = evalVolumeFraction(this->volumeFraction, e, gp->giveNaturalCoordinates(), tstep);

    M(0,0) = M(1,1) = M(2,2) = 2*m*vf;
    M(3,3) = M(4,4) = M(5,5) = m*vf;
    MB.beProductOf(M,B);
    answer.beTProductOf(B, MB);
}

void BTmuVfBTerm::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    FloatMatrix A;
    FloatArray r;
    this->evaluate_lin(A, cell, gp, tstep);
    cell.getUnknownVector(r, this->field, VM_Velocity, tstep);
    answer.beTProductOf(A, r);
}

void BTmuVfBTerm::getDimensions(Element& cell) const  {
    
}
void BTmuVfBTerm::initializeCell(Element& cell) const  {}

// NTmuVfSNTerm 

NTmuVfSNTerm::NTmuVfSNTerm (const Variable &testField, const Variable& unknownField, const Variable& volumeFraction) : Term(testField, unknownField), volumeFraction(volumeFraction) {}


void NTmuVfSNTerm::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const {
    
  FloatMatrix S, SI, SIN, Nd;
    FloatArray N;

    this->field.interpolation.evalN(N, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    Nd.beNMatrixOf (N, testField.size);
    e.giveCrossSection()->giveMaterial(gp)->giveCharacteristicMatrix(S, Permeability,  gp, tstep);
    double m = e.giveCrossSection()->giveMaterial(gp)->giveCharacteristicValue(FluidViscosity, gp, tstep);
    double vf = evalVolumeFraction(this->volumeFraction, e, gp->giveNaturalCoordinates(), tstep);
    SI.beInverseOf(S);
    SI.times(m*vf);
    
    SIN.beProductOf(SI,Nd);    
    answer.beTProductOf(Nd, SIN); 
}

void NTmuVfSNTerm::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    FloatMatrix A;
    FloatArray r;
    this->evaluate_lin(A, cell, gp, tstep);
    cell.getUnknownVector(r, this->field, VM_TotalIntrinsic, tstep);
    answer.beTProductOf(A, r);
}

void NTmuVfSNTerm::getDimensions(Element& cell) const  {
    
}
void NTmuVfSNTerm::initializeCell(Element& cell) const  {}



// deltaBTNpTerm 

deltaBTNpTerm::deltaBTNpTerm (const Variable &testField, const Variable& unknownField) : Term(testField, unknownField) {}


void deltaBTNpTerm::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const {
    FloatArray Nvf, Np;
    FloatMatrix B;   
    deltaB(B, this->testField, this->testField.interpolation, e, gp->giveNaturalCoordinates(), gp->giveMaterialMode());
    this->field.interpolation.evalN(Np, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    FloatMatrix Npm(Np);
    answer.beTProductOf(B,Npm);
}

void deltaBTNpTerm::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    FloatMatrix A;
    FloatArray p;
    this->evaluate_lin(A, cell, gp, tstep);
    cell.getUnknownVector(p, this->field, VM_TotalIntrinsic, tstep);
    answer.beProductOf(A, p);
}

void deltaBTNpTerm::getDimensions(Element& cell) const  {
    
}
void deltaBTNpTerm::initializeCell(Element& cell) const  {}

} // end namespace oofem
