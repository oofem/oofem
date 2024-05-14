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
#include "crosssection.h"

namespace oofem {

BTSigTerm::BTSigTerm (const Variable &testField, const Variable& unknownField) : Term(testField, unknownField) {}


void BTSigTerm::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const {
    FloatMatrix D, B, DB;
    e.giveCrossSection()->giveMaterial(gp)->giveCharacteristicMatrix(D, StiffnessMatrix, gp, tstep);
    this->grad(B, this->field, this->field.interpolation, e, gp->giveNaturalCoordinates(), gp->giveMaterialMode());
    DB.beProductOf(D, B);
    //answer.plusProductSymmUpper(B, DB, 1.0);
    answer.beTProductOf(B,DB);
}

void BTSigTerm::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    FloatArray u, eps, sig;
    FloatMatrix B;
    cell.getUnknownVector(u, this->field, VM_TotalIntrinsic, tstep);
    this->grad(B, this->field, this->field.interpolation, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode());
    eps.beProductOf(B, u);
    cell.giveCrossSection()->giveMaterial(gp)->giveCharacteristicVector(sig, eps, InternalForcesVector, gp, tstep);
    answer.beTProductOf(B, sig);
}

void BTSigTerm::getDimensions(Element& cell) const  {
    //int nnodes = interpol.giveNumberOfNodes();
    //int ndofs = v.size;
    //return nnodes*ndofs;
}
void BTSigTerm::initializeCell(Element& cell) const  {}

void BTSigTerm::grad(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords, const MaterialMode mmode) const {
    FloatMatrix dndx;
    int nnodes = interpol.giveNumberOfNodes();
    int ndofs = v.size;
    // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
    interpol.evaldNdx(dndx, coords, FEIElementGeometryWrapper(&cell));

    if ((mmode == _3dUP) || (mmode == _3dUPV)) {
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
    }
}

//wTgNTfTerm class (H)


gNTfTerm::gNTfTerm (const Variable &testField, const Variable& unknownField) : Term(testField, unknownField) {}


void gNTfTerm::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const  {
    FloatMatrix D, B, DB;
    e.giveCrossSection()->giveMaterial(gp)->giveCharacteristicMatrix(D, PermeabilityMatrix, gp, tstep); // update
    this->grad(B, this->field, this->field.interpolation, e, gp->giveNaturalCoordinates());
    DB.beProductOf(D, B);
    answer.beTProductOf(B, DB);
}

void gNTfTerm::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    FloatArray p, gradp, fp;
    FloatMatrix B;
    cell.getUnknownVector(p, this->field, VM_TotalIntrinsic, tstep);
    this->grad(B, this->field, this->field.interpolation, cell, gp->giveNaturalCoordinates());
    gradp.beProductOf(B, p);
    cell.giveCrossSection()->giveMaterial(gp)->giveCharacteristicVector(fp, gradp, FluidMassBalancePressureContribution, gp, tstep); // update
    answer.beTProductOf(B, fp);
}

void gNTfTerm::getDimensions(Element& cell) const  {
    //int nnodes = interpol.giveNumberOfNodes();
    //int ndofs = v.size;
    //return nnodes*ndofs;
}
void gNTfTerm::initializeCell(Element& cell) const  {}

void gNTfTerm::grad(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords) const {
    FloatMatrix at;
    // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
    interpol.evaldNdx(at, coords, FEIElementGeometryWrapper(&cell));
    answer.beTranspositionOf(at);
}

// BTamN Term (Qp)

BTamNTerm::BTamNTerm (const Variable &testField, const Variable& unknownField) : Term(testField, unknownField) {}

void BTamNTerm::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const  {
    FloatMatrix B, mn;
    FloatArray m({1,1,1,0,0,0}), Np;
    this->field.interpolation.evalN(Np, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    m.times(e.giveCrossSection()->giveMaterial(gp)->giveCharacteristicValue(BiotConstant, gp, tstep));
    mn.beDyadicProductOf(m, Np);
    this->grad(B, this->testField, this->testField.interpolation, e, gp->giveNaturalCoordinates(), gp->giveMaterialMode());
    answer.beTProductOf(B, mn);
}

void BTamNTerm::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    FloatArray p;
    FloatMatrix Q;
    cell.getUnknownVector(p, this->field, VM_TotalIntrinsic, tstep);
    this->evaluate_lin (Q, cell,gp, tstep);
    answer.beProductOf(Q,p);
}

void BTamNTerm::getDimensions(Element& cell) const  {
    //int nnodes = interpol.giveNumberOfNodes();
    //int ndofs = v.size;
    //return nnodes*ndofs;
}
void BTamNTerm::initializeCell(Element& cell) const  {}

void BTamNTerm::grad(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords, const MaterialMode mmode) const {
 FloatMatrix dndx;
    int nnodes = interpol.giveNumberOfNodes();
    int ndofs = v.size;
    // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
    interpol.evaldNdx(dndx, coords, FEIElementGeometryWrapper(&cell));

    if (mmode == _3dUP) {
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
    } else if (mmode == _2dUP) {
        answer.resize(6, nnodes*ndofs);
        for (int i = 0; i< nnodes; i++) {
            answer(0, i*ndofs+0) = dndx(i, 0);
            answer(1, i*ndofs+1) = dndx(i, 1);

            answer(5, i*ndofs+0) = dndx(i, 1);
            answer(5, i*ndofs+1) = dndx(i, 0);
        }
    }
}

// NTamTBTerm Term (Q^T du/dt)

NTamTBTerm::NTamTBTerm (const Variable &testField, const Variable& unknownField) : Term(testField, unknownField) {}

void NTamTBTerm::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const  {
    FloatMatrix B, mb;
    FloatArray m({1,1,1,0,0,0}), Np;
    this->testField.interpolation.evalN(Np, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    m.times(e.giveCrossSection()->giveMaterial(gp)->giveCharacteristicValue(BiotConstant, gp, tstep));
    this->grad(B, this->field, this->field.interpolation, e, gp->giveNaturalCoordinates(), gp->giveMaterialMode());
    mb.beTProductOf(m, B);
    FloatMatrix Npm(Np);
    answer.beProductOf(Npm, mb);
}

void NTamTBTerm::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    FloatArray ut;
    FloatMatrix Q;
    cell.getUnknownVector(ut, this->field, VM_Velocity ,tstep);
    this->evaluate_lin (Q, cell,gp, tstep);
    answer.beProductOf(Q,ut);
}

void NTamTBTerm::getDimensions(Element& cell) const  {
    //int nnodes = interpol.giveNumberOfNodes();
    //int ndofs = v.size;
    //return nnodes*ndofs;
}
void NTamTBTerm::initializeCell(Element& cell) const  {}

void NTamTBTerm::grad(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords, const MaterialMode mmode) const {
 FloatMatrix dndx;
    int nnodes = interpol.giveNumberOfNodes();
    int ndofs = v.size;
    // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
    interpol.evaldNdx(dndx, coords, FEIElementGeometryWrapper(&cell));

    if (mmode == _3dUP) {
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
    } else if (mmode == _2dUP) {
        answer.resize(6, nnodes*ndofs);
        for (int i = 0; i< nnodes; i++) {
            answer(0, i*ndofs+0) = dndx(i, 0);
            answer(1, i*ndofs+1) = dndx(i, 1);

            answer(5, i*ndofs+0) = dndx(i, 1);
            answer(5, i*ndofs+1) = dndx(i, 0);
        }
    }
}


// NTcN Term (S(dp/dt))

NTcN::NTcN (const Variable &testField, const Variable& unknownField) : Term(testField, unknownField) {}

void NTcN::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const  {
    FloatArray Np;
    this->field.interpolation.evalN(Np, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    answer.beDyadicProductOf(Np, Np);
    answer.times(e.giveCrossSection()->giveMaterial(gp)->giveCharacteristicValue(CompressibilityCoefficient, gp, tstep));
}

void NTcN::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    FloatArray p;
    FloatMatrix S;
    cell.getUnknownVector(p, this->field, VM_Velocity, tstep);
    this->evaluate_lin (S, cell,gp, tstep);
    answer.beProductOf(S,p);
}

void NTcN::getDimensions(Element& cell) const  {
    //int nnodes = interpol.giveNumberOfNodes();
    //int ndofs = v.size;
    //return nnodes*ndofs;
}
void NTcN::initializeCell(Element& cell) const  {}

// A external flux term $S=(N)^T f$, where $f$ is functor evaluating the flux. 
NTf_Surface::NTf_Surface (const Variable &testField, const NTfFunctor& f, int surf) : Term (testField, testField), f(f), isurf(surf) {}

void 
NTf_Surface::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tStep) const {
  FloatArray nvec, flux;
  FloatMatrix N;
  const FloatArray& lc = gp->giveNaturalCoordinates();

  this->f.evaluate(flux, lc, cell, this->testField, tStep);
  
  this->testField.interpolation.boundarySurfaceEvalN(nvec, isurf, lc, FEIElementGeometryWrapper(&cell));
  N.beNMatrixOf(nvec, testField.size);
  answer.beTProductOf(N, flux);
  
}

// A external flux term $S=(N)^T f$, where $f$ is functor evaluating the flux. 
NTf_Edge::NTf_Edge (const Variable &testField, const NTfFunctor& f, int surf) : Term (testField, testField), f(f), isurf(surf) {}

void 
NTf_Edge::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tStep) const {
  FloatArray nvec, flux;
  FloatMatrix N;
  const FloatArray& lc = gp->giveNaturalCoordinates();

  this->f.evaluate(flux, lc, cell, this->testField, tStep);
  
  this->testField.interpolation.boundaryEdgeEvalN(nvec, isurf, lc, FEIElementGeometryWrapper(&cell));
  N.beNMatrixOf(nvec, testField.size);
  answer.beTProductOf(N, flux);
  
}


} // end namespace oofem
