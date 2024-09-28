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
#include "termlibrary4.h"
#include "element.h"
#include "mathfem.h"

namespace oofem {

NTN::NTN (const Variable &testField, const Variable& unknownField) : Term(testField, unknownField) {}

void NTN::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const  {
    FloatArray Np;
    this->field.interpolation.evalN(Np, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    answer.beDyadicProductOf(Np, Np);
}

void NTN::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    FloatArray p;
    FloatMatrix S;
    cell.getUnknownVector(p, this->field, VM_Velocity, tstep);
    this->evaluate_lin (S, cell,gp, tstep);
    answer.beProductOf(S,p);
}

void NTN::getDimensions(Element& cell) const  {
    //int nnodes = interpol.giveNumberOfNodes();
    //int ndofs = v.size;
    //return nnodes*ndofs;
}
void NTN::initializeCell(Element& cell) const  {}


dnTaN::dnTaN (const Variable &testField, const Variable& unknownField) : Term(testField, unknownField) {}

void dnTaN::evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const  {
    FloatArray n;
    FloatMatrix dndx;
    this->testField.interpolation.evaldNdx(dndx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    this->field.interpolation.evalN(n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(&e));
    FloatMatrix nm(n, true), anm;
    FloatMatrix a({sqrt(0.5),sqrt(0.5)}, false); // assume 2D for now
    anm.beProductOf(a, nm);
    answer.beProductOf(dndx, anm);
}

void dnTaN::evaluate (FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  {
    FloatArray p;
    FloatMatrix S;
    cell.getUnknownVector(p, this->field, VM_Velocity, tstep);
    this->evaluate_lin (S, cell,gp, tstep);
    answer.beProductOf(S,p);
}

void dnTaN::getDimensions(Element& cell) const  {
    //int nnodes = interpol.giveNumberOfNodes();
    //int ndofs = v.size;
    //return nnodes*ndofs;
}
void dnTaN::initializeCell(Element& cell) const  {}







} // end namespace oofem
