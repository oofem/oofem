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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "linesurfacetension.h"
#include "node.h"
#include "domain.h"
#include "material.h"
#include "mathfem.h"
#include "fei2dlinelin.h"

namespace oofem {
FEI2dLineLin LineSurfaceTension :: interp(1, 2);

LineSurfaceTension :: LineSurfaceTension(int n, Domain *aDomain) : FMElement (n, aDomain)
{
    this->numberOfDofMans = 2;
    this->integrationRulesArray = 0;
}

LineSurfaceTension :: ~LineSurfaceTension () {}

IRResultType LineSurfaceTension :: initializeFrom(InputRecord *ir)
{
    return Element :: initializeFrom(ir);
}

FEInterpolation * LineSurfaceTension :: giveInterpolation()
{
    return &this->interp;
}

double LineSurfaceTension :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray c;
    c = *this->giveNode(1)->giveCoordinates();
    c.add(*this->giveNode(2)->giveCoordinates());
    c.times(0.5);
    c.subtract(coords);
    return c.computeNorm();
}

int LineSurfaceTension :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode,
        TimeStep *tStep, const FloatArray &gcoords, FloatArray &answer)
{
    FloatArray lcoords;
    if (!this->computeLocalCoordinates(lcoords, gcoords)) {
        answer.resize(0);
        return false;
    }
    this->EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(mode, tStep, lcoords, answer);
    return true;
}

void LineSurfaceTension :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
        TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    FloatArray n;
    this->interp.evalN(n, lcoords, FEIElementGeometryWrapper(this));

    answer.resize(2);
    answer.zero();
    for (int i = 1; i < n.giveSize(); i++) {
        answer(0) += n.at(i)*this->giveNode(i)->giveDofWithID(V_u)->giveUnknown(EID_MomentumBalance, mode, tStep);
        answer(1) += n.at(i)*this->giveNode(i)->giveDofWithID(V_v)->giveUnknown(EID_MomentumBalance, mode, tStep);
    }
}

void LineSurfaceTension :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
    answer.setValues(2, V_u, V_v);
}

void LineSurfaceTension :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
// Returns the mask for this element. Should i handle EquationID gracefully?
{
    if (ut == EID_MomentumBalance || ut == EID_MomentumBalance_ConservationEquation) {
        if (this->giveDomain()->giveDomainType() == _2dIncompressibleFlow) {
            answer.setValues(2, V_u, V_v);
        } else {
            answer.setValues(2, D_u, D_v);
        }
    } else {
        answer.resize(0);
    }
}

void LineSurfaceTension :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode, TimeStep *tStep)
{
    if( mtrx == ExternalForcesVector )  {
        this->computeLoadVector(answer, mode, tStep);
    } else if ( mtrx == InternalForcesVector ){
        answer.resize(0);
    } else {
        OOFEM_ERROR2("giveCharacteristicVector: Unknown CharType (%s).",__CharTypeToString(mtrx));
    }
}

void LineSurfaceTension :: giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep)
{
    if ( mtrx == StiffnessMatrix )  {
        this->computeTangent(answer, tStep);
    } else {
        OOFEM_ERROR2("giveCharacteristicMatrix: Unknown CharType (%s)",__CharTypeToString(mtrx));
    }
}

void LineSurfaceTension :: computeLoadVector(FloatArray &answer, ValueModeType mode, TimeStep *tStep)
{
    domainType dt = this->giveDomain()->giveDomainType();

    Node *node1, *node2;
    FloatArray v1, v2, v;
    double x1, x2, length, t, gamma_s;

    gamma_s = this->giveMaterial()->give('g',NULL);

    node1 = giveNode(1);
    node2 = giveNode(2);

    v.beDifferenceOf(*node2->giveCoordinates(), *node1->giveCoordinates());
    length = v.computeNorm();
    v.times(1.0/length);

    x1 = node1->giveCoordinate(1);
    x2 = node1->giveCoordinate(2);

    answer.resize(4);
    answer.zero();

    if (dt == _3dAxisymmMode) {
        OOFEM_CLASS_ERROR("Not tested");
        t = M_PI*(x1+x2);
        answer.at(1) = v.at(1)*t - length;
        answer.at(2) = v.at(2)*t;
        answer.at(3) = -v.at(1)*t - length;
        answer.at(4) = -v.at(2)*t;
    } else {
        //t = this->giveDomain()->giveCrossSection(1)->give(CS_Thickness);
        t = 1.0;
        // In this case simple linear case, the boundary term would cancel out the edge term,
        // so it can be simplified to check when the boundary is not there.
        answer.at(1) = v.at(1)*t;
        answer.at(2) = v.at(2)*t;
        answer.at(3) = -v.at(1)*t;
        answer.at(4) = -v.at(2)*t;
    }

    answer.times(gamma_s);
}

void LineSurfaceTension :: computeTangent(FloatMatrix &answer, TimeStep *tStep)
{
#if 1
    ///@todo Not sure if it's a good idea to use this tangent.
    answer.resize(4,4);
    answer.zero();
#else
    domainType dt = this->giveDomain()->giveDomainType();
    int ndofs = this->computeNumberOfDofs(EID_MomentumBalance);
    Node *node1, *node2;
    double x1, x2, y1, y2, dx, dy, vx, vy, length, width, gamma_s;

    gamma_s = this->giveMaterial()->give('g',NULL);

    node1 = giveNode(1);
    node2 = giveNode(2);

    x1 = node1->giveCoordinate(1);
    x2 = node2->giveCoordinate(1);

    y1 = node1->giveCoordinate(2);
    y2 = node2->giveCoordinate(2);

    dx = x2-x1;
    dy = y2-y1;

    length = sqrt(dx*dx + dy*dy);

    vx = dx/length;
    vy = dy/length;

    FloatArray Ah(4);
    Ah.at(1) = -vx;
    Ah.at(2) = -vy;
    Ah.at(3) = vx;
    Ah.at(4) = vy;

    FloatMatrix NpTNp(4,4);
    NpTNp.zero();
    NpTNp.at(1,1) = 1;
    NpTNp.at(2,2) = 1;
    NpTNp.at(3,3) = 1;
    NpTNp.at(4,4) = 1;
    NpTNp.at(1,3) = -1;
    NpTNp.at(2,4) = -1;
    NpTNp.at(3,1) = -1;
    NpTNp.at(4,2) = -1;

    answer.resize(ndofs,ndofs);
    answer.zero();
    if (dt == _3dAxisymmMode) {
        OOFEM_WARNING("Not tested");
        FloatArray Bh(4);
        Bh.zero();
        Bh.at(1) = 1;
        Bh.at(3) = 1;

        // It was simpler to write this in index notation.
        // Also using 0-based, to reduce typing
        double rJinv = (x1+x2)/length;
        answer.zero();
        for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) {
            answer(i,j) = M_PI*(Ah(i)*Bh(j) + Bh(i)*Ah(j)+ rJinv*(NpTNp(i,j) - Ah(i)*Ah(j)));
        }
    }
    else {
        width = 1;
        answer.beDyadicProductOf(Ah,Ah);
        answer.add(NpTNp);
        answer.times(width/length);
    }

    answer.times(gamma_s);
#endif
}

// Some extension Interfaces to follow:

Interface *LineSurfaceTension :: giveInterface(InterfaceType it)
{
    switch ( it ) {
        case SpatialLocalizerInterfaceType:
            return static_cast< SpatialLocalizerInterface * >(this);
        case EIPrimaryUnknownMapperInterfaceType:
            return static_cast< EIPrimaryUnknownMapperInterface * >(this);
        default:
            return FMElement :: giveInterface(it);
    }
}

} // end namespace oofem
