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
 *               Copyright (C) 1993 - 2010   Borek Patzak
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
#include "domain.h"
#include <math.h>

namespace oofem {

LineSurfaceTension :: LineSurfaceTension(int n, Domain *aDomain) : FMElement (n, aDomain)
{
    numberOfDofMans  = 2;
}
LineSurfaceTension :: ~LineSurfaceTension () {}

IRResultType LineSurfaceTension :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;

    this->Element :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, gamma_s, IFT_SurfaceTension_gamma_s, "gamma_s");


    if (this->boundarySides.giveSize() > 0) { // Why giveSize for IntArray, all other arrays use size()
        if (this->boundarySides.giveSize() != 2) {
            OOFEM_ERROR("Boundary sides must be 2 if exists");
        }
        if (this->boundarySides.at(1)) {
            this->bflag1 = true;
        }
        if (this->boundarySides.at(2)) {
            this->bflag2 = true;
        }
    }

    return IRRT_OK;
}

void LineSurfaceTension :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
// Returns the mask for this element. Should i handle EquationID gracefully?
{
    if (ut == EID_MomentumBalance || ut == EID_MomentumBalance_ConservationEquation) {
        answer.resize(2);
        if (this->giveDomain()->giveDomainType() == _2dIncompressibleFlow) {
            answer.at(1) = V_u;
            answer.at(2) = V_v;
        } else {
            answer.at(1) = D_u;
            answer.at(1) = D_v;
        }
    } else
        answer.resize(0);
}

void LineSurfaceTension :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode, TimeStep *tStep)
{
    if( mtrx == LoadVector )  {
        this->computeLoadVector(answer, mode, tStep);
    } else if ( mtrx == NodalInternalForcesVector ){
        answer.resize(this->computeNumberOfDofs(EID_MomentumBalance));
        answer.zero();
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
    double x1, x2, y1, y2, dx, dy, vx, vy, length, t;

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

    answer.resize(4);
    answer.zero();

    if (dt == _3dAxisymmMode) {
        OOFEM_CLASS_ERROR("Not tested");
        t = M_PI*(x1+x2);
        answer.at(1) = vx*t - length;
        answer.at(2) = vy*t;
        answer.at(3) = -vx*t - length;
        answer.at(4) = -vy*t;
        if (this->bflag1) {
            t = 2*M_PI*x1;
            answer.at(1) -= vx*t;
            answer.at(2) -= vy*t;
        }
        if (this->bflag2) {
            t = 2*M_PI*x2;
            answer.at(3) += vx*t;
            answer.at(4) += vy*t;
        }
    }
    else {
        //t = this->giveDomain()->giveCrossSection(1)->give(CS_Thickness); // TODO: Should i use this?
        t = 1;
        // In this case, the bondary term would cancel out the edge term,
        // so it can be simplified to check when the boundary is not there.
        if (!this->bflag1) {
            answer.at(1) = vx*t;
            answer.at(2) = vy*t;
        }
        if (!this->bflag2) {
            answer.at(3) = -vx*t;
            answer.at(4) = -vy*t;
        }
    }

    answer.times(this->gamma_s);
}

void LineSurfaceTension :: computeTangent(FloatMatrix &answer, TimeStep *tStep)
{
    // TODO: Not sure if it's a good idea to use this tangent.
    answer.resize(4,4);
    answer.zero();

    /*
    domainType dt = this->giveDomain()->giveDomainType();
    int ndofs = this->computeNumberOfDofs(EID_MomentumBalance);
    Node *node1, *node2;
    double x1, x2, y1, y2, dx, dy, vx, vy, length, width;

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
        if (this->bflag1) {
            rJinv = 2*x1/length;
            for (int i = 0; i < 2; i++) for (int j = 0; j < 4; j++) {
                answer(i,j) += -M_PI*(2*Ah(i)*Bh(j) + rJinv*(NpTNp(i,j) - Ah(i)*Ah(j)));
            }
        }
        if (this->bflag2) {
            rJinv = 2*x2/length;
            for (int i = 2; i < 4; i++) for (int j = 0; j < 4; j++) {
                answer(i,j) += -M_PI*(2*Ah(i)*Bh(j) + rJinv*(NpTNp(i,j) - Ah(i)*Ah(j)));
            }
        }
    }
    else {
        width = 1;
        answer.beDyadicProductOf(Ah,Ah);
        answer.add(NpTNp);
        answer.times(width/length);

        // In this special case, the terms will just cancel out on each node respectively.
        if (this->bflag1) {
            answer.at(1,1) = answer.at(1,2) = answer.at(1,3) = answer.at(1,4) = 0;
            answer.at(2,1) = answer.at(2,2) = answer.at(2,3) = answer.at(2,4) = 0;
        }
        if (this->bflag2) {
            answer.at(3,1) = answer.at(3,2) = answer.at(3,3) = answer.at(3,4) = 0;
            answer.at(4,1) = answer.at(4,2) = answer.at(4,3) = answer.at(4,4) = 0;
        }
    }

    answer.times(this->gamma_s);*/
}

} // end namespace oofem
