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

#include "fei2dquadconst.h"
#include "mathfem.h"
#include "flotmtrx.h"
#include "flotarry.h"

namespace oofem {
void
FEI2dQuadConst :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(1);
    answer.at(1) = 1;
}

void
FEI2dQuadConst :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(1, 2);
    answer.at(1, 1) = 0;
    answer.at(1, 2) = 0;
}

void
FEI2dQuadConst :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(2);

    answer.at(1) = 0.25*( cellgeo.giveVertexCoordinates(1)->at(xind) +
                    cellgeo.giveVertexCoordinates(2)->at(xind) +
                    cellgeo.giveVertexCoordinates(3)->at(xind) +
                    cellgeo.giveVertexCoordinates(4)->at(xind) );
    answer.at(2) = 0.25*( cellgeo.giveVertexCoordinates(1)->at(yind) +
                    cellgeo.giveVertexCoordinates(2)->at(yind) +
                    cellgeo.giveVertexCoordinates(3)->at(yind) +
                    cellgeo.giveVertexCoordinates(4)->at(yind) );
}


int
FEI2dQuadConst :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo)
{
    answer.resize(2);
    answer.at(1)=answer.at(2) = 0.0;

    return 1;
}


double
FEI2dQuadConst :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI2dQuadConst :: giveTransformationJacobian: not implemented");
    return 0.0;
}


void
FEI2dQuadConst :: edgeEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI2dQuadConst :: edgeEvalN: not implemented");
}

void
FEI2dQuadConst :: edgeEvaldNds(FloatArray &answer, int iedge,
                           const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI2dQuadConst :: edgeEvaldNds: not implemented");
}

void
FEI2dQuadConst :: edgeLocal2global(FloatArray &answer, int iedge,
                               const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI2dQuadConst :: edgeLocal2global: not implemented");
}


double
FEI2dQuadConst :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    return 0.5 * this->edgeComputeLength(edgeNodes, cellgeo);
}


void
FEI2dQuadConst :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{
    int aNode = 0, bNode = 0;
    edgeNodes.resize(2);

    if ( iedge == 1 ) { // edge between nodes 1 2
        aNode = 1;
        bNode = 2;
    } else if ( iedge == 2 ) { // edge between nodes 2 3
        aNode = 2;
        bNode = 3;
    } else if ( iedge == 3 ) { // edge between nodes 2 3
        aNode = 3;
        bNode = 1;
    } else {
        OOFEM_ERROR2("FEI2dQuadConst :: computeEdgeMapping: wrong egde number (%d)", iedge);
    }

    edgeNodes.at(1) = aNode;
    edgeNodes.at(2) = bNode;

    //OOFEM_ERROR("FEI2dQuadConst :: computeLocalEdgeMapping: not implemented");

}

double
FEI2dQuadConst :: edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo)
{
    double dx, dy;
    int nodeA, nodeB;

    nodeA   = edgeNodes.at(1);
    nodeB   = edgeNodes.at(2);

    dx      = cellgeo.giveVertexCoordinates(nodeB)->at(xind) - cellgeo.giveVertexCoordinates(nodeA)->at(xind);
    dy      = cellgeo.giveVertexCoordinates(nodeB)->at(yind) - cellgeo.giveVertexCoordinates(nodeA)->at(yind);
    return ( sqrt(dx * dx + dy * dy) );
}
} // end namespace oofem
