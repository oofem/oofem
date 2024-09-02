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

#include "fei1dlin.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "floatmatrixf.h"
#include "floatarrayf.h"

namespace oofem {
double
FEI1dLin :: giveLength(const FEICellGeometry &cellgeo) const
{
    return fabs( cellgeo.giveVertexCoordinates(2).at(cindx) - cellgeo.giveVertexCoordinates(1).at(cindx) );
}

FloatArrayF<2>
FEI1dLin :: evalN(double ksi) 
{
    return {( 1. - ksi ) * 0.5, ( 1. + ksi ) * 0.5};
}

void
FEI1dLin :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double ksi = lcoords.at(1);
    answer.resize(2);

    answer.at(1) = ( 1. - ksi ) * 0.5;
    answer.at(2) = ( 1. + ksi ) * 0.5;
}

std::pair<double, FloatMatrixF<1,2>>
FEI1dLin :: evaldNdx(const FEICellGeometry &cellgeo) const
{
    double l = cellgeo.giveVertexCoordinates(2).at(cindx) - cellgeo.giveVertexCoordinates(1).at(cindx);
    return {0.5 * l, {-1.0 / l, 1.0 / l}};
}

double
FEI1dLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double l = cellgeo.giveVertexCoordinates(2).at(cindx) - cellgeo.giveVertexCoordinates(1).at(cindx);
    answer.resize(2, 1);

    answer.at(1, 1) = -1.0 / l;
    answer.at(2, 1) =  1.0 / l;
    return 0.5 * l;
}

void
FEI1dLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    FloatArray n;
    answer.resize(1);

    this->evalN(n, lcoords, cellgeo);
    answer.at(1) = n.at(1) * cellgeo.giveVertexCoordinates(1).at(cindx) +
                   n.at(2) * cellgeo.giveVertexCoordinates(2).at(cindx);
}

int
FEI1dLin :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo) const
{
    double x1 = cellgeo.giveVertexCoordinates(1).at(cindx);
    double x2 = cellgeo.giveVertexCoordinates(2).at(cindx);
    double ksi = ( 2.0 * coords.at(1) - ( x1 + x2 ) ) / ( x2 - x1 );
    answer.resize(1);
    answer.at(1) = clamp(ksi, -1., 1.);
    return fabs(ksi) <= 1.0;
}

double
FEI1dLin :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    return 0.5 * ( cellgeo.giveVertexCoordinates(2).at(cindx) - cellgeo.giveVertexCoordinates(1).at(cindx) );
}

IntArray FEI1dLin :: boundaryEdgeGiveNodes(int boundary, Element_Geometry_Type egt) const
{
    return {1, 2};
}

void FEI1dLin :: boundaryEdgeEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    this->evalN(answer, lcoords, cellgeo);
}

double FEI1dLin :: boundaryEdgeGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const 
{
    return this->giveTransformationJacobian(lcoords, cellgeo);
}

void FEI1dLin :: boundaryEdgeLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    this->local2global(answer, lcoords, cellgeo);
}



} // end namespace oofem
