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

#include "discsegintegrationrule.h"

#include "gausspoint.h"
#include "fei1dlin.h"

namespace oofem {
DiscontinuousSegmentIntegrationRule :: DiscontinuousSegmentIntegrationRule(int n, Element *e, const std :: vector< Line > &iSegments, const FloatArray &iXS, const FloatArray &iXE) :
    GaussIntegrationRule(n, e),
    mSegments(iSegments),
    mXS(iXS), mXE(iXE)
{}

DiscontinuousSegmentIntegrationRule :: ~DiscontinuousSegmentIntegrationRule() {}

int DiscontinuousSegmentIntegrationRule :: SetUpPointsOnLine(int iNumPointsPerSeg, MaterialMode mode)
{
    int numPointsTot = iNumPointsPerSeg * mSegments.size();
    int pointsPassed = 0;

    ////////////////////////////////////////////
    // Allocate Gauss point array
    FloatArray coords_xi, weights;
    this->giveLineCoordsAndWeights(iNumPointsPerSeg, coords_xi, weights);
    this->gaussPoints.resize(numPointsTot);
    ////////////////////////////////////////////

    double totalLength = mXS.distance(mXE);

    std :: vector< FloatArray >newGPCoord;

    // Loop over line segments
    for ( size_t i = 0; i < mSegments.size(); i++ ) {
        for ( int j = 0; j < iNumPointsPerSeg; j++ ) {
            FloatArray global;
            GaussPoint * &gp = this->gaussPoints [ pointsPassed ];

            gp = new GaussPoint(this, pointsPassed + 1, {coords_xi.at(j + 1)}, weights.at(j + 1), mode);

            const FloatArray &coord = gp->giveNaturalCoordinates();

            global.resize( mXS.giveSize() );
            for ( int m = 1; m <= mXS.giveSize(); m++ ) {
                global.at(m) = 0.5 * ( ( 1.0 - coord.at(1) ) * mSegments [ i ].giveVertex(1).at(m) + ( 1.0 + coord.at(1) ) * mSegments [ i ].giveVertex(2).at(m) );
            }

            newGPCoord.push_back(global);


            // Local coordinate along the line segment
            double xi = 2.0 * ( global.distance(mXS) / totalLength - 0.5 );
            gp->setNaturalCoordinates({ xi });

            gp->setSubPatchCoordinates({ xi });
            gp->setGlobalCoordinates(global);

            gp->setWeight(1.0 * gp->giveWeight() * mSegments [ i ].giveLength() / totalLength);  // update integration weight

            pointsPassed++;
        }
    }

    return this->giveNumberOfIntegrationPoints();
}
} /* namespace oofem */
