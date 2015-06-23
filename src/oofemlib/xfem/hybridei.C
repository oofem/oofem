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

#include "xfem/hybridei.h"
#include "xfemmanager.h"
#include "node.h"
#include "domain.h"
#include "floatmatrix.h"
#include "classfactory.h"

#include <string>


namespace oofem {
REGISTER_EnrichmentItem(HybridEI)

HybridEI :: HybridEI(int n, XfemManager *xm, Domain *aDomain) :
    GeometryBasedEI(n, xm, aDomain)
{}

HybridEI :: ~HybridEI()
{}

void HybridEI :: evalLevelSetNormal(double &oLevelSet, const FloatArray &iGlobalCoord, const FloatArray &iN, const IntArray &iNodeInd) const
{
    interpLevelSet(oLevelSet, iN, iNodeInd);
}

void HybridEI :: evalLevelSetTangential(double &oLevelSet, const FloatArray &iGlobalCoord, const FloatArray &iN, const IntArray &iNodeInd) const
{
    interpLevelSetTangential(oLevelSet, iN, iNodeInd);
}

void HybridEI :: evalGradLevelSetNormal(FloatArray &oGradLevelSet, const FloatArray &iGlobalCoord, const FloatMatrix &idNdX, const IntArray &iNodeInd) const
{
    interpGradLevelSet(oGradLevelSet, idNdX, iNodeInd);
}

void HybridEI :: interpLevelSet(double &oLevelSet, const FloatArray &iN, const IntArray &iNodeInd) const
{
    oLevelSet = 0.0;
    for ( int i = 1; i <= iN.giveSize(); i++ ) {
        double levelSetNode = 0.0;
        const FloatArray &nodePos = this->giveDomain()->giveNode(iNodeInd [ i - 1 ])->giveNodeCoordinates();
        if ( evalLevelSetNormalInNode(levelSetNode, iNodeInd [ i - 1 ], nodePos) ) {
            oLevelSet += iN.at(i) * levelSetNode;
        }
    }
}

void HybridEI :: interpLevelSetTangential(double &oLevelSet, const FloatArray &iN, const IntArray &iNodeInd) const
{
    oLevelSet = 0.0;
    for ( int i = 1; i <= iN.giveSize(); i++ ) {
        double levelSetNode = 0.0;
        const FloatArray &nodePos = this->giveDomain()->giveNode(iNodeInd [ i - 1 ])->giveNodeCoordinates();
        if ( evalLevelSetTangInNode(levelSetNode, iNodeInd [ i - 1 ], nodePos) ) {
            oLevelSet += iN.at(i) * levelSetNode;
        }
    }
}

void HybridEI :: interpGradLevelSet(FloatArray &oGradLevelSet, const FloatMatrix &idNdX, const IntArray &iNodeInd) const
{
    int dim = idNdX.giveNumberOfColumns();

    if ( oGradLevelSet.giveSize() != dim ) {
        oGradLevelSet.resize(dim);
    }

    oGradLevelSet.zero();

    for ( int i = 1; i <= idNdX.giveNumberOfRows(); i++ ) {
        for ( int j = 1; j <= dim; j++ ) {
            double levelSetNode = 0.0;
            const FloatArray &nodePos = this->giveDomain()->giveNode(iNodeInd [ i - 1 ])->giveNodeCoordinates();
            if ( evalLevelSetNormalInNode(levelSetNode, iNodeInd [ i - 1 ], nodePos) ) {
                oGradLevelSet.at(j) += idNdX.at(i, j) * levelSetNode;
            }
        }
    }
}
} /* namespace oofem */
