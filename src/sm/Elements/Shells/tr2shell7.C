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

#include "../sm/Elements/Shells/tr2shell7.h"
#include "../sm/Materials/structuralms.h"
#include "node.h"
#include "load.h"
#include "mathfem.h"
#include "domain.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "fei3dtrquad.h"
#include "boundaryload.h"
#include "vtkxmlexportmodule.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(Tr2Shell7);

FEI3dTrQuad Tr2Shell7 :: interpolation;

IntArray Tr2Shell7 :: orderingDofTypes = {1, 2, 3, 8, 9, 10, 15, 16, 17, 22, 23, 24, 29, 30, 31, 36, 37, 38,
                        4, 5, 6, 11, 12, 13, 18, 19, 20, 25, 26, 27, 32, 33, 34, 39, 40, 41,
                        7, 14, 21, 28, 35, 42};
IntArray Tr2Shell7 :: orderingNodes = {1, 2, 3, 19, 20, 21, 37, 4, 5, 6, 22, 23, 24, 38, 7, 8, 9, 25, 26, 27, 39,
                       10, 11, 12, 28, 29, 30, 40, 13, 14, 15, 31, 32, 33, 41, 16, 17, 18,
                       34, 35, 36, 42};
IntArray Tr2Shell7 :: orderingEdgeNodes = {1, 2, 3, 10, 11, 12, 19, 4, 5, 6, 13, 14, 15, 20, 7, 8, 9, 16, 17, 18, 21};


Tr2Shell7 :: Tr2Shell7(int n, Domain *aDomain) : Shell7Base(n, aDomain)
{
    this->numberOfDofMans = 6;
}

const IntArray &
Tr2Shell7 :: giveOrderingDofTypes() const
{
    return this->orderingDofTypes;
}

const IntArray &
Tr2Shell7 :: giveOrderingNodes() const
{
    return this->orderingNodes;
}
const IntArray &
Tr2Shell7 :: giveOrderingEdgeNodes() const
{
    return this->orderingEdgeNodes;
}

FEInterpolation *Tr2Shell7 :: giveInterpolation() const { return & interpolation; }



void
Tr2Shell7 :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        int nPointsTri  = 6;   // points in the plane
        
        // Layered cross section for bulk integration
        //@todo - must use a cast here since check consistency has not been called yet
        LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( Tr2Shell7 :: giveCrossSection() );
        if ( layeredCS == NULL ) {
            OOFEM_ERROR("Tr2Shell7 only supports layered cross section");
        }
        this->numberOfGaussPoints = layeredCS->giveNumberOfLayers() * nPointsTri * layeredCS->giveNumIntegrationPointsInLayer();
        layeredCS->setupLayeredIntegrationRule(integrationRulesArray, this, nPointsTri);

    }
}




void
Tr2Shell7 :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    if ( iEdge == 1 ) {        // edge between nodes 1-4-2
        answer = {1, 2, 3, 8, 9, 10, 22, 23, 24,  4, 5, 6, 11, 12, 13, 25, 26, 27,   7, 14, 28};
    } else if ( iEdge == 2 ) { // edge between nodes 2-5-3
        answer = {  8, 9, 10, 15, 16, 17, 29, 30, 31,   11, 12, 13, 18, 19, 20, 32, 33, 34,   14, 21, 35};
    } else if ( iEdge == 3 ) { // edge between nodes 3-6-1
        answer = {  15, 16, 17, 1, 2, 3, 36, 37, 38,   18, 19, 20, 4, 5, 6, 39, 40, 41,   21, 7, 42};
    } else {
        OOFEM_ERROR("wrong edge number");
    }
}


void
Tr2Shell7 :: giveSurfaceDofMapping(IntArray &answer, int iSurf) const
{
    answer.resize(42);
    for ( int i = 1; i <= 42; i++ ) {
        answer.at(i) = i;
    }
}


double
Tr2Shell7 :: computeAreaAround(GaussPoint *gp, double xi)
{
    FloatArray G1, G2, temp;
    FloatMatrix Gcov;
    FloatArray lcoords(3);
    lcoords.at(1) = gp->giveNaturalCoordinate(1);
    lcoords.at(2) = gp->giveNaturalCoordinate(2);
    lcoords.at(3) = xi;
    this->evalInitialCovarBaseVectorsAt(lcoords, Gcov);
    G1.beColumnOf(Gcov, 1);
    G2.beColumnOf(Gcov, 2);
    temp.beVectorProductOf(G1, G2);
    double detJ = temp.computeNorm();
    return detJ *gp->giveWeight();
}



double
Tr2Shell7 :: computeVolumeAroundLayer(GaussPoint *gp, int layer)
{
    double detJ;
    FloatMatrix Gcov;
    FloatArray lcoords;
    lcoords = gp->giveNaturalCoordinates();
    this->evalInitialCovarBaseVectorsAt(lcoords, Gcov);
    detJ = Gcov.giveDeterminant() * 0.5 * this->layeredCS->giveLayerThickness(layer);
    return detJ *gp->giveWeight();
}


void
Tr2Shell7 :: compareMatrices(const FloatMatrix &matrix1, const FloatMatrix &matrix2, FloatMatrix &answer)
{
    int ndofs = 42;
    answer.resize(ndofs, ndofs);
    for ( int i = 1; i <= ndofs; i++ ) {
        for ( int j = 1; j <= 18; j++ ) {
            if ( fabs( matrix1.at(i, j) ) > 1.0e-12 ) {
                double diff = ( matrix1.at(i, j) - matrix2.at(i, j) );
                double relDiff =  diff / matrix1.at(i, j);
                if ( fabs(relDiff) < 1.0e-4 ) {
                    answer.at(i, j) = 0.0;
                } else if ( fabs(diff) < 1.0e3 ) {
                    answer.at(i, j) = 0.0;
                } else {
                    answer.at(i, j) = relDiff;
                }
            } else {
                answer.at(i, j) = -1.0;
            }
        }
    }
}
} // end namespace oofem
