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

#include "tr2shell7.h"
#include "node.h"
#include "load.h"
#include "structuralms.h"
#include "mathfem.h"
#include "domain.h"
#include "equationid.h"
#include "gaussintegrationrule.h"
#include "gausspnt.h"
#include "fei3dtrquad.h"
#include "boundaryload.h"
#include "vtkxmlexportmodule.h"

namespace oofem {
FEI3dTrQuad Tr2Shell7 :: interpolation;

IntArray Tr2Shell7 :: ordering_phibar(18);
IntArray Tr2Shell7 :: ordering_m(18);
IntArray Tr2Shell7 :: ordering_gam(6);
IntArray Tr2Shell7 :: ordering_all(42);
IntArray Tr2Shell7 :: ordering_gr(42);
IntArray Tr2Shell7 :: ordering_gr_edge(21);
bool Tr2Shell7 :: __initialized = Tr2Shell7 :: initOrdering();



Tr2Shell7 :: Tr2Shell7(int n, Domain *aDomain) : Shell7Base(n, aDomain)
{
    this->numberOfDofMans = 6;
}

const IntArray &
Tr2Shell7 :: giveOrdering(SolutionField fieldType) const
{
    if ( fieldType == Midplane ) {
        return this->ordering_phibar;
    } else if ( fieldType == Director ) {
        return this->ordering_m;
    } else if ( fieldType == InhomStrain ) {
        return this->ordering_gam;
    } else if ( fieldType == All ) {
        return this->ordering_all;
    } else if ( fieldType == AllInv ) {
        return this->ordering_gr;
    } else /*if ( fieldType == EdgeInv )*/ {
        return this->ordering_gr_edge;
    }
}


void
Tr2Shell7 :: giveLocalNodeCoords(FloatArray &nodeLocalXiCoords, FloatArray &nodeLocalEtaCoords)
{
    nodeLocalXiCoords.setValues(6, 1., 0., 0., .5, 0., .5);      // corner nodes then midnodes, uncertain of node numbering
    nodeLocalEtaCoords.setValues(6, 0., 1., 0., .5, .5, 0.);
}


FEInterpolation *Tr2Shell7 :: giveInterpolation() { return & interpolation; }



void
Tr2Shell7 :: computeGaussPoints()
{
    if ( !integrationRulesArray ) {
        int nPointsTri  = 6;   // points in the plane
        int nPointsEdge = 2;   // edge integration
        specialIntegrationRulesArray = new IntegrationRule * [ 3 ];

        // Midplane and thickness

        // Midplane (Mass matrix integrated analytically through the thickness)
        specialIntegrationRulesArray [ 1 ] = new GaussIntegrationRule(1, this);
        specialIntegrationRulesArray [ 1 ]->SetUpPointsOnWedge(nPointsTri, 1, _3dMat); //@todo replce with triangle which has a xi3-coord


        // Edge
        specialIntegrationRulesArray [ 2 ] = new GaussIntegrationRule(1, this);
        specialIntegrationRulesArray [ 2 ]->SetUpPointsOnLine(nPointsEdge, _3dMat);


        // Layered cross section for bulk integration
        //@todo - must use a cast here since check consistency has not been called yet
        LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( Tr2Shell7 :: giveCrossSection() );
        if ( layeredCS == NULL ) {
            OOFEM_ERROR("Tr2Shell7 only supports layered cross section");
        }
        this->numberOfIntegrationRules = layeredCS->giveNumberOfLayers();
        this->numberOfGaussPoints = layeredCS->giveNumberOfLayers()*nPointsTri*layeredCS->giveNumIntegrationPointsInLayer();
        layeredCS->setupLayeredIntegrationRule(integrationRulesArray, this, nPointsTri);
        //layeredCS->printYourself();

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
        //answer.setValues(21, 1, 2, 3, 4, 5, 6, 7,   8, 9, 10, 11, 12, 13, 14,   22, 23, 24, 25, 26, 27, 28);
        answer.setValues(21, 1, 2, 3, 8, 9, 10, 22, 23, 24,  4, 5, 6, 11, 12, 13, 25, 26, 27,   7, 14, 28);
    } else if ( iEdge == 2 ) { // edge between nodes 2-5-3
        //answer.setValues(21,   8, 9, 10, 11, 12, 13, 14,   15, 16, 17, 18, 19, 20, 21,   29, 30, 31, 32, 33, 34, 35 );
        answer.setValues(21,   8, 9, 10, 15, 16, 17, 29, 30, 31,   11, 12, 13, 18, 19, 20, 32, 33, 34,   14, 21, 35);
    } else if ( iEdge == 3 ) { // edge between nodes 3-6-1
        //answer.setValues(21,   15, 16, 17, 18, 19, 20, 21,   1, 2, 3, 4, 5, 6, 7,   36, 37, 38, 39, 40, 41, 42);
        answer.setValues(21,   15, 16, 17, 1, 2, 3, 36, 37, 38,   18, 19, 20, 4, 5, 6, 39, 40, 41,   21, 7, 42);
    } else {
        _error("giveEdgeDofMapping: wrong edge number");
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
Tr2Shell7 :: computeAreaAround(GaussPoint *gp)
{
    FloatArray G1, G2, temp;
    FloatMatrix Gcov;
    this->evalInitialCovarBaseVectorsAt(gp, Gcov);
    G1.beColumnOf(Gcov,1);
    G2.beColumnOf(Gcov,2);
    temp.beVectorProductOf(G1, G2);
    double detJ = temp.computeNorm();
    return detJ * gp->giveWeight() * 0.5;
}



double
Tr2Shell7 :: computeVolumeAroundLayer(GaussPoint *gp, int layer)
{
    double detJ;
    FloatMatrix Gcov;
    this->evalInitialCovarBaseVectorsAt(gp, Gcov);
    detJ = Gcov.giveDeterminant() * 0.5 * this->layeredCS->giveLayerThickness(layer);
    return detJ * gp->giveWeight();
}


void
Tr2Shell7 :: compareMatrices(const FloatMatrix &matrix1, const FloatMatrix &matrix2, FloatMatrix &answer)
{
    int ndofs = 42;
    answer.resize(ndofs, ndofs);
    for ( int i = 1; i <= ndofs; i++ ) {
        for ( int j = 1; j <= 18; j++ ) {
            if ( abs( matrix1.at(i, j) ) > 1.0e-12 ) {
                double diff = ( matrix1.at(i, j) - matrix2.at(i, j) );
                double relDiff =  diff / matrix1.at(i, j);
                if ( abs(relDiff) < 1.0e-4 ) {
                    answer.at(i, j) = 0.0;
                } else if ( abs(diff) < 1.0e3 ) {
                    answer.at(i, j) = 0.0;
                } else {
                    answer.at(i, j) = relDiff;
                }
            } else{
                answer.at(i, j) = -1.0;
            }
        }
    }
}




} // end namespace oofem
