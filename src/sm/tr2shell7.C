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

#include "tr2shell7.h"
#include "node.h"
#include "load.h"
#include "structuralms.h"
#include "mathfem.h"
#include "domain.h"
#include "equationid.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "fei3dtrquad.h"
#include "boundaryload.h"
#include "vtkxmlexportmodule.h"
#include "classfactory.h"

namespace oofem {

REGISTER_Element( Tr2Shell7 );

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


FEInterpolation *Tr2Shell7 :: giveInterpolation() const { return & interpolation; }



void
Tr2Shell7 :: computeGaussPoints()
{
    if ( !integrationRulesArray ) {
        int nPointsTri = 6;                     // points in the plane
        int nPointsEdge = 2;
        integrationRulesArray = new IntegrationRule * [ 3 ];

        // Midplane and thickness

        // Midplane only (Mass matrix integrated analytically through the thickness)
        integrationRulesArray [ 1 ] = new GaussIntegrationRule(1, this);
        integrationRulesArray [ 1 ]->SetUpPointsOnWedge(nPointsTri, 1, _3dMat);


        // Edge
        integrationRulesArray [ 2 ] = new GaussIntegrationRule(1, this);
        integrationRulesArray [ 2 ]->SetUpPointsOnLine(nPointsEdge, _3dMat);


        // Layered cross section
        LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( Tr2Shell7 :: giveCrossSection() );
        if ( layeredCS == NULL ) {
            OOFEM_ERROR("Tr2Shell7 only supports layered cross section");
        }

        int numberOfLayers = layeredCS->giveNumberOfLayers();
        layerIntegrationRulesArray = new IntegrationRule * [ numberOfLayers ];

        // may need to extend this to handle Newton-Cotes integration in the thickness direction
        for ( int i = 1; i <= numberOfLayers; i++ ) {
            layerIntegrationRulesArray [ i - 1 ] = new GaussIntegrationRule(1, this);
            layerIntegrationRulesArray [ i - 1 ]->SetUpPointsOnWedge(nPointsTri, layeredCS->giveNumIntegrationPointsInLayer(), _3dMat);
        }

        layeredCS->mapLayerGpCoordsToShellCoords(layeredCS, layerIntegrationRulesArray);
        layeredCS->printYourself();
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


// Integration

double
Tr2Shell7 :: computeVolumeAround(GaussPoint *gp)
{
    FloatArray G1, G2, G3, temp;
    double detJ;
    this->evalInitialCovarBaseVectorsAt(gp, G1, G2, G3);
    temp.beVectorProductOf(G1, G2);
    detJ = temp.dotProduct(G3) * 0.5 * this->giveCrossSection()->give(CS_Thickness);
    return detJ * gp->giveWeight();
}

double
Tr2Shell7 :: computeAreaAround(GaussPoint *gp)
{
    FloatArray G1, G2, G3, temp;
    this->evalInitialCovarBaseVectorsAt(gp, G1, G2, G3);
    temp.beVectorProductOf(G1, G2);
    double detJ = temp.computeNorm();
    return detJ * gp->giveWeight() * 0.5;
}




double
Tr2Shell7 :: computeVolumeAroundLayer(GaussPoint *gp, int layer)
{
    FloatArray G1, G2, G3, temp;
    double detJ;
    this->evalInitialCovarBaseVectorsAt(gp, G1, G2, G3);
    temp.beVectorProductOf(G1, G2);
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveCrossSection() );
    detJ = temp.dotProduct(G3) * 0.5 * layeredCS->giveLayerThickness(layer);
    return detJ * gp->giveWeight();
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
            } else{
                answer.at(i, j) = -1.0;
            }
        }
    }
}


void 
Tr2Shell7 :: vtkGiveFictiousNodeCoords(FloatArray nodeCoords[15], int layer)
{
    // compute fictious node coords
    FloatArray nodeLocalXiCoords, nodeLocalEtaCoords;
    this->giveLocalNodeCoords( nodeLocalXiCoords, nodeLocalEtaCoords);

    for ( int i = 1; i <= 3; i++ ){
        FloatArray coords, localCoords(3);
        localCoords.at(1) = nodeLocalXiCoords.at(i);
        localCoords.at(2) = nodeLocalEtaCoords.at(i);

        localCoords.at(3) = -1.0;
        this->vtkEvalInitialGlobalCoordinateAt(localCoords, layer, coords);
        nodeCoords[i-1].resize(3); 
        nodeCoords[i-1].at(1) = coords.at(1);
        nodeCoords[i-1].at(2) = coords.at(2);
        nodeCoords[i-1].at(3) = coords.at(3);

        localCoords.at(3) = 1.0;
        this->vtkEvalInitialGlobalCoordinateAt(localCoords, layer, coords);
        nodeCoords[i+3-1].resize(3); 
        nodeCoords[i+3-1].at(1) = coords.at(1);
        nodeCoords[i+3-1].at(2) = coords.at(2);
        nodeCoords[i+3-1].at(3) = coords.at(3);
    }

    for ( int i = 1; i <= 3; i++ ){
        FloatArray coords, localCoords(3);
        localCoords.at(1) = nodeLocalXiCoords.at(i+3);
        localCoords.at(2) = nodeLocalEtaCoords.at(i+3);

        localCoords.at(3) = -1.0;
        this->vtkEvalInitialGlobalCoordinateAt(localCoords, layer, coords);
        nodeCoords[i+6-1].resize(3); 
        nodeCoords[i+6-1].at(1) = coords.at(1);
        nodeCoords[i+6-1].at(2) = coords.at(2);
        nodeCoords[i+6-1].at(3) = coords.at(3);

        localCoords.at(3) = 1.0;
        this->vtkEvalInitialGlobalCoordinateAt(localCoords, layer, coords);
        nodeCoords[i+6+3-1].resize(3); 
        nodeCoords[i+6+3-1].at(1) = coords.at(1);
        nodeCoords[i+6+3-1].at(2) = coords.at(2);
        nodeCoords[i+6+3-1].at(3) = coords.at(3);
    }
    

    //Middle nodes - linear variation in the thickness direction -> mean value of top and bottom
        for ( int i = 1; i <= 3; i++ ){
        nodeCoords[12+i-1].resize(3);
        nodeCoords[12+i-1].at(1) = 0.5 * ( nodeCoords[i-1].at(1) + nodeCoords[i+3-1].at(1) );
        nodeCoords[12+i-1].at(2) = 0.5 * ( nodeCoords[i-1].at(2) + nodeCoords[i+3-1].at(2) );
        nodeCoords[12+i-1].at(3) = 0.5 * ( nodeCoords[i-1].at(3) + nodeCoords[i+3-1].at(3) );
        }
}


void 
Tr2Shell7 :: vtkGiveUpdatedFictiousNodeCoords(FloatArray nodeCoords[15], int layer, TimeStep *tStep)
{
    // compute fictious node coords
    FloatArray nodeLocalXiCoords, nodeLocalEtaCoords;
    this->giveLocalNodeCoords( nodeLocalXiCoords, nodeLocalEtaCoords);

    for ( int i = 1; i <= 3; i++ ){
        FloatArray coords, localCoords(3);
        localCoords.at(1) = nodeLocalXiCoords.at(i);
        localCoords.at(2) = nodeLocalEtaCoords.at(i);

        localCoords.at(3) = -1.0;
        this->vtkEvalUpdatedGlobalCoordinateAt(localCoords, layer, coords, tStep);
        nodeCoords[i-1].resize(3); 
        nodeCoords[i-1].at(1) = coords.at(1);
        nodeCoords[i-1].at(2) = coords.at(2);
        nodeCoords[i-1].at(3) = coords.at(3);

        localCoords.at(3) = 1.0;
        this->vtkEvalUpdatedGlobalCoordinateAt(localCoords, layer, coords, tStep);
        nodeCoords[i+3-1].resize(3); 
        nodeCoords[i+3-1].at(1) = coords.at(1);
        nodeCoords[i+3-1].at(2) = coords.at(2);
        nodeCoords[i+3-1].at(3) = coords.at(3);
    }

    for ( int i = 1; i <= 3; i++ ){
        FloatArray coords, localCoords(3);
        localCoords.at(1) = nodeLocalXiCoords.at(i+3);
        localCoords.at(2) = nodeLocalEtaCoords.at(i+3);

        localCoords.at(3) = -1.0;
        this->vtkEvalUpdatedGlobalCoordinateAt(localCoords, layer, coords, tStep);
        nodeCoords[i+6-1].resize(3); 
        nodeCoords[i+6-1].at(1) = coords.at(1);
        nodeCoords[i+6-1].at(2) = coords.at(2);
        nodeCoords[i+6-1].at(3) = coords.at(3);

        localCoords.at(3) = 1.0;
        this->vtkEvalUpdatedGlobalCoordinateAt(localCoords, layer, coords, tStep);
        nodeCoords[i+6+3-1].resize(3); 
        nodeCoords[i+6+3-1].at(1) = coords.at(1);
        nodeCoords[i+6+3-1].at(2) = coords.at(2);
        nodeCoords[i+6+3-1].at(3) = coords.at(3);
    }
    

    //Middle nodes - linear variation in the thickness direction -> mean value of top and bottom
        for ( int i = 1; i <= 3; i++ ){
        nodeCoords[12+i-1].resize(3);
        nodeCoords[12+i-1].at(1) = 0.5 * ( nodeCoords[i-1].at(1) + nodeCoords[i+3-1].at(1) );
        nodeCoords[12+i-1].at(2) = 0.5 * ( nodeCoords[i-1].at(2) + nodeCoords[i+3-1].at(2) );
        nodeCoords[12+i-1].at(3) = 0.5 * ( nodeCoords[i-1].at(3) + nodeCoords[i+3-1].at(3) );
        }

}


} // end namespace oofem
