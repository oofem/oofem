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

#include "tr2shell7xfem.h"
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
#include "classfactory.h"

#include "tr2shell7.h"

namespace oofem {

REGISTER_Element( Tr2Shell7XFEM );

FEI3dTrQuad Tr2Shell7XFEM :: interpolation;

IntArray Tr2Shell7XFEM :: ordering_all(42);
IntArray Tr2Shell7XFEM :: ordering_gr(42);
IntArray Tr2Shell7XFEM :: ordering_gr_edge(21);
bool Tr2Shell7XFEM :: __initialized = Tr2Shell7XFEM :: initOrdering();



Tr2Shell7XFEM :: Tr2Shell7XFEM(int n, Domain *aDomain) : Shell7BaseXFEM(n, aDomain)
{
    this->numberOfDofMans = 6;
}

const IntArray &
Tr2Shell7XFEM :: giveOrdering(SolutionField fieldType) const
{
    if ( fieldType == All ) {
        return this->ordering_all;
    } else if ( fieldType == AllInv ) {
        return this->ordering_gr;
    } else /*if ( fieldType == EdgeInv  )*/ {
        return this->ordering_gr_edge;
    }
}


void 
Tr2Shell7XFEM :: giveLocalNodeCoords(FloatArray &nodeLocalXiCoords, FloatArray &nodeLocalEtaCoords)
{
    nodeLocalXiCoords.setValues( 6, 1., 0., 0., .5, 0., .5); // corner nodes then midnodes, uncertain of node numbering
    nodeLocalEtaCoords.setValues(6, 0., 1., 0., .5, .5, 0.);
}


FEInterpolation* Tr2Shell7XFEM :: giveInterpolation() const { return &interpolation; }


void 
Tr2Shell7XFEM :: computeGaussPoints()
{
    if ( !integrationRulesArray ) {
        //int nPointsTri  = 6;   // points in the plane
        int nPointsTri  = 25;   // points in the plane
        int nPointsEdge = 2;   // edge integration

        // need to check if interface has failed but also need to update the integration rule later
        XfemManager *xMan = this->giveDomain()->giveXfemManager();
        for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) { 
            Delamination *dei =  dynamic_cast< Delamination * >( xMan->giveEnrichmentItem(i) ); 
            if (dei) {
                int numberOfInterfaces = this->layeredCS->giveNumberOfLayers()-1;
                czIntegrationRulesArray = new IntegrationRule * [ numberOfInterfaces ];
                for ( int i = 0; i < numberOfInterfaces; i++ ) {
                    czIntegrationRulesArray [ i ] = new GaussIntegrationRule(1, this);
                    czIntegrationRulesArray [ i ]->SetUpPointsOnTriangle(nPointsTri, _3dInterface);
                    //czIntegrationRulesArray [ i ]->SetUpPointsOnTriangle(3, _3dInterface);
                }
            }
        }
        

        specialIntegrationRulesArray = new IntegrationRule * [ 3 ];

        // Midplane (Mass matrix integrated analytically through the thickness)
        specialIntegrationRulesArray [ 1 ] = new GaussIntegrationRule(1, this);
        specialIntegrationRulesArray [ 1 ]->SetUpPointsOnTriangle(nPointsTri, _3dMat); 

        // Edge
        specialIntegrationRulesArray [ 2 ] = new GaussIntegrationRule(1, this);
        specialIntegrationRulesArray [ 2 ]->SetUpPointsOnLine(nPointsEdge, _3dMat);
        
        // Layered cross section for bulk integration
        this->numberOfIntegrationRules = this->layeredCS->giveNumberOfLayers();
        this->numberOfGaussPoints = this->layeredCS->giveNumberOfLayers()*nPointsTri*this->layeredCS->giveNumIntegrationPointsInLayer();
        this->layeredCS->setupLayeredIntegrationRule(integrationRulesArray, this, nPointsTri);

        // Thickness integration for stress recovery
        specialIntegrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this);
        specialIntegrationRulesArray [ 0 ]->SetUpPointsOnLine(this->layeredCS->giveNumIntegrationPointsInLayer(), _3dMat);




    }

}


void
Tr2Shell7XFEM :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    if ( iEdge == 1 )        { // edge between nodes 1-4-2
        answer.setValues(21, 1, 2, 3, 8, 9, 10, 22, 23, 24,  4, 5, 6, 11, 12, 13, 25, 26, 27,   7, 14, 28);

    } else if ( iEdge == 2 ) { // edge between nodes 2-5-3
        answer.setValues(21,   8, 9, 10, 15, 16, 17, 29, 30, 31,   11, 12, 13, 18, 19, 20, 32, 33, 34,   14, 21, 35 );

    } else if ( iEdge == 3 ) { // edge between nodes 3-6-1
        answer.setValues(21,   15, 16, 17, 1, 2, 3, 36, 37, 38,   18, 19, 20, 4, 5, 6, 39, 40, 41,   21, 7, 42);
    } else {
        _error("giveEdgeDofMapping: wrong edge number");
    }
}


void 
Tr2Shell7XFEM :: giveSurfaceDofMapping(IntArray &answer, int iSurf) const
{
    answer.resize(42);
    for( int i = 1; i <= 42; i++){    
        answer.at(i) = i;
    }
}


// Integration

double 
Tr2Shell7XFEM :: computeAreaAround(GaussPoint *gp, double xi)
{
    FloatArray G1, G2, temp;
    FloatMatrix Gcov;
    FloatArray lCoords(3);
    lCoords.at(1) = gp->giveCoordinate(1);
    lCoords.at(2) = gp->giveCoordinate(2);
    lCoords.at(3) = xi;
    this->evalInitialCovarBaseVectorsAt(lCoords, Gcov);
    G1.beColumnOf(Gcov,1);
    G2.beColumnOf(Gcov,2);
    temp.beVectorProductOf(G1, G2);
    double detJ = temp.computeNorm();
    return detJ * gp->giveWeight() ;
}


double 
Tr2Shell7XFEM :: computeVolumeAroundLayer(GaussPoint *gp, int layer)
{
    double detJ;
    FloatMatrix Gcov;
    FloatArray lcoords;
    lcoords = *gp->giveCoordinates();
    this->evalInitialCovarBaseVectorsAt(lcoords, Gcov);
    detJ = Gcov.giveDeterminant() * 0.5 * this->layeredCS->giveLayerThickness(layer);
    return detJ * gp->giveWeight();
}



} // end namespace oofem

