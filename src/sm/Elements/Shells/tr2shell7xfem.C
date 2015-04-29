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

#include "../sm/Elements/Shells/tr2shell7xfem.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/xfem/enrichmentitems/crack.h"
#include "../sm/xfem/enrichmentitems/shellcrack.h"
#include "node.h"
#include "load.h"
#include "mathfem.h"
#include "domain.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "fei3dtrquad.h"
#include "boundaryload.h"
#include "classfactory.h"
#include "xfem/patchintegrationrule.h"
#include "xfem/XFEMDebugTools.h"
#include <string>
#include <sstream>


namespace oofem {

REGISTER_Element( Tr2Shell7XFEM );

FEI3dTrQuad Tr2Shell7XFEM :: interpolation;

IntArray Tr2Shell7XFEM :: orderingDofTypes = {1, 2, 3, 8, 9, 10, 15, 16, 17, 22, 23, 24, 29, 30, 31, 36, 37, 38,
                        4, 5, 6, 11, 12, 13, 18, 19, 20, 25, 26, 27, 32, 33, 34, 39, 40, 41,
                        7, 14, 21, 28, 35, 42};
IntArray Tr2Shell7XFEM :: orderingNodes = {1, 2, 3, 19, 20, 21, 37, 4, 5, 6, 22, 23, 24, 38, 7, 8, 9, 25, 26, 27, 39,
                       10, 11, 12, 28, 29, 30, 40, 13, 14, 15, 31, 32, 33, 41, 16, 17, 18,
                       34, 35, 36, 42};
IntArray Tr2Shell7XFEM :: orderingEdgeNodes = {1, 2, 3, 10, 11, 12, 19, 4, 5, 6, 13, 14, 15, 20, 7, 8, 9, 16, 17, 18, 21};


Tr2Shell7XFEM :: Tr2Shell7XFEM(int n, Domain *aDomain) : Shell7BaseXFEM(n, aDomain)
{
    this->numberOfDofMans = 6;
}

const IntArray &
Tr2Shell7XFEM :: giveOrderingDofTypes() const
{
    return this->orderingDofTypes;
}

const IntArray &
Tr2Shell7XFEM :: giveOrderingNodes() const
{
    return this->orderingNodes;
}
const IntArray &
Tr2Shell7XFEM :: giveOrderingEdgeNodes() const
{
    return this->orderingEdgeNodes;
}


void 
Tr2Shell7XFEM :: giveLocalNodeCoords(FloatArray &nodeLocalXiCoords, FloatArray &nodeLocalEtaCoords)
{
    nodeLocalXiCoords = {1., 0., 0., .5, 0., .5};  // corner nodes then midnodes, uncertain of node numbering
    nodeLocalEtaCoords = {0., 1., 0., .5, .5, 0.};
}


FEInterpolation* Tr2Shell7XFEM :: giveInterpolation() const { return &interpolation; }


void 
Tr2Shell7XFEM :: computeGaussPoints()
{

    this->xMan = this->giveDomain()->giveXfemManager();

    if ( integrationRulesArray.size() == 0 ) {  
        //if( this->xMan->isElementEnriched(this) ) {
            this->updateIntegrationRuleMultiCrack();
        //}
    }

    int nPointsTri  = 6;   // points in the plane
    if ( integrationRulesArray.size() == 0 ) {
        // Cohesive zone
        for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { 
            int numberOfInterfaces = this->layeredCS->giveNumberOfLayers()-1;
            czIntegrationRulesArray.resize( numberOfInterfaces );
            for ( int j = 0; j < numberOfInterfaces; j++ ) {
                czIntegrationRulesArray [ j ].reset( new GaussIntegrationRule(1, this) );
                czIntegrationRulesArray [ j ]->SetUpPointsOnTriangle(nPointsTri, _3dInterface);
            }
        
        }
    }


}




bool Tr2Shell7XFEM :: updateIntegrationRuleMultiCrack()
{
    int nPointsTri  = 6;   // points in the plane

    
    int numberOfLayers     = this->layeredCS->giveNumberOfLayers();
    int numPointsThickness = this->layeredCS->giveNumIntegrationPointsInLayer();
    double totalThickness  = this->layeredCS->computeIntegralThick();
    int numEI = this->xMan->giveNumberOfEnrichmentItems();
    std :: vector< std :: vector< FloatArray > >pointPartitions;

    integrationRulesArray.resize(numberOfLayers);
    this->crackSubdivisions.resize(numberOfLayers);     // Store the subdivisions for each layer (empty otherwise)
    this->numSubDivisionsArray.resize(numberOfLayers);  // number of subdivision for each - set this to one for layers not subdivided

    bool createdRule;
    for ( int i = 0; i < numberOfLayers; i++ ) {
        createdRule = false;
        
        double zMid_i = this->layeredCS->giveLayerMidZ(i+1); // global z-coord
        double xiMid_i = 1.0 - 2.0 * ( totalThickness - this->layeredCS->giveMidSurfaceZcoordFromBottom() - zMid_i ) / totalThickness; // local xi-coord
        
        if ( this->xMan->isElementEnriched(this) ) {
        
            for ( int eiIndex = 1; eiIndex <= numEI; eiIndex++ ) {
                EnrichmentItem *ei = this->xMan->giveEnrichmentItem(eiIndex);
                if ( dynamic_cast< ShellCrack*> (ei) ) {

                    // Determine if the crack goes through the current layer
                    if( this->evaluateHeavisideGamma(xiMid_i, static_cast< ShellCrack* >(ei)) > 0) {

                        // Get the points describing each subdivision of the element
                        double startXi, endXi;
                        bool intersection = false;
                        pointPartitions.resize(0);
                        this->XfemElementInterface_prepareNodesForDelaunay(pointPartitions, startXi, endXi, eiIndex, intersection);

                        if ( intersection ) {
                            for ( int j = 0; j < int( pointPartitions.size() ); j++ ) {
                                // Triangulate the subdivisions
                                this->XfemElementInterface_partitionElement(this->crackSubdivisions [ i ], pointPartitions [ j ]);
                            }

                     
                            integrationRulesArray [ i ].reset( new PatchIntegrationRule(i + 1, this, this->crackSubdivisions [ i ]) );
                            int nPointsTriSubTri = 3; 
                            integrationRulesArray [ i ]->SetUpPointsOnWedge(nPointsTriSubTri, numPointsThickness, _3dMat);         
                            this->numSubDivisionsArray [ i ] = this->crackSubdivisions [ i ].size();
                            createdRule = true;         
                            continue;
                        }
                    }            
                }
            }
        }
            
        if( !createdRule ) {
            integrationRulesArray [ i ].reset( new LayeredIntegrationRule(i + 1, this) );
            integrationRulesArray [ i ]->SetUpPointsOnWedge(nPointsTri, numPointsThickness, _3dMat);
            this->numSubDivisionsArray [ i ] = 1;                 
        }
        
    }
    
    this->layeredCS->mapLayerGpCoordsToShellCoords(integrationRulesArray);
    
    // Cohesive zone
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { 
        //Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) ); 
        int numberOfInterfaces = this->layeredCS->giveNumberOfLayers()-1;
        czIntegrationRulesArray.resize(numberOfInterfaces);
        for ( int j = 0; j < numberOfInterfaces; j++ ) {
            czIntegrationRulesArray [ j ].reset( new GaussIntegrationRule(1, this) );
            czIntegrationRulesArray [ j ]->SetUpPointsOnTriangle(nPointsTri, _3dInterface);
        }
    }


    return true;
}


double 
Tr2Shell7XFEM :: computeArea()
{
  
  //FEIElementGeometryWrapper(this);
  return this->interpolation.giveArea( FEIElementGeometryWrapper(this) );
}

void
Tr2Shell7XFEM :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    if ( iEdge == 1 )        { // edge between nodes 1-4-2
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
    lCoords.at(1) = gp->giveNaturalCoordinate(1);
    lCoords.at(2) = gp->giveNaturalCoordinate(2);
    lCoords.at(3) = xi;
    this->evalInitialCovarBaseVectorsAt(lCoords, Gcov);
    G1.beColumnOf(Gcov,1);
    G2.beColumnOf(Gcov,2);
    temp.beVectorProductOf(G1, G2);
    double detJ = temp.computeNorm();
    return detJ *gp->giveWeight();
}


double 
Tr2Shell7XFEM :: computeVolumeAroundLayer(GaussPoint *gp, int layer)
{
    double detJ;
    FloatMatrix Gcov;
    FloatArray lcoords;
    lcoords = gp->giveNaturalCoordinates();
    this->evalInitialCovarBaseVectorsAt(lcoords, Gcov);
    detJ = Gcov.giveDeterminant() * 0.5 * this->layeredCS->giveLayerThickness(layer);
    return detJ *gp->giveWeight();
}



} // end namespace oofem

