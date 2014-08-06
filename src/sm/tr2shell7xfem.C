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
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "fei3dtrquad.h"
#include "boundaryload.h"
#include "classfactory.h"
#include "tr2shell7.h"
#include "xfem/patchintegrationrule.h"
#include "xfem/XFEMDebugTools.h"
#include "xfem/enrichmentdomain.h"
#include "xfem/enrichmentitems/crack.h"
#include "xfem/enrichmentitems/shellcrack.h"
#include <string>
#include <sstream>


namespace oofem {

REGISTER_Element( Tr2Shell7XFEM );

FEI3dTrQuad Tr2Shell7XFEM :: interpolation;

IntArray Tr2Shell7XFEM :: ordering_all = {1, 2, 3, 8, 9, 10, 15, 16, 17, 22, 23, 24, 29, 30, 31, 36, 37, 38,
                        4, 5, 6, 11, 12, 13, 18, 19, 20, 25, 26, 27, 32, 33, 34, 39, 40, 41,
                        7, 14, 21, 28, 35, 42};
IntArray Tr2Shell7XFEM :: ordering_gr {1, 2, 3, 19, 20, 21, 37, 4, 5, 6, 22, 23, 24, 38, 7, 8, 9, 25, 26, 27, 39,
                       10, 11, 12, 28, 29, 30, 40, 13, 14, 15, 31, 32, 33, 41, 16, 17, 18,
                       34, 35, 36, 42};
IntArray Tr2Shell7XFEM :: ordering_gr_edge = {1, 2, 3, 10, 11, 12, 19, 4, 5, 6, 13, 14, 15, 20, 7, 8, 9, 16, 17, 18, 21};


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
    nodeLocalXiCoords = {1., 0., 0., .5, 0., .5};  // corner nodes then midnodes, uncertain of node numbering
    nodeLocalEtaCoords = {0., 1., 0., .5, .5, 0.};
}


FEInterpolation* Tr2Shell7XFEM :: giveInterpolation() const { return &interpolation; }


void 
Tr2Shell7XFEM :: computeGaussPoints()
{

    this->xMan = this->giveDomain()->giveXfemManager();

    if ( integrationRulesArray.size() == 0 ) {  
        if( this->xMan->isElementEnriched(this) ) {
            //this->updateIntegrationRule();
            this->updateIntegrationRuleMultiCrack();
        }
    }
    
    if ( integrationRulesArray.size() == 0 ) {

        int nPointsTri  = 6;   // points in the plane
        int nPointsEdge = 2;   // edge integration            

        // Cohesive zone
        for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { 
            //Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) ); 
            //if (dei) {
                int numberOfInterfaces = this->layeredCS->giveNumberOfLayers()-1;
                czIntegrationRulesArray.resize( numberOfInterfaces );
                for ( int i = 0; i < numberOfInterfaces; i++ ) {
                    czIntegrationRulesArray [ i ] = new GaussIntegrationRule(1, this);
                    czIntegrationRulesArray [ i ]->SetUpPointsOnTriangle(nPointsTri, _3dInterface);
                }
            //}
        }

        // Layered cross section for bulk integration
        //this->numberOfIntegrationRules = this->layeredCS->giveNumberOfLayers();
        this->numberOfGaussPoints = this->layeredCS->giveNumberOfLayers()*nPointsTri*this->layeredCS->giveNumIntegrationPointsInLayer();
        this->layeredCS->setupLayeredIntegrationRule(integrationRulesArray, this, nPointsTri);
        
        specialIntegrationRulesArray.resize(3);

        // Midplane (Mass matrix integrated analytically through the thickness)
        specialIntegrationRulesArray [ 1 ] = new GaussIntegrationRule(1, this);
        specialIntegrationRulesArray [ 1 ]->SetUpPointsOnTriangle(nPointsTri, _3dMat); 

        // Edge
        specialIntegrationRulesArray [ 2 ] = new GaussIntegrationRule(1, this);
        specialIntegrationRulesArray [ 2 ]->SetUpPointsOnLine(nPointsEdge, _3dMat);
        
        // Thickness integration for stress recovery
        specialIntegrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this);
        specialIntegrationRulesArray [ 0 ]->SetUpPointsOnLine(this->layeredCS->giveNumIntegrationPointsInLayer(), _3dMat);

    }

}


bool Tr2Shell7XFEM :: updateIntegrationRule()
{
    bool partitionSucceeded = false;

    if ( this->xMan->isElementEnriched(this) ) {


        std :: vector< std :: vector< FloatArray > >pointPartitions;
        //std :: vector< Triangle > allTri;

        int numEI = this->xMan->giveNumberOfEnrichmentItems();
        for ( int eiIndex = 1; eiIndex <= numEI; eiIndex++ ) {
            EnrichmentItem *ei = this->xMan->giveEnrichmentItem(eiIndex);
            if ( dynamic_cast< Crack*> (ei) ) {

                // Get the points describing each subdivision of the element
                double startXi, endXi;
                bool intersection = false;
                this->XfemElementInterface_prepareNodesForDelaunay(pointPartitions, startXi, endXi, eiIndex, intersection);

                if ( intersection ) {
                    // Use XfemElementInterface_partitionElement to subdivide the element
                    for ( int i = 0; i < int( pointPartitions.size() ); i++ ) {
                        // Triangulate the subdivisions
                        this->XfemElementInterface_partitionElement(this->allTri, pointPartitions [ i ]);
                    }

                    partitionSucceeded = true;
                }
            
            }
        }

        ////////////////////////////////////////
        // When we reach this point, we have a triangulation that is adapted to all
        // cracks passing through the element. Therefore, we can set up integration
        // points on each triangle.

        if ( this->xMan->giveVtkDebug() ) {
            std :: stringstream str3;
            int elIndex = this->giveGlobalNumber();
            str3 << "TriEl" << elIndex << ".vtk";
            std :: string name3 = str3.str();

            XFEMDebugTools :: WriteTrianglesToVTK(name3, allTri);
            XFEMDebugTools :: WriteTrianglesToVTK(name3, this->allTri);
        }
        
        // Create integrationrule based on a 'wedge patch'
        int nPointsTri  = 6;   // points in the plane

        if ( this->allTri.size() == 0 ) { // No subdivision, return and create iRule as normal
            return partitionSucceeded;

        } else { // create iRule according to subdivision
            int numberOfLayers     = this->layeredCS->giveNumberOfLayers();
            int numPointsThickness = this->layeredCS->giveNumIntegrationPointsInLayer();

            //integrationRulesArray = new IntegrationRule * [ numberOfLayers ];
            integrationRulesArray.resize(numberOfLayers);
            for ( int i = 0; i < numberOfLayers; i++ ) {
                integrationRulesArray [ i ] = new PatchIntegrationRule(1, this, this->allTri);
                integrationRulesArray [ i ]->SetUpPointsOnWedge(nPointsTri, numPointsThickness, _3dMat);             
            }
            this->layeredCS->mapLayerGpCoordsToShellCoords(integrationRulesArray);
            
        }
    }

    int nPointsTri  = 6;   // points in the plane  

    // Cohesive zone
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { 
        //Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) ); 
            int numberOfInterfaces = this->layeredCS->giveNumberOfLayers()-1;
            //czIntegrationRulesArray = new IntegrationRule * [ numberOfInterfaces ];
            czIntegrationRulesArray.resize(numberOfInterfaces);
            for ( int i = 0; i < numberOfInterfaces; i++ ) {
                czIntegrationRulesArray [ i ] = new GaussIntegrationRule(1, this);
                czIntegrationRulesArray [ i ]->SetUpPointsOnTriangle(nPointsTri, _3dInterface);
            }
    }

    return partitionSucceeded;
}




bool Tr2Shell7XFEM :: updateIntegrationRuleMultiCrack()
{
    bool partitionSucceeded = false;
    int nPointsTri  = 6;   // points in the plane

    //if ( this->xMan->isElementEnriched(this) ) {
    if ( 0 ) {

        int numberOfLayers     = this->layeredCS->giveNumberOfLayers();
        int numPointsThickness = this->layeredCS->giveNumIntegrationPointsInLayer();
        double totalThickness  = this->layeredCS->computeIntegralThick();
        int numEI = this->xMan->giveNumberOfEnrichmentItems();
        std :: vector< std :: vector< FloatArray > >pointPartitions;


        //integrationRulesArray = new IntegrationRule * [ numberOfLayers ];
        integrationRulesArray.resize(numberOfLayers);
        this->crackSubdivisions.resize(numberOfLayers);
        for ( int i = 0; i < numberOfLayers; i++ ) {
            double zMid_i = this->layeredCS->giveLayerMidZ(i+1); // global z-coord
            double xiMid_i = 1.0 - 2.0 * ( totalThickness - this->layeredCS->giveMidSurfaceZcoordFromBottom() - zMid_i ) / totalThickness; // local z-coord

            for ( int eiIndex = 1; eiIndex <= numEI; eiIndex++ ) {
                
                EnrichmentItem *ei = this->xMan->giveEnrichmentItem(eiIndex);
                if ( dynamic_cast< ShellCrack*> (ei) ) {

                    // Determine if the crack goes through the current layer
                    if( this->evaluateHeavisideGamma(xiMid_i, ei) > 0) {

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

                            partitionSucceeded = true;
                        }

                        integrationRulesArray [ i ] = new PatchIntegrationRule(1, this, this->crackSubdivisions [ i ]);
                        integrationRulesArray [ i ]->SetUpPointsOnWedge(nPointsTri, numPointsThickness, _3dMat);             
                    }
                }
            }                 
        }
    
        this->layeredCS->mapLayerGpCoordsToShellCoords(integrationRulesArray);
    }
        
    // Cohesive zone
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { 
        Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) ); 
            int numberOfInterfaces = this->layeredCS->giveNumberOfLayers()-1;
            //czIntegrationRulesArray = new IntegrationRule * [ numberOfInterfaces ];
            czIntegrationRulesArray.resize(numberOfInterfaces);
            for ( int i = 0; i < numberOfInterfaces; i++ ) {
                czIntegrationRulesArray [ i ] = new GaussIntegrationRule(1, this);
                czIntegrationRulesArray [ i ]->SetUpPointsOnTriangle(nPointsTri, _3dInterface);
            }
    }

    return partitionSucceeded;
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
    lcoords = * gp->giveNaturalCoordinates();
    this->evalInitialCovarBaseVectorsAt(lcoords, Gcov);
    detJ = Gcov.giveDeterminant() * 0.5 * this->layeredCS->giveLayerThickness(layer);
    return detJ *gp->giveWeight();
}



} // end namespace oofem

