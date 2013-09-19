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

#include "xfemelementinterface.h"
#include "enrichmentitem.h"
#include "engngm.h"
#include "gausspoint.h"
#include "materialmode.h"
#include "fei2dquadlin.h"
//#include "patch.h"
#include "patchintegrationrule.h"
#include "delaunay.h"
#include "xfemmanager.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "enrichmentdomain.h"

#include "XFEMDebugTools.h"
#include <string>
#include <sstream>

namespace oofem {
void XfemElementInterface :: XfemElementInterface_createEnrBmatrixAt(FloatMatrix &oAnswer, GaussPoint &iGP, Element &iEl)
{
    const int dim = 2;
    const int nDofMan = iEl.giveNumberOfDofManagers();

    FloatMatrix dNdx;
    FloatArray N;
    FEInterpolation *interp = iEl.giveInterpolation();
    interp->evaldNdx( dNdx, * iGP.giveCoordinates(), FEIElementGeometryWrapper(& iEl) );
    interp->evalN( N, * iGP.giveCoordinates(), FEIElementGeometryWrapper(& iEl) );

    const IntArray &elNodes = iEl.giveDofManArray();

    // Compute global coordinates of Gauss point
    FloatArray globalCoord;
    globalCoord.setValues(2, 0.0, 0.0);

    for ( int i = 1; i <= nDofMan; i++ ) {
        DofManager *dMan = iEl.giveDofManager(i);
        globalCoord.at(1) += N.at(i) * dMan->giveCoordinate(1);
        globalCoord.at(2) += N.at(i) * dMan->giveCoordinate(2);
    }


    // Standard FE part of B-matrix
    std :: vector< FloatMatrix > Bc(nDofMan);
    for ( int i = 1; i <= nDofMan; i++ ) {
        FloatMatrix &BNode = Bc [ i - 1 ];
        BNode.resize(3, 2);
        BNode.zero();
        BNode.at(1, 1) = dNdx.at(i, 1);
        BNode.at(2, 2) = dNdx.at(i, 2);
        BNode.at(3, 1) = dNdx.at(i, 2);
        BNode.at(3, 2) = dNdx.at(i, 1);
    }


    // XFEM part of B-matrix
    XfemManager *xMan = iEl.giveDomain()->giveXfemManager();


    std :: vector< FloatMatrix > Bd(nDofMan); // One Bd per node

    int counter = nDofMan * dim;

    for ( int j = 1; j <= nDofMan; j++ ) {

        DofManager *dMan = iEl.giveDofManager(j);

    	// Compute the total number of enrichments for node j
    	int numEnrNode = 0;
        for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
            EnrichmentItem *ei = xMan->giveEnrichmentItem(i);
            if ( ei->isDofManEnriched(* dMan) ) {
            	numEnrNode += ei->giveNumDofManEnrichments(* dMan);
            }
        }

        FloatMatrix &BdNode = Bd [ j - 1 ];
        BdNode.resize(3, numEnrNode * dim);
        BdNode.zero();


        int globalNodeInd = dMan->giveGlobalNumber();

        int nodeEnrCounter = 0;

        for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
        	EnrichmentItem *ei = xMan->giveEnrichmentItem(i);

        	double levelSetGP = 0.0;
        	ei->interpLevelSet(levelSetGP, N, elNodes);

        	FloatArray gradLevelSetGP;
        	ei->interpGradLevelSet(gradLevelSetGP, dNdx, elNodes);


        	if ( ei->isDofManEnriched(* dMan) ) {

                int numEnr = ei->giveNumDofManEnrichments(* dMan);

                // Enrichment function derivative in Gauss point
                std :: vector< FloatArray >efgpD;
                ei->evaluateEnrFuncDerivAt(efgpD, globalCoord, levelSetGP, gradLevelSetGP, globalNodeInd);

                // Enrichment function in Gauss Point
                std :: vector< double >efGP;
                ei->evaluateEnrFuncAt(efGP, globalCoord, levelSetGP, globalNodeInd);


                const FloatArray &nodePos = * ( dMan->giveCoordinates() );

                double levelSetNode  = 0.0;
                ei->evalLevelSetNormalInNode( levelSetNode, dMan->giveGlobalNumber() );

                std :: vector< double >efNode;
                ei->evaluateEnrFuncAt(efNode, nodePos, levelSetNode, globalNodeInd);


				for(int k = 0; k < numEnr; k++) {
					// matrix to be added anytime a node is enriched
					// Creates nabla*(ef*N)
					FloatArray grad_ef_N;
					grad_ef_N.resize(dim);
					for ( int p = 1; p <= dim; p++ ) {
						grad_ef_N.at(p) = dNdx.at(j, p) * ( efGP[k] - efNode[k] ) + N.at(j) * efgpD[k].at(p);
					}

					BdNode.at(1, nodeEnrCounter+1) = grad_ef_N.at(1);
					BdNode.at(2, nodeEnrCounter+2) = grad_ef_N.at(2);
					BdNode.at(3, nodeEnrCounter+1) = grad_ef_N.at(2);
					BdNode.at(3, nodeEnrCounter+2) = grad_ef_N.at(1);

					nodeEnrCounter += 2;
					counter += 2;
				}
			}

		}
    }


    // Create the total B-matrix by appending each contribution to B after one another.
    oAnswer.resize(3, counter);
    oAnswer.zero();
    int column = 1;
    for ( int i = 0; i < nDofMan; i++ ) {
        oAnswer.setSubMatrix(Bc [ i ], 1, column);
        column += 2;
        if ( Bd [ i ].isNotEmpty() ) {
            oAnswer.setSubMatrix(Bd [ i ], 1, column);

            column += Bd [ i ].giveNumberOfColumns();
        }
    }
}

void XfemElementInterface :: XfemElementInterface_partitionElement(std::vector< Triangle > &oTriangles, const std :: vector< FloatArray > &iPoints)
{
    Delaunay dl;
    dl.triangulate(iPoints, oTriangles);
}

void XfemElementInterface :: XfemElementInterface_updateIntegrationRule()
{
    XfemManager *xMan = this->element->giveDomain()->giveXfemManager();
    if ( xMan->isElementEnriched(element) ) {

        std :: vector< std :: vector< FloatArray > >pointPartitions;
        std :: vector< Triangle >allTri;

        // Get the points describing each subdivision of the element
        this->XfemElementInterface_prepareNodesForDelaunay(pointPartitions);

        for ( int i = 0; i < int( pointPartitions.size() ); i++ ) {
        	// Triangulate the subdivisions
            this->XfemElementInterface_partitionElement( allTri, pointPartitions [ i ]);
        }

#if XFEM_DEBUG_VTK > 0
        std :: stringstream str3;
        int elIndex = this->element->giveGlobalNumber();
        str3 << "TriEl" << elIndex << ".vtk";
        std :: string name3 = str3.str();

        XFEMDebugTools :: WriteTrianglesToVTK(name3, allTri);
#endif


        int ruleNum = 1;
        AList< IntegrationRule >irlist;
        IntegrationRule *intRule = new PatchIntegrationRule(ruleNum, element, allTri);

        MaterialMode matMode = element->giveMaterialMode();
        intRule->SetUpPointsOnTriangle( xMan->giveNumGpPerTri() , matMode);

        irlist.put(1, intRule);
        element->setIntegrationRules(& irlist);
    }
}

void XfemElementInterface :: XfemElementInterface_prepareNodesForDelaunay(std :: vector< std :: vector< FloatArray > > &oPointPartitions)
{
    XfemManager *xMan = this->element->giveDomain()->giveXfemManager();

    std :: vector< FloatArray >intersecPoints;
    std :: vector< int >intersecEdgeInd;

    // TODO:    Can we do this recursively to achieve proper splitting
    //			when several enrichment items interact with the
    //			same element?
    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
        xMan->giveEnrichmentItem(i)->computeIntersectionPoints(intersecPoints, intersecEdgeInd, element);
    }


    if ( intersecPoints.size() == 2 ) {
        // The element is completely cut in two.
        // Therefore, we create two subpartitions:
        // one on each side of the interface.
        oPointPartitions.resize(2);

        for ( int i = 1; i <= int(intersecPoints.size()); i++ ) {
        	oPointPartitions[0].push_back(intersecPoints[i-1]);
            oPointPartitions[1].push_back(intersecPoints[i-1]);
        }


        // Check on which side of the interface each node is located.
        const double &x1 = intersecPoints [ 0 ].at(1);
        const double &x2 = intersecPoints [ 1 ].at(1);
        const double &y1 = intersecPoints [ 0 ].at(2);
        const double &y2 = intersecPoints [ 1 ].at(2);

        for ( int i = 1; i <= this->element->giveNumberOfDofManagers(); i++ ) {
        	const double &x = element->giveDofManager(i)->giveCoordinates()->at(1);
            const double &y = element->giveDofManager(i)->giveCoordinates()->at(2);
            double det = ( x1 - x ) * ( y2 - y ) - ( x2 - x ) * ( y1 - y );
            FloatArray *node = element->giveDofManager(i)->giveCoordinates();

            if ( det > 0.0 ) {
            	oPointPartitions[0].push_back(*node);
            } else {
            	oPointPartitions[1].push_back(*node);
            }
        }
    }
    else if( intersecPoints.size() == 1 )
    {
    	// TODO: For now, assume that the number of element edges is
    	// equal to the number of nodes.
    	int nNodes = this->element->giveNumberOfNodes();
    	std::vector<FloatArray> edgeCoords, nodeCoords;

        FloatArray tipCoord;
        int dim = element->giveDofManager(1)->giveCoordinates()->giveSize();
        tipCoord.resize(dim);

        bool foundTip = false;
        for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
            EnrichmentItem *ei = xMan->giveEnrichmentItem(i);

            if ( ei->giveElementTipCoord( tipCoord, element->giveNumber() ) ) {
                foundTip = true;
                break;
            }
        }

        if(foundTip)
        {

			for(int i = 1; i <= nNodes; i++)
			{
				// Store edge points
				if( i == intersecEdgeInd[0] )
				{
					// Take the intersection point ...
					edgeCoords.push_back(intersecPoints[0]);
				}
				else
				{
					// ... or the center of the edge.
					IntArray bNodes;
					this->element->giveInterpolation()->boundaryGiveNodes(bNodes, i);

					int nsLoc = bNodes.at( 1 );
					int neLoc = bNodes.at( bNodes.giveSize() );

					const FloatArray &coordS = *(element->giveDofManager(nsLoc)->giveCoordinates());
					const FloatArray &coordE = *(element->giveDofManager(neLoc)->giveCoordinates());

					FloatArray coordEdge;
					coordEdge = 0.5*coordS + 0.5*coordE;
					edgeCoords.push_back(coordEdge);
				}

				// Store node coords
				const FloatArray &coord = *(element->giveDofManager(i)->giveCoordinates());
				nodeCoords.push_back( coord );

			}

			oPointPartitions.resize((2*nNodes));

			// Divide into subdomains
			for(int i = 1; i <= nNodes; i++)
			{
				////////////////
				// Take edge center or intersection point
				oPointPartitions[2*i-1].push_back(edgeCoords[i-1]);

				// Take crack tip position
				oPointPartitions[2*i-1].push_back(tipCoord);

				// Take node
				oPointPartitions[2*i-1].push_back( *(element->giveDofManager(i)->giveCoordinates()) );

				////////////////
				// Take edge center or intersection point
				oPointPartitions[2*i-2].push_back(edgeCoords[i-1]);

				// Take next node
				if(i == nNodes)
				{
					oPointPartitions[2*i-2].push_back( *(element->giveDofManager(1)->giveCoordinates()) );
				}
				else
				{
					oPointPartitions[2*i-2].push_back( *(element->giveDofManager(i+1)->giveCoordinates()) );
				}

				// Take crack tip position
				oPointPartitions[2*i-2].push_back(tipCoord);

			}
        } // If a tip was found
        else
        {
        	oPointPartitions.resize(1);

            for ( int i = 1; i <= this->element->giveNumberOfDofManagers(); i++ )
            {
                const FloatArray &nodeCoord = *element->giveDofManager(i)->giveCoordinates();
                oPointPartitions[0].push_back(nodeCoord);
            }
        }
    }

}

void XfemElementInterface :: recomputeGaussPoints() {

	bool recompute = false;


	// Do checks to determine if the Gauss points need to be recomputed
	// For now, we choose to always recompute cut elements.
    XfemManager *xMan = element->giveDomain()->giveXfemManager();

    for(int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++) {
    	std::vector<FloatArray> intersecPoints;
    	EnrichmentItem *ei = xMan->giveEnrichmentItem(i);

        std::vector< int > intersecEdgeInd;
    	ei->computeIntersectionPoints(intersecPoints, intersecEdgeInd, element);
    	int numIntersecPoints = intersecPoints.size();

        if ( numIntersecPoints > 0 )
        {
        	recompute = true;
        }

    }


	if( recompute ) {

		// Fetch old Gauss points


		// Create new partitioning (and delete old Gauss points)

        this->XfemElementInterface_updateIntegrationRule();

		// Map Gauss point variables
		// (area weighted least squares?)

	}

}

} // end namespace oofem
