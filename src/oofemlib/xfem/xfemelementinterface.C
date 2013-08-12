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


    // TODO: Maybe we can allow for several enrichment items to
    std :: vector< FloatMatrix > Bd(nDofMan); // One Bd per node

    int counter = nDofMan * dim;
    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
        EnrichmentItem *ei = xMan->giveEnrichmentItem(i);


        double levelSetGP = 0.0;
        ei->interpLevelSet(levelSetGP, N, elNodes);

        FloatArray gradLevelSetGP;
        ei->interpGradLevelSet(gradLevelSetGP, dNdx, elNodes);



        for ( int j = 1; j <= nDofMan; j++ ) {
            DofManager *dMan = iEl.giveDofManager(j);

            // TODO: Make sure that this checks for enrichments of both bulk and tip
            if ( ei->isDofManEnriched(* dMan) ) {
                // Different nodes can be enriched by different enrichment functions.
                // Therefore, we have to compute enrichment functions inside this loop.

                int numEnr = ei->giveNumDofManEnrichments(* dMan);

                FloatMatrix &BdNode = Bd [ j - 1 ];
                BdNode.resize(3, numEnr * dim);
                BdNode.zero();


                int globalNodeInd = dMan->giveGlobalNumber();


                // Enrichment function derivative in Gauss point
                std :: vector< FloatArray >efgpD;
                ei->evaluateEnrFuncDerivAt(efgpD, globalCoord, levelSetGP, gradLevelSetGP, globalNodeInd);
                //ei->evaluateEnrFuncDerivAt(efgpD, (*iGP.giveCoordinates()), levelSetGP, gradLevelSetGP, globalNodeInd);

                // Enrichment function in Gauss Point
                std :: vector< double >efGP;
                ei->evaluateEnrFuncAt(efGP, globalCoord, levelSetGP, globalNodeInd);
                //ei->evaluateEnrFuncAt(efGP, *(iGP.giveCoordinates()), levelSetGP, globalNodeInd);


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

					BdNode.at(1, 2*k+1) = grad_ef_N.at(1);
					BdNode.at(2, 2*k+2) = grad_ef_N.at(2);
					BdNode.at(3, 2*k+1) = grad_ef_N.at(2);
					BdNode.at(3, 2*k+2) = grad_ef_N.at(1);

//					BdNode.at(1, 1) = grad_ef_N.at(1);
//					BdNode.at(2, 2) = grad_ef_N.at(2);
//					BdNode.at(3, 1) = grad_ef_N.at(2);
//					BdNode.at(3, 2) = grad_ef_N.at(1);

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

void XfemElementInterface :: XfemElementInterface_partitionElement(AList< Triangle > *answer, std :: vector< FloatArray > &together)
{
    Delaunay dl;
    dl.triangulate(together, answer);
}

void XfemElementInterface :: XfemElementInterface_updateIntegrationRule()
{
    XfemManager *xMan = this->element->giveDomain()->giveXfemManager();
    if ( xMan->isElementEnriched(element) ) {
        AList< Triangle >triangles;
        AList< Triangle >triangles2;
        // all the points coming into triangulation

        std :: vector< std :: vector< FloatArray > >pointPartitions;

        this->XfemElementInterface_prepareNodesForDelaunay(pointPartitions);

        for ( int i = 0; i < int( pointPartitions.size() ); i++ ) {
            this->XfemElementInterface_partitionElement(& triangles, pointPartitions [ i ]);
        }




        int elIndex = element->giveGlobalNumber();

        ///////////////////////////////////////////
        // Write splitted elements to vtk.
        std :: stringstream str1;
        str1 << "TriEl" << elIndex << "Side1.vtk";
        std :: string name1 = str1.str();

        XFEMDebugTools :: WriteTrianglesToVTK(name1, triangles);

        std :: stringstream str2;
        str2 << "TriEl" << elIndex << "Side2.vtk";
        std :: string name2 = str2.str();

        XFEMDebugTools :: WriteTrianglesToVTK(name2, triangles2);

        std :: vector< Triangle >allTri;


        //        for ( int i = 1; i <= triangles2.giveSize(); i++ ) {
        //            int sz = triangles.giveSize();
        //            triangles.put( sz + 1, triangles2.at(i) );
        //            triangles2.unlink(i);
        //        }


        //        for(int i = 1; i <= triangles.giveSize(); i++) {
        //            allTri.push_back( *(triangles.at(i)) );
        //        }

        for ( int i = 1; i <= triangles2.giveSize(); i++ ) {
            int sz = triangles.giveSize();
            //                triangles.put( sz + 1, triangles2.at(i) );
            //                triangles2.unlink(i);

            triangles.put( sz + 1, new Triangle( * ( triangles2.at(i) ) ) );
        }



        std :: stringstream str3;
        str3 << "TriEl" << elIndex << ".vtk";
        std :: string name3 = str3.str();

        XFEMDebugTools :: WriteTrianglesToVTK(name3, triangles);


        for ( int i = 1; i <= triangles.giveSize(); i++ ) {
            Triangle t( *( triangles.at(i) ) );
            allTri.push_back(t);
        }

        int ruleNum = 1;
        //        int numGPPerTri = 3;
        int numGPPerTri = 12;
        //        int numGPPerTri = 25;
        AList< IntegrationRule >irlist;
        IntegrationRule *intRule = new PatchIntegrationRule(ruleNum, element, allTri);

        MaterialMode matMode = element->giveMaterialMode();
        intRule->SetUpPointsOnTriangle(numGPPerTri, matMode);

        irlist.put(1, intRule);
        element->setIntegrationRules(& irlist);

/*
        AList< IntegrationRule >irlist;
        
        for ( int i = 1; i <= triangles.giveSize(); i++ ) {

            Patch *patch = new TrianglePatch(element);
            for ( int j = 1; j <= triangles.at(i)->giveVertices()->giveSize(); j++ ) {
                FloatArray *nCopy = new FloatArray( *triangles.at( i )->giveVertex(j) );
                patch->setVertex(nCopy);
            }

            PatchIntegrationRule *pir = new PatchIntegrationRule(i, element, patch);
            int pointNr = 3;
            MaterialMode matMode = element->giveMaterialMode();
            pir->SetUpPointsOnTriangle(pointNr, matMode);
            irlist.put(i, pir);
        }

        element->setIntegrationRules(& irlist);
*/
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

        // here the intersection points are copied in order to be put into two groups
        for ( int i = 1; i <= int(intersecPoints.size()); i++ )
        {

        	oPointPartitions[0].push_back(intersecPoints[i-1]);
//        	int sz = answer1->giveSize();
//            answer1->put( sz + 1, new FloatArray(intersecPoints[i-1]) );

            FloatArray *ip = &intersecPoints[i-1];
//            FloatArray *ipCopy = new FloatArray(*ip);

            oPointPartitions[1].push_back(*ip);
//            int sz2 = answer2->giveSize();
//            answer2->put(sz2 + 1, ipCopy);
        }

        // here the group is determined
        double x1 = intersecPoints [ 0 ].at(1);
        double x2 = intersecPoints [ 1 ].at(1);
        double y1 = intersecPoints [ 0 ].at(2);
        double y2 = intersecPoints [ 1 ].at(2);
        for ( int i = 1; i <= this->element->giveNumberOfDofManagers(); i++ ) {
            double x = element->giveDofManager(i)->giveCoordinates()->at(1);
            double y = element->giveDofManager(i)->giveCoordinates()->at(2);
            double det = ( x1 - x ) * ( y2 - y ) - ( x2 - x ) * ( y1 - y );
            FloatArray *node = element->giveDofManager(i)->giveCoordinates();
//            FloatArray *nodesCopy = new FloatArray(*node);
            if ( det > 0.00001 )
            {
            	oPointPartitions[0].push_back(*node);
//                int sz = answer1->giveSize();
//                answer1->put(sz + 1, nodesCopy);
            } else if ( det < ( -1 ) * 0.00001 )
            {
            	oPointPartitions[1].push_back(*node);
//                int sz = answer2->giveSize();
//                answer2->put(sz + 1, nodesCopy);
            }
        }
    }
    else if( intersecPoints.size() == 1 )
    {
    	// TODO: For now, assume that the number of edges is
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
/*
				// Take edge center or intersection point
				oPointPartitions[i-1].push_back(edgeCoords[i-1]);

				// Take next node
				if(i == nNodes)
				{
					oPointPartitions[i-1].push_back( *(element->giveDofManager(1)->giveCoordinates()) );
				}
				else
				{
					oPointPartitions[i-1].push_back( *(element->giveDofManager(i+1)->giveCoordinates()) );
				}

				// Take center on next edge
				if(i == nNodes)
				{
					oPointPartitions[i-1].push_back(edgeCoords[0]);
				}
				else
				{
					oPointPartitions[i-1].push_back(edgeCoords[i]);
				}

				// Take crack tip position
				oPointPartitions[i-1].push_back(tipCoord);
*/

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
//                double x = element->giveDofManager(i)->giveCoordinates()->at(1);
//                double y = element->giveDofManager(i)->giveCoordinates()->at(2);

                FloatArray *node = element->giveDofManager(i)->giveCoordinates();
                FloatArray *nodesCopy = new FloatArray(*node);

                oPointPartitions[0].push_back(*nodesCopy);
            }
        }
    }

    // nodes of an element are copied to a different memory location
    // so that the whole container of points for triangulation can be dealt with
    // more easily (e.g. deleted)

    //    for ( int i = 1; i <= intersecPoints.giveSize(); i++ ) {
    //        intersecPoints.unlink(i);
    //    }
}
} // end namespace oofem
