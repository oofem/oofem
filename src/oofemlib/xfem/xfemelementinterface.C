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

#include "xfemelementinterface.h"
#include "enrichmentitem.h"
#include "geometrybasedei.h"
#include "engngm.h"
#include "gausspoint.h"
#include "materialmode.h"
#include "fei2dquadlin.h"
#include "patchintegrationrule.h"
#include "delaunay.h"
#include "xfemmanager.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "dynamicinputrecord.h"
#include "mathfem.h"

#include "XFEMDebugTools.h"
#include <string>
#include <sstream>

namespace oofem {
XfemElementInterface :: XfemElementInterface(Element *e) :
    Interface(),
    element(e),
    mUsePlaneStrain(false)
{
    mpCZIntegrationRules.clear();
}

XfemElementInterface :: ~XfemElementInterface()
{}

void XfemElementInterface :: XfemElementInterface_createEnrBmatrixAt(FloatMatrix &oAnswer, GaussPoint &iGP, Element &iEl)
{
    ComputeBOrBHMatrix(oAnswer, iGP, iEl, false);
}

void XfemElementInterface :: XfemElementInterface_createEnrBHmatrixAt(FloatMatrix &oAnswer, GaussPoint &iGP, Element &iEl)
{
    ComputeBOrBHMatrix(oAnswer, iGP, iEl, true);
}

void XfemElementInterface :: ComputeBOrBHMatrix(FloatMatrix &oAnswer, GaussPoint &iGP, Element &iEl, bool iComputeBH)
{
    /*
     * Computes the B or BH matrix.
     * iComputeBH = true implies that BH is computed,
     * while B is computed if iComputeBH = false.
     */
    const int dim = 2;
    const int nDofMan = iEl.giveNumberOfDofManagers();

    int shearInd = 3, numRows = 3;
    if ( mUsePlaneStrain ) {
        shearInd = 4;
        numRows = 4;
    }

    if ( iComputeBH ) {
        numRows++;
    }

    FloatMatrix dNdx;
    FloatArray N;
    FEInterpolation *interp = iEl.giveInterpolation();
    const FEIElementGeometryWrapper geomWrapper(& iEl);
    interp->evaldNdx(dNdx, iGP.giveNaturalCoordinates(), geomWrapper);
    interp->evalN(N, iGP.giveNaturalCoordinates(), geomWrapper);

    const IntArray &elNodes = iEl.giveDofManArray();

    // Compute global coordinates of Gauss point
    FloatArray globalCoord = {
        0.0, 0.0
    };

    for ( int i = 1; i <= nDofMan; i++ ) {
        const Node *node = iEl.giveNode(i);
        const FloatArray &nodeCoord = node->giveNodeCoordinates();
        globalCoord.at(1) += N.at(i) * nodeCoord [ 0 ];
        globalCoord.at(2) += N.at(i) * nodeCoord [ 1 ];
    }


    // Standard FE part of B-matrix
    std :: vector< FloatMatrix >Bc(nDofMan);
    for ( int i = 1; i <= nDofMan; i++ ) {
        FloatMatrix &BNode = Bc [ i - 1 ];
        BNode.resize(numRows, 2);
        BNode.at(1, 1)                  = dNdx.at(i, 1);
        BNode.at(2, 2)                  = dNdx.at(i, 2);
        BNode.at(shearInd, 1)   = dNdx.at(i, 2);

        if ( iComputeBH ) {
            BNode.at(shearInd + 1, 2)   = dNdx.at(i, 1);
        } else {
            BNode.at(shearInd, 2)   = dNdx.at(i, 1);
        }
    }


    // XFEM part of B-matrix
    double enrDofsScaleFactor = 1.0;
    XfemManager *xMan = NULL;
    if ( iEl.giveDomain()->hasXfemManager() ) {
        xMan = iEl.giveDomain()->giveXfemManager();
        enrDofsScaleFactor = xMan->giveEnrDofScaleFactor();
    }

    std :: vector< FloatMatrix >Bd(nDofMan);   // One Bd per node

    int counter = nDofMan * dim;

    int numEnrNode = 0;

    for ( int j = 1; j <= nDofMan; j++ ) {
        DofManager *dMan = iEl.giveDofManager(j);
        const Node *node = iEl.giveNode(j);

        // Compute the total number of enrichments for node j
        if ( iEl.giveDomain()->hasXfemManager() ) {
            numEnrNode = XfemElementInterface_giveNumDofManEnrichments(* dMan, * xMan);
        }

        if ( numEnrNode > 0 ) {
            FloatMatrix &BdNode = Bd [ j - 1 ];
            BdNode.resize(numRows, numEnrNode * dim);


            const int globalNodeInd = dMan->giveGlobalNumber();

            int nodeEnrCounter = 0;

            const std :: vector< int > &nodeEiIndices = xMan->giveNodeEnrichmentItemIndices(globalNodeInd);
            for ( size_t i = 0; i < nodeEiIndices.size(); i++ ) {
                EnrichmentItem *ei = xMan->giveEnrichmentItem(nodeEiIndices [ i ]);

                if ( ei->isDofManEnriched(* dMan) ) {
                    int numEnr = ei->giveNumDofManEnrichments(* dMan);

                    // Enrichment function derivative in Gauss point
                    std :: vector< FloatArray >efgpD;
                    ei->evaluateEnrFuncDerivAt(efgpD, globalCoord, iGP.giveNaturalCoordinates(), globalNodeInd, * element, N, dNdx, elNodes);
                    // Enrichment function in Gauss Point
                    std :: vector< double >efGP;
                    ei->evaluateEnrFuncAt(efGP, globalCoord, iGP.giveNaturalCoordinates(), globalNodeInd, * element, N, elNodes);


                    const FloatArray &nodePos = node->giveNodeCoordinates();

                    double levelSetNode  = 0.0;
                    ei->evalLevelSetNormalInNode(levelSetNode, globalNodeInd, nodePos);

                    std :: vector< double >efNode;
                    FloatArray nodeNaturalCoord;
                    iEl.computeLocalCoordinates(nodeNaturalCoord, nodePos);
                    ei->evaluateEnrFuncInNode(efNode, * node);

                    for ( int k = 0; k < numEnr; k++ ) {
                        // matrix to be added anytime a node is enriched
                        // Creates nabla*(ef*N)
                        FloatArray grad_ef_N;
                        grad_ef_N.resize(dim);
                        for ( int p = 1; p <= dim; p++ ) {
                            grad_ef_N.at(p) = dNdx.at(j, p) * ( efGP [ k ] - efNode [ k ] ) + N.at(j) * efgpD [ k ].at(p);
                        }

                        BdNode.at(1, nodeEnrCounter + 1)                  = grad_ef_N.at(1);
                        BdNode.at(2, nodeEnrCounter + 2)                  = grad_ef_N.at(2);
                        BdNode.at(shearInd, nodeEnrCounter + 1)   = grad_ef_N.at(2);

                        if ( iComputeBH ) {
                            BdNode.at(shearInd + 1, nodeEnrCounter + 2)   = grad_ef_N.at(1);
                        } else {
                            BdNode.at(shearInd, nodeEnrCounter + 2)   = grad_ef_N.at(1);
                        }

                        nodeEnrCounter += 2;
                        counter += 2;
                    }
                }
            }
        }
    }


    // Create the total B-matrix by appending each contribution to B after one another.
    oAnswer.resize(numRows, counter);

    int column = 1;
    for ( int i = 0; i < nDofMan; i++ ) {
        oAnswer.setSubMatrix(Bc [ i ], 1, column);
        column += 2;
        if ( Bd [ i ].isNotEmpty() ) {
            Bd[i].times(enrDofsScaleFactor);
            oAnswer.setSubMatrix(Bd [ i ], 1, column);

            column += Bd [ i ].giveNumberOfColumns();
        }
    }
}

void XfemElementInterface :: XfemElementInterface_createEnrNmatrixAt(FloatMatrix &oAnswer, const FloatArray &iLocCoord, Element &iEl, bool iSetDiscontContribToZero)
{
    std :: vector< int >elNodes;

    int numElNodes = iEl.giveNumberOfDofManagers();

    for ( int i = 0; i < numElNodes; i++ ) {
        elNodes.push_back(i + 1);
    }

    XfemElementInterface_createEnrNmatrixAt(oAnswer, iLocCoord, iEl, elNodes, iSetDiscontContribToZero);
}

void XfemElementInterface :: XfemElementInterface_createEnrNmatrixAt(FloatMatrix &oAnswer, const FloatArray &iLocCoord, Element &iEl, const std :: vector< int > &iLocNodeInd, bool iSetDiscontContribToZero)
{
    const int dim = 2;
    const int nDofMan = iEl.giveNumberOfDofManagers();

    FloatArray Nc;
    FEInterpolation *interp = iEl.giveInterpolation();
    interp->evalN( Nc, iLocCoord, FEIElementGeometryWrapper(& iEl) );

    const IntArray &elNodes = iEl.giveDofManArray();

    // Compute global coordinates of Gauss point
    FloatArray globalCoord(2);
    globalCoord.zero();

    for ( int i = 1; i <= nDofMan; i++ ) {
        DofManager *dMan = iEl.giveDofManager(i);
        globalCoord.at(1) += Nc.at(i) * dMan->giveCoordinate(1);
        globalCoord.at(2) += Nc.at(i) * dMan->giveCoordinate(2);
    }


    // XFEM part of N-matrix
    XfemManager *xMan = iEl.giveDomain()->giveXfemManager();
    double enrDofsScaleFactor = xMan->giveEnrDofScaleFactor();

    int counter = iLocNodeInd.size() * dim;

    std :: vector< std :: vector< double > >Nd( iLocNodeInd.size() );
    for ( int j = 1; j <= int( iLocNodeInd.size() ); j++ ) {
        DofManager *dMan = iEl.giveDofManager(iLocNodeInd [ j - 1 ]);
        Node *node = dynamic_cast< Node * >( dMan );

        // Compute the total number of enrichments for node j
        int numEnrNode = XfemElementInterface_giveNumDofManEnrichments(* dMan, * xMan);
        std :: vector< double > &NdNode = Nd [ j - 1 ];
        NdNode.assign(numEnrNode, 0.0);


        int globalNodeInd = dMan->giveGlobalNumber();

        size_t nodeCounter = 0;

        int placeInArray = element->giveDomain()->giveDofManPlaceInArray(globalNodeInd);
        const std :: vector< int > &nodeEiIndices = xMan->giveNodeEnrichmentItemIndices(placeInArray);
        for ( size_t i = 0; i < nodeEiIndices.size(); i++ ) {
            EnrichmentItem *ei = xMan->giveEnrichmentItem(nodeEiIndices [ i ]);

            if ( ei->isDofManEnriched(* dMan) ) {
                int numEnr = ei->giveNumDofManEnrichments(* dMan);


                // Enrichment function in Gauss Point
                std :: vector< double >efGP;
                ei->evaluateEnrFuncAt(efGP, globalCoord, iLocCoord, globalNodeInd, iEl, Nc, elNodes);


                const FloatArray &nodePos = * ( dMan->giveCoordinates() );

                std :: vector< double >efNode;

                FloatArray nodePosLocCoord;
                iEl.computeLocalCoordinates(nodePosLocCoord, nodePos);
                ei->evaluateEnrFuncInNode(efNode, * node);


                for ( int k = 0; k < numEnr; k++ ) {
                    if ( iSetDiscontContribToZero ) {
                        NdNode [ nodeCounter ] = 0.0;
                    } else   {
                        NdNode [ nodeCounter ] = ( efGP [ k ] - efNode [ k ] ) * Nc.at(j);
                    }
                    counter++;
                    nodeCounter++;
                }
            }
        }
    }

    int numN = iLocNodeInd.size();

    for ( int j = 1; j <= int( iLocNodeInd.size() ); j++ ) {
        numN += Nd [ j - 1 ].size();
    }

    FloatArray NTot;
    NTot.resize(numN);
    NTot.zero();
    int column = 1;

    for ( int i = 1; i <= int( iLocNodeInd.size() ); i++ ) {
        NTot.at(column) = Nc.at(iLocNodeInd [ i - 1 ]);
        column++;

        const std :: vector< double > &NdNode = Nd [ i - 1 ];
        for ( size_t j = 1; j <= NdNode.size(); j++ ) {
            NTot.at(column) = NdNode [ j - 1 ]*enrDofsScaleFactor;
            column++;
        }
    }

    oAnswer.beNMatrixOf(NTot, 2);
}

int XfemElementInterface :: XfemElementInterface_giveNumDofManEnrichments(const DofManager &iDMan, XfemManager &iXMan) const
{
    int numEnrNode = 0;
    int globalNodeInd = iDMan.giveGlobalNumber();
    int placeInArray = element->giveDomain()->giveDofManPlaceInArray(globalNodeInd);
    const std :: vector< int > &nodeEiIndices = iXMan.giveNodeEnrichmentItemIndices(placeInArray);
    for ( size_t i = 0; i < nodeEiIndices.size(); i++ ) {
        EnrichmentItem *ei = iXMan.giveEnrichmentItem(nodeEiIndices [ i ]);
        if ( ei->isDofManEnriched(iDMan) ) {
            numEnrNode += ei->giveNumDofManEnrichments(iDMan);
        }
    }

    return numEnrNode;
}

void XfemElementInterface :: XfemElementInterface_partitionElement(std :: vector< Triangle > &oTriangles, const std :: vector< FloatArray > &iPoints)
{
    Delaunay dl;
    dl.triangulate(iPoints, oTriangles);
}



bool XfemElementInterface :: XfemElementInterface_updateIntegrationRule()
{
    bool partitionSucceeded = false;

    XfemManager *xMan = this->element->giveDomain()->giveXfemManager();
    if ( xMan->isElementEnriched(element) ) {
        MaterialMode matMode = element->giveMaterialMode();

        bool firstIntersection = true;

        std :: vector< std :: vector< FloatArray > >pointPartitions;
        std :: vector< Triangle >allTri;

        std :: vector< int >enrichingEIs;
        int elPlaceInArray = xMan->giveDomain()->giveElementPlaceInArray( element->giveGlobalNumber() );
        xMan->giveElementEnrichmentItemIndices(enrichingEIs, elPlaceInArray);


        for ( size_t p = 0; p < enrichingEIs.size(); p++ ) {
            int eiIndex = enrichingEIs [ p ];

            if ( firstIntersection ) {
                // Get the points describing each subdivision of the element
                double startXi, endXi;
                bool intersection = false;
                this->XfemElementInterface_prepareNodesForDelaunay(pointPartitions, startXi, endXi, eiIndex, intersection);

                if ( intersection ) {
                    firstIntersection = false;

                    // Use XfemElementInterface_partitionElement to subdivide the element
                    for ( int i = 0; i < int ( pointPartitions.size() ); i++ ) {
                        // Triangulate the subdivisions
                        this->XfemElementInterface_partitionElement(allTri, pointPartitions [ i ]);
                    }

                    partitionSucceeded = true;
                }
            } // if(firstIntersection)
            else {
                // Loop over triangles
                std :: vector< Triangle >allTriCopy;
                for ( size_t triIndex = 0; triIndex < allTri.size(); triIndex++ ) {
                    // Call alternative version of XfemElementInterface_prepareNodesForDelaunay
                    std :: vector< std :: vector< FloatArray > >pointPartitionsTri;
                    double startXi, endXi;
                    bool intersection = false;
                    XfemElementInterface_prepareNodesForDelaunay(pointPartitionsTri, startXi, endXi, allTri [ triIndex ], eiIndex, intersection);

                    if ( intersection ) {
                        // Use XfemElementInterface_partitionElement to subdivide triangle j
                        for ( int i = 0; i < int ( pointPartitionsTri.size() ); i++ ) {
                            this->XfemElementInterface_partitionElement(allTriCopy, pointPartitionsTri [ i ]);
                        }
                    } else {
                        allTriCopy.push_back(allTri [ triIndex ]);
                    }
                }

                allTri = allTriCopy;
            }
        }

        ////////////////////////////////////////
        // When we reach this point, we have a
        // triangulation that is adapted to all
        // cracks passing through the element.
        // Therefore, we can set up integration
        // points on each triangle.

        if ( xMan->giveVtkDebug() ) {
            std :: stringstream str3;
            int elIndex = this->element->giveGlobalNumber();
            str3 << "TriEl" << elIndex << ".vtk";
            std :: string name3 = str3.str();

            if ( allTri.size() > 0 ) {
                XFEMDebugTools :: WriteTrianglesToVTK(name3, allTri);
            }
        }


        int ruleNum = 1;
        if ( partitionSucceeded ) {
            std :: vector< std :: unique_ptr< IntegrationRule > >intRule;
            intRule.emplace_back( new PatchIntegrationRule(ruleNum, element, allTri) );
            intRule [ 0 ]->SetUpPointsOnTriangle(xMan->giveNumGpPerTri(), matMode);
            element->setIntegrationRules( std :: move(intRule) );
        }
    }

    return partitionSucceeded;
}

void XfemElementInterface :: XfemElementInterface_prepareNodesForDelaunay(std :: vector< std :: vector< FloatArray > > &oPointPartitions, double &oCrackStartXi, double &oCrackEndXi, int iEnrItemIndex, bool &oIntersection)
{
    int dim = element->giveDofManager(1)->giveCoordinates()->giveSize();

    FloatArray elCenter( element->giveDofManager(1)->giveCoordinates()->giveSize() );
    elCenter.zero();
    std :: vector< const FloatArray * >nodeCoord;
    for ( int i = 1; i <= this->element->giveNumberOfDofManagers(); i++ ) {
        nodeCoord.push_back( element->giveDofManager(i)->giveCoordinates() );
        elCenter.add( * ( element->giveDofManager(i)->giveCoordinates() ) );
    }
    elCenter.times( 1.0 / double( element->giveNumberOfDofManagers() ) );

    XfemManager *xMan = this->element->giveDomain()->giveXfemManager();
    GeometryBasedEI *ei = dynamic_cast< GeometryBasedEI * >( xMan->giveEnrichmentItem(iEnrItemIndex) );

    if ( ei == NULL ) {
        oIntersection = false;
        return;
    }

    std :: vector< FloatArray >intersecPoints;
    std :: vector< int >intersecEdgeInd;

    std :: vector< double >minDistArcPos;
    ei->computeIntersectionPoints(intersecPoints, intersecEdgeInd, element, minDistArcPos);

    for ( size_t i = 0; i < intersecPoints.size(); i++ ) {
        intersecPoints [ i ].resizeWithValues(dim);
    }


    if ( intersecPoints.size() == 2 ) {
        // The element is completely cut in two.
        // Therefore, we create two subpartitions:
        // one on each side of the interface.
        oPointPartitions.resize(2);

        putPointsInCorrectPartition(oPointPartitions, intersecPoints, nodeCoord);

        // Export start and end points of
        // the intersection line.
        oCrackStartXi   = std :: min(minDistArcPos [ 0 ], minDistArcPos [ 1 ]);
        oCrackEndXi     = std :: max(minDistArcPos [ 0 ], minDistArcPos [ 1 ]);

        oIntersection = true;
        return;
    } else if ( intersecPoints.size() == 1 ) {
        std :: vector< FloatArray >edgeCoords;

        FloatArray tipCoord;
        tipCoord.resize(dim);

        bool foundTip = false;
        double tipArcPos = -1.0;

        if ( ei->giveElementTipCoord(tipCoord, tipArcPos, * element, elCenter) ) {
            foundTip = true;
            tipCoord.resizeWithValues(dim);
        }
        int nEdges = this->element->giveInterpolation()->giveNumberOfEdges();
        if ( foundTip ) {
            oPointPartitions.resize( ( nEdges + 1 ) );

            // Divide into subdomains
            int triPassed = 0;
            for ( int i = 1; i <= nEdges; i++ ) {
                IntArray bNodes;
                this->element->giveInterpolation()->boundaryGiveNodes(bNodes, i);
                int nsLoc = bNodes.at(1);
                int neLoc = bNodes.at( bNodes.giveSize() );

                const FloatArray &coordS = * ( element->giveDofManager(nsLoc)->giveCoordinates() );
                const FloatArray &coordE = * ( element->giveDofManager(neLoc)->giveCoordinates() );

                if ( i == intersecEdgeInd [ 0 ] ) {
                    oPointPartitions [ triPassed ].push_back(tipCoord);
                    oPointPartitions [ triPassed ].push_back(intersecPoints [ 0 ]);
                    oPointPartitions [ triPassed ].push_back(coordE);
                    triPassed++;

                    oPointPartitions [ triPassed ].push_back(tipCoord);
                    oPointPartitions [ triPassed ].push_back(coordS);
                    oPointPartitions [ triPassed ].push_back(intersecPoints [ 0 ]);
                    triPassed++;
                } else   {
                    oPointPartitions [ triPassed ].push_back(tipCoord);
                    oPointPartitions [ triPassed ].push_back(coordS);
                    oPointPartitions [ triPassed ].push_back(coordE);
                    triPassed++;
                }
            }

            // Export start and end points of
            // the intersection line.
            oCrackStartXi   = std :: min(minDistArcPos [ 0 ], tipArcPos);
            oCrackEndXi     = std :: max(minDistArcPos [ 0 ], tipArcPos);
        }             // If a tip was found
        else {
            //            printf( "Warning: no tip found in element %d with only one edge intersection.\n", element->giveGlobalNumber() );

            oPointPartitions.resize(1);

            for ( int i = 1; i <= this->element->giveNumberOfDofManagers(); i++ ) {
                const FloatArray &nodeCoord = * element->giveDofManager(i)->giveCoordinates();
                oPointPartitions [ 0 ].push_back(nodeCoord);
            }

            // test Jim
            // Add first intersection point
            oPointPartitions [ 0 ].push_back(intersecPoints [ 0 ]);

            // want to add the extrapolated intersection point
            //FloatArray test;
            //test.setValues(2, 0.0, 0.4);
            //oPointPartitions [ 0 ].push_back( test );

            // Export start and end points of
            // the intersection line.
            oCrackStartXi   = minDistArcPos [ 0 ];
            oCrackEndXi = tipArcPos;
        }

        oIntersection = true;


        //oPointPartitions.resize(0);
        //oIntersection = true;
        //false;
        return;
    }

    oIntersection = false;
}

void XfemElementInterface :: XfemElementInterface_prepareNodesForDelaunay(std :: vector< std :: vector< FloatArray > > &oPointPartitions, double &oCrackStartXi, double &oCrackEndXi, const Triangle &iTri, int iEnrItemIndex, bool &oIntersection)
{
    int dim = element->giveDofManager(1)->giveCoordinates()->giveSize();

    FloatArray elCenter( iTri.giveVertex(1).giveSize() );
    elCenter.zero();
    std :: vector< const FloatArray * >nodeCoord;
    for ( int i = 1; i <= 3; i++ ) {
        nodeCoord.push_back( & iTri.giveVertex(i) );
        elCenter.add( iTri.giveVertex(i) );
    }
    elCenter.times(1.0 / 3.0);


    XfemManager *xMan = this->element->giveDomain()->giveXfemManager();
    GeometryBasedEI *ei = dynamic_cast< GeometryBasedEI * >( xMan->giveEnrichmentItem(iEnrItemIndex) );

    if ( ei == NULL ) {
        oIntersection = false;
        return;
    }

    std :: vector< FloatArray >intersecPoints;
    std :: vector< int >intersecEdgeInd;

    std :: vector< double >minDistArcPos;
    ei->computeIntersectionPoints(intersecPoints, intersecEdgeInd, element, iTri, minDistArcPos);
    for ( size_t i = 0; i < intersecPoints.size(); i++ ) {
        intersecPoints [ i ].resizeWithValues(dim);
    }

    if ( intersecPoints.size() == 2 ) {
        // The element is completely cut in two.
        // Therefore, we create two subpartitions:
        // one on each side of the interface.

        oPointPartitions.resize(2);

        putPointsInCorrectPartition(oPointPartitions, intersecPoints, nodeCoord);

        // Export start and end points of
        // the intersection line.
        oCrackStartXi   = std :: min(minDistArcPos [ 0 ], minDistArcPos [ 1 ]);
        oCrackEndXi     = std :: max(minDistArcPos [ 0 ], minDistArcPos [ 1 ]);

        oIntersection = true;
        return;
    } else if ( intersecPoints.size() == 1 ) {
        int nEdges = 3;
        std :: vector< FloatArray >edgeCoords, nodeCoords;

        FloatArray tipCoord;
        int dim = element->giveDofManager(1)->giveCoordinates()->giveSize();
        tipCoord.resize(dim);

        bool foundTip = false;
        double tipArcPos = -1.0;

        if ( ei->giveElementTipCoord(tipCoord, tipArcPos, * element, elCenter) ) {
            tipCoord.resizeWithValues(dim);
            foundTip = true;
        }

        if ( foundTip ) {
            oPointPartitions.resize( ( nEdges + 1 ) );

            // Divide into subdomains
            int triPassed = 0;
            for ( int i = 1; i <= nEdges; i++ ) {
                const FloatArray &coordS = iTri.giveVertex(i);

                int endInd = i + 1;
                if ( i == nEdges ) {
                    endInd = 1;
                }
                const FloatArray &coordE = iTri.giveVertex(endInd);

                if ( i == intersecEdgeInd [ 0 ] ) {
                    oPointPartitions [ triPassed ].push_back(tipCoord);
                    oPointPartitions [ triPassed ].push_back(intersecPoints [ 0 ]);
                    oPointPartitions [ triPassed ].push_back(coordE);
                    triPassed++;

                    oPointPartitions [ triPassed ].push_back(tipCoord);
                    oPointPartitions [ triPassed ].push_back(coordS);
                    oPointPartitions [ triPassed ].push_back(intersecPoints [ 0 ]);
                    triPassed++;
                } else   {
                    oPointPartitions [ triPassed ].push_back(tipCoord);
                    oPointPartitions [ triPassed ].push_back(coordS);
                    oPointPartitions [ triPassed ].push_back(coordE);
                    triPassed++;
                }
            }

            // Export start and end points of
            // the intersection line.
            oCrackStartXi   = std :: min(minDistArcPos [ 0 ], tipArcPos);
            oCrackEndXi     = std :: max(minDistArcPos [ 0 ], tipArcPos);
        }             // If a tip was found
        else {
            oPointPartitions.resize(1);

            for ( int i = 1; i <= 3; i++ ) {
                const FloatArray &nodeCoord = iTri.giveVertex(i);
                oPointPartitions [ 0 ].push_back(nodeCoord);
            }

            // Export start and end points of
            // the intersection line.
            oCrackStartXi   = minDistArcPos [ 0 ];
            oCrackEndXi = tipArcPos;
        }

        oIntersection = true;
        return;
    }

    oIntersection = false;
}

void XfemElementInterface :: putPointsInCorrectPartition(std :: vector< std :: vector< FloatArray > > &oPointPartitions,
                                                         const std :: vector< FloatArray > &iIntersecPoints,
                                                         const std :: vector< const FloatArray * > &iNodeCoord) const
{
    for ( size_t i = 0; i < iIntersecPoints.size(); i++ ) {
        oPointPartitions [ 0 ].push_back(iIntersecPoints [ i ]);
        oPointPartitions [ 1 ].push_back(iIntersecPoints [ i ]);
    }

    // Check on which side of the interface each node is located.
    const double &x1 = iIntersecPoints [ 0 ].at(1);
    const double &x2 = iIntersecPoints [ 1 ].at(1);
    const double &y1 = iIntersecPoints [ 0 ].at(2);
    const double &y2 = iIntersecPoints [ 1 ].at(2);

    for ( size_t i = 1; i <= iNodeCoord.size(); i++ ) {
        const double &x = iNodeCoord [ i - 1 ]->at(1);
        const double &y = iNodeCoord [ i - 1 ]->at(2);
        double det = ( x1 - x ) * ( y2 - y ) - ( x2 - x ) * ( y1 - y );

        if ( det > 0.0 ) {
            oPointPartitions [ 0 ].push_back(* iNodeCoord [ i - 1 ]);
        } else {
            oPointPartitions [ 1 ].push_back(* iNodeCoord [ i - 1 ]);
        }
    }
}

void XfemElementInterface :: partitionEdgeSegment(int iBndIndex, std :: vector< Line > &oSegments, std :: vector< FloatArray > &oIntersectionPoints, const double &iTangDistPadding)
{
    const double levelSetTol2 = 1.0e-12;
    //    const double gammaPadding = 0.001;

    XfemManager *xMan = this->element->giveDomain()->giveXfemManager();

    FEInterpolation *interp = element->giveInterpolation(); // Geometry interpolation
    IntArray edgeNodes;
    FEInterpolation2d *interp2d = dynamic_cast< FEInterpolation2d * >( interp );
    if ( interp2d == NULL ) {
        OOFEM_ERROR("In XfemElementInterface :: partitionEdgeSegment: failed to cast to FEInterpolation2d.\n")
    }
    interp2d->computeLocalEdgeMapping(edgeNodes, iBndIndex);

    // Fetch start and end points.
    const FloatArray &xS = * ( element->giveDofManager( edgeNodes.at(1) )->giveCoordinates() );
    const FloatArray &xE = * ( element->giveDofManager( edgeNodes.at(2) )->giveCoordinates() );

    // The point of departure is the original edge segment.
    // This segment will be subdivided as many times as necessary.
    Line seg1(xS, xE);
    //    oSegments.clear();
    oSegments.push_back(seg1);


    // Loop over enrichment items
    int numEI = xMan->giveNumberOfEnrichmentItems();
    for ( int eiIndex = 1; eiIndex <= numEI; eiIndex++ ) {
        EnrichmentItem *ei = xMan->giveEnrichmentItem(eiIndex);

        std :: vector< Line >newSegments;

        // Loop over segments
        size_t numSeg = oSegments.size();
        for ( size_t segInd = 0; segInd < numSeg; segInd++ ) {
            // Check if the segment is cut by the current enrichment item

            const FloatArray &seg_xS = oSegments [ segInd ].giveVertex(1);
            const FloatArray &seg_xE = oSegments [ segInd ].giveVertex(2);


            // Local coordinates of vertices
            FloatArray xiS;
            bool evaluationSucceeded = true;
            if ( !element->computeLocalCoordinates(xiS, seg_xS) ) {
                //TODO: Check for numerical round-off error
                //                printf("xiS: "); xiS.printYourself();
                //                evaluationSucceeded = false;
            }
            FloatArray xiE;
            if ( !element->computeLocalCoordinates(xiE, seg_xE) ) {
                //                printf("xiE: "); xiE.printYourself();
                //                evaluationSucceeded = false;
            }

            const IntArray &elNodes = element->giveDofManArray();
            FloatArray Ns, Ne;
            interp->evalN( Ns, xiS, FEIElementGeometryWrapper(element) );
            interp->evalN( Ne, xiE, FEIElementGeometryWrapper(element) );

            double phiS         = 0.0, phiE     = 0.0;
            double gammaS       = 0.0, gammaE   = 0.0;


            for ( int i = 1; i <= Ns.giveSize(); i++ ) {
                const FloatArray &nodePos = * ( element->giveNode(i)->giveCoordinates() );
                double phiNode = 0.0;
                if ( !ei->evalLevelSetNormalInNode(phiNode, elNodes [ i - 1 ], nodePos) ) {
                    evaluationSucceeded = false;
                }

                double gammaNode = 0.0;
                if ( !ei->evalLevelSetTangInNode(gammaNode, elNodes [ i - 1 ], nodePos) ) {
                    evaluationSucceeded = false;
                }

                phiS += Ns.at(i) * phiNode;
                gammaS += Ns.at(i) * gammaNode;

                phiE += Ne.at(i) * phiNode;
                gammaE += Ne.at(i) * gammaNode;
            }


            if ( phiS * phiE < levelSetTol2 && evaluationSucceeded ) {
                double xi = EnrichmentItem :: calcXiZeroLevel(phiS, phiE);
                double gamma = 0.5 * ( 1.0 - xi ) * gammaS + 0.5 * ( 1.0 + xi ) * gammaE;

                // If we are inside in tangential direction
                if ( gamma > -iTangDistPadding ) {
                    // If so, subdivide it ...

                    // Compute global coordinates of the intersection point
                    int nDim = std :: min( seg_xS.giveSize(), seg_xE.giveSize() );
                    FloatArray p;
                    p.resize(nDim);

                    for ( int i = 1; i <= nDim; i++ ) {
                        ( p.at(i) ) = 0.5 * ( 1.0 - xi ) * ( ( seg_xS.at(i) ) ) + 0.5 * ( 1.0 + xi ) * ( ( seg_xE.at(i) ) );
                    }

                    Line segA(seg_xS, p);
                    newSegments.push_back(segA);
                    Line segB(p, seg_xE);
                    newSegments.push_back(segB);

                    // Export the intersection point
                    oIntersectionPoints.push_back(p);
                } else {
                    newSegments.push_back(oSegments [ segInd ]);
                }
            } else {
                // ... else keep the segment.
                newSegments.push_back(oSegments [ segInd ]);
            }
        }

        oSegments = newSegments;
    }
}

MaterialMode XfemElementInterface :: giveMaterialMode()
{
    if ( mUsePlaneStrain ) {
        return _PlaneStrain;
    } else {
        return _PlaneStress;
    }
}

void XfemElementInterface :: updateYourselfCZ(TimeStep *tStep)
{
    size_t numSeg = mpCZIntegrationRules.size();

    for ( size_t i = 0; i < numSeg; i++ ) {
        if ( mpCZIntegrationRules [ i ] != NULL ) {
            mpCZIntegrationRules [ i ]->updateYourself(tStep);
        }
    }
}

void XfemElementInterface :: computeDisplacementJump(oofem :: GaussPoint &iGP, oofem :: FloatArray &oJump, const oofem :: FloatArray &iSolVec, const oofem :: FloatMatrix &iNMatrix)
{
    const int dim = 2;
    oJump.resize(dim);
    oJump.beProductOf(iNMatrix, iSolVec);
}

void XfemElementInterface :: computeNCohesive(FloatMatrix &oN, GaussPoint &iGP, int iEnrItemIndex, const std :: vector< int > &iTouchingEnrItemIndices)
{
    const int dim = 2;

    FloatArray Nc, globalCoord, localCoord;
    globalCoord = iGP.giveGlobalCoordinates();
    element->computeLocalCoordinates(localCoord, globalCoord);
    FEInterpolation *interp = element->giveInterpolation();
    interp->evalN( Nc, localCoord, FEIElementGeometryWrapper(element) );

    const int nDofMan = element->giveNumberOfDofManagers();

    // XFEM part of N-matrix
    XfemManager *xMan = element->giveDomain()->giveXfemManager();
    double enrDofsScaleFactor = xMan->giveEnrDofScaleFactor();


    int counter = nDofMan * dim;

    std :: vector< std :: vector< double > >Nd(nDofMan);

    for ( int j = 1; j <= nDofMan; j++ ) {
        DofManager *dMan = element->giveDofManager(j);

        // Compute the total number of enrichments for node j
        int numEnrNode = XfemElementInterface_giveNumDofManEnrichments(* dMan, * xMan);
        std :: vector< double > &NdNode = Nd [ j - 1 ];
        NdNode.assign(numEnrNode, 0.0);


        int globalNodeInd = dMan->giveGlobalNumber();


        int ndNodeInd = 0;
        const std :: vector< int > &nodeEiIndices = xMan->giveNodeEnrichmentItemIndices(globalNodeInd);
        for ( size_t i = 0; i < nodeEiIndices.size(); i++ ) {
            EnrichmentItem *ei = xMan->giveEnrichmentItem(nodeEiIndices [ i ]);

            GeometryBasedEI *geoEI = dynamic_cast< GeometryBasedEI * >( ei );

            if ( geoEI != NULL ) {
                if ( geoEI->isDofManEnriched(* dMan) ) {
                    int numEnr = geoEI->giveNumDofManEnrichments(* dMan);

                    std :: vector< double >efJumps;
                    bool gpLivesOnCurrentCrack = ( nodeEiIndices [ i ] == iEnrItemIndex );

                    bool gpLivesOnInteractingCrack = false;
                    for ( int touchingEIIndex : iTouchingEnrItemIndices ) {
                        if ( nodeEiIndices [ i ] == touchingEIIndex ) {
                            gpLivesOnInteractingCrack = true;
                        }
                    }

                    if ( nodeEiIndices [ i ] == iEnrItemIndex || gpLivesOnInteractingCrack ) {
                        geoEI->evaluateEnrFuncJumps(efJumps, globalNodeInd, iGP, gpLivesOnCurrentCrack);
                    }

                    for ( int k = 0; k < numEnr; k++ ) {
                        if ( nodeEiIndices [ i ] == iEnrItemIndex || gpLivesOnInteractingCrack ) {
                            NdNode [ ndNodeInd ] = efJumps [ k ] * Nc.at(j);
                        } else {
                            NdNode [ ndNodeInd ] = 0.0;
                        }

                        counter++;
                        ndNodeInd++;
                    }
                }
            }
        }
    }

    int numN = nDofMan;

    for ( int j = 1; j <= nDofMan; j++ ) {
        numN += Nd [ j - 1 ].size();
    }

    FloatArray NTot;
    NTot.resize(numN);
    NTot.zero();
    int column = 1;

    for ( int i = 1; i <= nDofMan; i++ ) {
        //        NTot.at(column) = Nc.at(i); // We do not want the continuous part.
        column++;

        const std :: vector< double > &NdNode = Nd [ i - 1 ];
        for ( size_t j = 1; j <= NdNode.size(); j++ ) {
            NTot.at(column) = NdNode [ j - 1 ]*enrDofsScaleFactor;
            column++;
        }
    }

    oN.beNMatrixOf(NTot, 2);
}
} // end namespace oofem
