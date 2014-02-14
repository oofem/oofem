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
#include "engngm.h"
#include "gausspoint.h"
#include "materialmode.h"
#include "fei2dquadlin.h"
#include "patchintegrationrule.h"
#include "delaunay.h"
#include "xfemmanager.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "enrichmentdomain.h"
#include "dynamicinputrecord.h"

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"
#include "structuralelement.h"

#include "XFEMDebugTools.h"
#include <string>
#include <sstream>
#include <math.h>

// TODO: Remove need for these includes in base class. /ES
#include "xfem/enrichmentitems/crack.h"

namespace oofem {
XfemElementInterface :: XfemElementInterface(Element *e) :
    Interface(),
    element(e),
    mpCZMat(NULL),
    mCZMaterialNum(-1),
    mCSNumGaussPoints(4),
    mUsePlaneStrain(false)
{
    mpCZIntegrationRules.clear();
}

XfemElementInterface :: ~XfemElementInterface()
{
    size_t numCZRules = mpCZIntegrationRules.size();

    for ( size_t i = 0; i < numCZRules; i++ ) {
        if ( mpCZIntegrationRules [ i ] != NULL ) {
            delete mpCZIntegrationRules [ i ];
            mpCZIntegrationRules [ i ] = NULL;
        }
    }

    mpCZIntegrationRules.clear();
}

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

    if(iComputeBH){
    	numRows++;
    }

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
        BNode.resize(numRows, 2);
        BNode.zero();
        BNode.at(1, 1)                  = dNdx.at(i, 1);
        BNode.at(2, 2)                  = dNdx.at(i, 2);
        BNode.at(shearInd, 1)   = dNdx.at(i, 2);

        if(iComputeBH){
        	BNode.at(shearInd+1	, 2)   = dNdx.at(i, 1);
        }
        else{
            BNode.at(shearInd	, 2)   = dNdx.at(i, 1);
        }
    }


    // XFEM part of B-matrix
    XfemManager *xMan = NULL;
	if( iEl.giveDomain()->hasXfemManager() ) {
		xMan = iEl.giveDomain()->giveXfemManager();
	}

    std :: vector< FloatMatrix > Bd(nDofMan);  // One Bd per node

    int counter = nDofMan * dim;

    for ( int j = 1; j <= nDofMan; j++ ) {
        DofManager *dMan = iEl.giveDofManager(j);

        // Compute the total number of enrichments for node j
        int numEnrNode = 0;

        if( xMan != NULL ) {
			for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
				EnrichmentItem *ei = xMan->giveEnrichmentItem(i);
				if ( ei->isDofManEnriched(* dMan) ) {
					numEnrNode += ei->giveNumDofManEnrichments(* dMan);
				}
			}
        }

        if ( numEnrNode > 0 ) {
            FloatMatrix &BdNode = Bd [ j - 1 ];
            BdNode.resize(numRows, numEnrNode * dim);
            BdNode.zero();


            int globalNodeInd = dMan->giveGlobalNumber();

            int nodeEnrCounter = 0;

            for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
                EnrichmentItem *ei = xMan->giveEnrichmentItem(i);

                double levelSetGP = 0.0;
                ei->interpLevelSet(levelSetGP, N, elNodes);

                FloatArray gradLevelSetGP(dim);
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

                        if(iComputeBH){
                        	BdNode.at(shearInd+1	, nodeEnrCounter + 2)   = grad_ef_N.at(1);
                        }
                        else {
                        	BdNode.at(shearInd		, nodeEnrCounter + 2)   = grad_ef_N.at(1);
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

void XfemElementInterface :: XfemElementInterface_createEnrNmatrixAt(FloatMatrix &oAnswer, const FloatArray &iLocCoord, Element &iEl)
{
    const int dim = 2;
    const int nDofMan = iEl.giveNumberOfDofManagers();

    FloatArray Nc;
    FEInterpolation *interp = iEl.giveInterpolation();
    interp->evalN( Nc, iLocCoord, FEIElementGeometryWrapper(& iEl) );

    const IntArray &elNodes = iEl.giveDofManArray();

    // Compute global coordinates of Gauss point
    FloatArray globalCoord;
    globalCoord.setValues(2, 0.0, 0.0);

    for ( int i = 1; i <= nDofMan; i++ ) {
        DofManager *dMan = iEl.giveDofManager(i);
        globalCoord.at(1) += Nc.at(i) * dMan->giveCoordinate(1);
        globalCoord.at(2) += Nc.at(i) * dMan->giveCoordinate(2);
    }


    // XFEM part of N-matrix
    XfemManager *xMan = iEl.giveDomain()->giveXfemManager();


    std :: vector< FloatMatrix > Bd(nDofMan);  // One Bd per node

    int counter = nDofMan * dim;

    std :: vector< std :: vector< double > > Nd(nDofMan);

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

        std :: vector< double > &NdNode = Nd [ j - 1 ];
        NdNode.assign(numEnrNode, 0.0);


        int globalNodeInd = dMan->giveGlobalNumber();

        size_t nodeCounter = 0;

        for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
            EnrichmentItem *ei = xMan->giveEnrichmentItem(i);

            double levelSetGP = 0.0;
            ei->interpLevelSet(levelSetGP, Nc, elNodes);


            if ( ei->isDofManEnriched(* dMan) ) {
                int numEnr = ei->giveNumDofManEnrichments(* dMan);


                // Enrichment function in Gauss Point
                std :: vector< double >efGP;
                ei->evaluateEnrFuncAt(efGP, globalCoord, levelSetGP, globalNodeInd);


                const FloatArray &nodePos = * ( dMan->giveCoordinates() );

                double levelSetNode  = 0.0;
                ei->evalLevelSetNormalInNode( levelSetNode, dMan->giveGlobalNumber() );

                std :: vector< double >efNode;
                ei->evaluateEnrFuncAt(efNode, nodePos, levelSetNode, globalNodeInd);


                for ( int k = 0; k < numEnr; k++ ) {
                    NdNode [ nodeCounter ] = ( efGP [ k ] - efNode [ k ] ) * Nc.at(j);
                    counter++;
                    nodeCounter++;
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
        NTot.at(column) = Nc.at(i);
        column++;

        const std :: vector< double > &NdNode = Nd [ i - 1 ];
        for ( size_t j = 1; j <= NdNode.size(); j++ ) {
            NTot.at(column) = NdNode [ j - 1 ];
            column++;
        }
    }

    oAnswer.beNMatrixOf(NTot, 2);
}

void XfemElementInterface :: XfemElementInterface_partitionElement(std :: vector< Triangle > &oTriangles, const std :: vector< FloatArray > &iPoints)
{
    Delaunay dl;
    dl.triangulate(iPoints, oTriangles);
}

bool XfemElementInterface :: XfemElementInterface_updateIntegrationRule()
{
    bool partitionSucceeded = false;


    if ( mpCZMat != NULL ) {
        for ( size_t i = 0; i < mpCZIntegrationRules.size(); i++ ) {
            if ( mpCZIntegrationRules [ i ] != NULL ) {
                delete mpCZIntegrationRules [ i ];
            }
        }

        mpCZIntegrationRules.clear();
        mCZEnrItemIndices.clear();
    }

    XfemManager *xMan = this->element->giveDomain()->giveXfemManager();
    if ( xMan->isElementEnriched(element) ) {
        if ( mpCZMat == NULL && mCZMaterialNum > 0 ) {
            initializeCZMaterial();
        }


        MaterialMode matMode = element->giveMaterialMode();

        bool firstIntersection = true;

        std :: vector< std :: vector< FloatArray > >pointPartitions;
        std :: vector< Triangle >allTri;

        int numEI = xMan->giveNumberOfEnrichmentItems();
        for ( int eiIndex = 1; eiIndex <= numEI; eiIndex++ ) {

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


                    if ( mpCZMat != NULL ) {
                        Crack *crack = dynamic_cast<Crack*>( xMan->giveEnrichmentItem(eiIndex) );
                        if(crack == NULL) {
                        	OOFEM_ERROR("Error in XfemElementInterface :: XfemElementInterface_updateIntegrationRule(): Cohesive zones are only available for cracks.\n")
                        }

                        // We have xi_s and xi_e. Fetch sub polygon.
                        std :: vector< FloatArray >crackPolygon;
                        crack->giveSubPolygon(crackPolygon, startXi, endXi);

                        /*
                         *                                              printf("crackPolygon: \n");
                         *                                              for(size_t i = 0; i < crackPolygon.size(); i++) {
                         *                                                      printf("i: %d x: %e y: %e\n", i, crackPolygon[i].at(1), crackPolygon[i].at(2) );
                         *                                              }
                         */

                        ///////////////////////////////////
                        // Add cohesive zone Gauss points
                        size_t numSeg = crackPolygon.size() - 1;

                        for ( size_t segIndex = 0; segIndex < numSeg; segIndex++ ) {
                            int czRuleNum = 1;
                            mpCZIntegrationRules.push_back( new GaussIntegrationRule(czRuleNum, element) );
                            mCZEnrItemIndices.push_back(eiIndex);
                            const FloatArray **coords = new const FloatArray * [ 2 ];
                            coords [ 0 ] = new FloatArray(crackPolygon [ segIndex      ]);
                            coords [ 1 ] = new FloatArray(crackPolygon [ segIndex + 1 ]);

                            // Compute crack normal
                            FloatArray crackTang;
                            crackTang.beDifferenceOf(crackPolygon [ segIndex + 1 ], crackPolygon [ segIndex         ]);
                            crackTang.normalize();

                            FloatArray crackNormal;
                            crackNormal.setValues( 2, -crackTang.at(2), crackTang.at(1) );

                            mpCZIntegrationRules [ segIndex ]->SetUpPointsOn2DEmbeddedLine(mCSNumGaussPoints, matMode, coords);

                            for ( int i = 0; i < mpCZIntegrationRules [ segIndex ]->giveNumberOfIntegrationPoints(); i++ ) {
                                double gw = mpCZIntegrationRules [ segIndex ]->getIntegrationPoint(i)->giveWeight();
                                double segLength = crackPolygon [ segIndex         ].distance(crackPolygon [ segIndex + 1 ]);
                                gw *= 0.5 * segLength;
                                GaussPoint &gp = * ( mpCZIntegrationRules [ segIndex ]->getIntegrationPoint(i) );
                                gp.setWeight(gw);

                                // Fetch material status and set normal
                                StructuralInterfaceMaterialStatus *ms = dynamic_cast< StructuralInterfaceMaterialStatus * >( mpCZMat->giveStatus(& gp) );
                                if ( ms == NULL ) {
                                    OOFEM_ERROR("In XfemElementInterface :: XfemElementInterface_updateIntegrationRule(): Failed to fetch material status.\n");
                                }

                                ms->letNormalBe(crackNormal);

                                // Give Gauss point reference to the enrichment item
                                // to simplify post processing.
                                crack->AppendCohesiveZoneGaussPoint(&gp);

                            }

                            delete coords [ 0 ];
                            delete coords [ 1 ];
                            delete [] coords;
                        }
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


                        // Add cohesive zone Gauss points

                        if ( mpCZMat != NULL ) {
                            Crack *crack = dynamic_cast<Crack*>( xMan->giveEnrichmentItem(eiIndex) );
                            if(crack == NULL) {
                            	OOFEM_ERROR("Error in XfemElementInterface :: XfemElementInterface_updateIntegrationRule(): Cohesive zones are only available for cracks.\n")
                            }

                            // We have xi_s and xi_e. Fetch sub polygon.
                            std :: vector< FloatArray >crackPolygon;
                            crack->giveSubPolygon(crackPolygon, startXi, endXi);

                            size_t numSeg = crackPolygon.size() - 1;

                            for ( size_t segIndex = 0; segIndex < numSeg; segIndex++ ) {
                                int czRuleNum = 1;
                                mpCZIntegrationRules.push_back( new GaussIntegrationRule(czRuleNum, element) );
                                size_t newRuleInd = mpCZIntegrationRules.size() - 1;
                                mCZEnrItemIndices.push_back(eiIndex);
                                const FloatArray **coords = new const FloatArray * [ 2 ];
                                coords [ 0 ] = new FloatArray(crackPolygon [ segIndex      ]);
                                coords [ 1 ] = new FloatArray(crackPolygon [ segIndex + 1 ]);

                                // Compute crack normal
                                FloatArray crackTang;
                                crackTang.beDifferenceOf(crackPolygon [ segIndex + 1 ], crackPolygon [ segIndex         ]);
                                crackTang.normalize();

                                FloatArray crackNormal;
                                crackNormal.setValues( 2, -crackTang.at(2), crackTang.at(1) );

                                mpCZIntegrationRules [ newRuleInd ]->SetUpPointsOn2DEmbeddedLine(mCSNumGaussPoints, matMode, coords);

                                for ( int i = 0; i < mpCZIntegrationRules [ newRuleInd ]->giveNumberOfIntegrationPoints(); i++ ) {
                                    double gw = mpCZIntegrationRules [ newRuleInd ]->getIntegrationPoint(i)->giveWeight();
                                    double segLength = crackPolygon [ segIndex         ].distance(crackPolygon [ segIndex + 1 ]);
                                    gw *= 0.5 * segLength;
                                    GaussPoint &gp = * ( mpCZIntegrationRules [ newRuleInd ]->getIntegrationPoint(i) );
                                    gp.setWeight(gw);

                                    // Fetch material status and set normal
                                    StructuralInterfaceMaterialStatus *ms = dynamic_cast< StructuralInterfaceMaterialStatus * >( mpCZMat->giveStatus(& gp) );
                                    if ( ms == NULL ) {
                                        OOFEM_ERROR("In XfemElementInterface :: XfemElementInterface_updateIntegrationRule(): Failed to fetch material status.\n");
                                    }

                                    ms->letNormalBe(crackNormal);

                                    // Give Gauss point reference to the enrichment item
                                    // to simplify post processing.
                                    crack->AppendCohesiveZoneGaussPoint(&gp);
                                }

                                delete coords [ 0 ];
                                delete coords [ 1 ];
                                delete [] coords;
                            }
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

            XFEMDebugTools :: WriteTrianglesToVTK(name3, allTri);
        }


        int ruleNum = 1;
        AList< IntegrationRule >irlist;
        IntegrationRule *intRule = new PatchIntegrationRule(ruleNum, element, allTri);

        intRule->SetUpPointsOnTriangle(xMan->giveNumGpPerTri(), matMode);

        irlist.put(1, intRule);
        if ( partitionSucceeded ) {
            element->setIntegrationRules(& irlist);
        }


        if ( xMan->giveVtkDebug() ) {
            ////////////////////////////////////////////////////////////////////////
            // Write CZ GP to VTK

            std :: vector< FloatArray >czGPCoord;

            for ( size_t czRulInd = 0; czRulInd < mpCZIntegrationRules.size(); czRulInd++ ) {
                for ( int i = 0; i < mpCZIntegrationRules [ czRulInd ]->giveNumberOfIntegrationPoints(); i++ ) {
                    czGPCoord.push_back( * ( mpCZIntegrationRules [ czRulInd ]->getIntegrationPoint(i)->giveCoordinates() ) );
                }
            }

            double time = 0.0;

            Element *el = element;

            Domain *dom = el->giveDomain();
            if ( dom != NULL ) {
                EngngModel *em = dom->giveEngngModel();
                if ( em != NULL ) {
                    TimeStep *ts = em->giveCurrentStep();
                    if ( ts != NULL ) {
                        time = ts->giveTargetTime();
                    }
                }
            }

            std :: stringstream str;
            int elIndex = this->element->giveGlobalNumber();
            str << "CZGaussPointsTime" << time << "El" << elIndex << ".vtk";
            std :: string name = str.str();

            XFEMDebugTools :: WritePointsToVTK(name, czGPCoord);
            ////////////////////////////////////////////////////////////////////////
        }
    }

    return partitionSucceeded;
}

void XfemElementInterface :: XfemElementInterface_prepareNodesForDelaunay(std :: vector< std :: vector< FloatArray > > &oPointPartitions, double &oCrackStartXi, double &oCrackEndXi, int iEnrItemIndex, bool &oIntersection)
{
    std :: vector< const FloatArray * >nodeCoord;
    for ( int i = 1; i <= this->element->giveNumberOfDofManagers(); i++ ) {
        nodeCoord.push_back( element->giveDofManager(i)->giveCoordinates() );
    }

    XfemManager *xMan = this->element->giveDomain()->giveXfemManager();
    EnrichmentItem *ei = xMan->giveEnrichmentItem(iEnrItemIndex);

    std :: vector< FloatArray >intersecPoints;
    std :: vector< int >intersecEdgeInd;

    std :: vector< double >minDistArcPos;
    ei->computeIntersectionPoints(intersecPoints, intersecEdgeInd, element, minDistArcPos);


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
        // TODO: For now, assume that the number of element edges is
        // equal to the number of nodes.
        int nNodes = this->element->giveNumberOfNodes();
        std :: vector< FloatArray >edgeCoords, nodeCoords;

        FloatArray tipCoord;
        int dim = element->giveDofManager(1)->giveCoordinates()->giveSize();
        tipCoord.resize(dim);

        bool foundTip = false;
        double tipArcPos = -1.0;

        if ( ei->giveElementTipCoord( tipCoord, tipArcPos, element->giveNumber() ) ) {
            foundTip = true;
        }

        if ( foundTip ) {
            for ( int i = 1; i <= nNodes; i++ ) {
                // Store edge points
                if ( i == intersecEdgeInd [ 0 ] ) {
                    // Take the intersection point ...
                    edgeCoords.push_back(intersecPoints [ 0 ]);
                } else {
                    // ... or the center of the edge.
                    IntArray bNodes;
                    this->element->giveInterpolation()->boundaryGiveNodes(bNodes, i);

                    int nsLoc = bNodes.at(1);
                    int neLoc = bNodes.at( bNodes.giveSize() );

                    const FloatArray &coordS = * ( element->giveDofManager(nsLoc)->giveCoordinates() );
                    const FloatArray &coordE = * ( element->giveDofManager(neLoc)->giveCoordinates() );

                    FloatArray coordEdge;
                    coordEdge = 0.5 * coordS + 0.5 * coordE;
                    edgeCoords.push_back(coordEdge);
                }

                // Store node coords
                const FloatArray &coord = * ( element->giveDofManager(i)->giveCoordinates() );
                nodeCoords.push_back(coord);
            }

            oPointPartitions.resize( ( 2 * nNodes ) );

            // Divide into subdomains
            for ( int i = 1; i <= nNodes; i++ ) {
                ////////////////
                // Take edge center or intersection point
                oPointPartitions [ 2 * i - 1 ].push_back(edgeCoords [ i - 1 ]);

                // Take crack tip position
                oPointPartitions [ 2 * i - 1 ].push_back(tipCoord);

                // Take node
                oPointPartitions [ 2 * i - 1 ].push_back( * ( element->giveDofManager(i)->giveCoordinates() ) );

                ////////////////
                // Take edge center or intersection point
                oPointPartitions [ 2 * i - 2 ].push_back(edgeCoords [ i - 1 ]);

                // Take next node
                if ( i == nNodes ) {
                    oPointPartitions [ 2 * i - 2 ].push_back( * ( element->giveDofManager(1)->giveCoordinates() ) );
                } else {
                    oPointPartitions [ 2 * i - 2 ].push_back( * ( element->giveDofManager(i + 1)->giveCoordinates() ) );
                }

                // Take crack tip position
                oPointPartitions [ 2 * i - 2 ].push_back(tipCoord);
            }

            // Export start and end points of
            // the intersection line.
            oCrackStartXi   = std :: min(minDistArcPos [ 0 ], tipArcPos);
            oCrackEndXi     = std :: max(minDistArcPos [ 0 ], tipArcPos);
        }             // If a tip was found
        else {
            printf( "Warning: no tip found in element %d with only one edge intersection.\n", element->giveGlobalNumber() );

            oPointPartitions.resize(1);

            for ( int i = 1; i <= this->element->giveNumberOfDofManagers(); i++ ) {
                const FloatArray &nodeCoord = * element->giveDofManager(i)->giveCoordinates();
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

void XfemElementInterface :: XfemElementInterface_prepareNodesForDelaunay(std :: vector< std :: vector< FloatArray > > &oPointPartitions, double &oCrackStartXi, double &oCrackEndXi, const Triangle &iTri, int iEnrItemIndex, bool &oIntersection)
{
    std :: vector< const FloatArray * >nodeCoord;
    for ( int i = 1; i <= 3; i++ ) {
        nodeCoord.push_back( & iTri.giveVertex(i) );
    }

    XfemManager *xMan = this->element->giveDomain()->giveXfemManager();
    EnrichmentItem *ei = xMan->giveEnrichmentItem(iEnrItemIndex);

    std :: vector< FloatArray >intersecPoints;
    std :: vector< int >intersecEdgeInd;

    std :: vector< double >minDistArcPos;
    ei->computeIntersectionPoints(intersecPoints, intersecEdgeInd, element, iTri, minDistArcPos);


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
        int nNodes = 3;
        std :: vector< FloatArray >edgeCoords, nodeCoords;

        FloatArray tipCoord;
        int dim = element->giveDofManager(1)->giveCoordinates()->giveSize();
        tipCoord.resize(dim);

        bool foundTip = false;
        double tipArcPos = -1.0;

        if ( ei->giveElementTipCoord(tipCoord, tipArcPos, element->giveNumber(), iTri) ) {
            foundTip = true;
        }

        if ( foundTip ) {
            for ( int i = 1; i <= nNodes; i++ ) {
                // Store edge points
                if ( i == intersecEdgeInd [ 0 ] ) {
                    // Take the intersection point ...
                    edgeCoords.push_back(intersecPoints [ 0 ]);
                } else {
                    // ... or the center of the edge.

                    FloatArray coordS, coordE;

                    // Global coordinates of vertices
                    switch ( i ) {
                    case 1:
                        coordS = * ( nodeCoord [ 0 ] );
                        coordE = * ( nodeCoord [ 1 ] );
                        break;
                    case 2:
                        coordS = * ( nodeCoord [ 1 ] );
                        coordE = * ( nodeCoord [ 2 ] );
                        break;

                    case 3:
                        coordS = * ( nodeCoord [ 2 ] );
                        coordE = * ( nodeCoord [ 0 ] );
                        break;
                    default:
                        break;
                    }


                    FloatArray coordEdge;
                    coordEdge = 0.5 * coordS + 0.5 * coordE;
                    edgeCoords.push_back(coordEdge);
                }

                // Store node coords
                const FloatArray &coord = iTri.giveVertex(i);
                nodeCoords.push_back(coord);
            }

            oPointPartitions.resize( ( 2 * nNodes ) );

            // Divide into subdomains
            for ( int i = 1; i <= nNodes; i++ ) {
                ////////////////
                // Take edge center or intersection point
                oPointPartitions [ 2 * i - 1 ].push_back(edgeCoords [ i - 1 ]);

                // Take crack tip position
                oPointPartitions [ 2 * i - 1 ].push_back(tipCoord);

                // Take node
                oPointPartitions [ 2 * i - 1 ].push_back( * ( element->giveDofManager(i)->giveCoordinates() ) );

                ////////////////
                // Take edge center or intersection point
                oPointPartitions [ 2 * i - 2 ].push_back(edgeCoords [ i - 1 ]);

                // Take next node
                if ( i == nNodes ) {
                    oPointPartitions [ 2 * i - 2 ].push_back( iTri.giveVertex(1) );
                } else {
                    oPointPartitions [ 2 * i - 2 ].push_back( iTri.giveVertex(i + 1) );
                }

                // Take crack tip position
                oPointPartitions [ 2 * i - 2 ].push_back(tipCoord);
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

void XfemElementInterface :: XfemElementInterface_computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    if( element->giveDomain()->hasXfemManager() ) {

		XfemManager *xMan = element->giveDomain()->giveXfemManager();
		int nEI = xMan->giveNumberOfEnrichmentItems();
		CrossSection *cs = NULL;

		for ( int i = 1; i <= nEI; i++ ) {
			EnrichmentItem &ei = * ( xMan->giveEnrichmentItem(i) );
			if ( ei.isMaterialModified(* gp, * element, cs) ) {
				StructuralCrossSection *structCS = dynamic_cast< StructuralCrossSection * >( cs );

				if ( structCS != NULL ) {
					structCS->giveCharMaterialStiffnessMatrix(answer, rMode, gp, tStep);
					return;
				} else {
					OOFEM_ERROR("XfemElementInterface :: XfemElementInterface_computeConstitutiveMatrixAt: failed to fetch StructuralMaterial\n");
				}
			}
		}
    }

    // If no enrichment modifies the material,
    // compute stiffness based on the bulk material.
    StructuralElement &structEl = dynamic_cast< StructuralElement & >( ( * element ) );
    structEl.StructuralElement :: computeConstitutiveMatrixAt(answer, rMode, gp, tStep);
}

void XfemElementInterface :: XfemElementInterface_computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    StructuralCrossSection *cs = dynamic_cast< StructuralCrossSection * >( element->giveCrossSection() );
    if ( cs == NULL ) {
        OOFEM_ERROR("XfemElementInterface :: XfemElementInterface_computeStressVector: cs == NULL.\n");
    }

    cs->giveRealStresses(answer, gp, strain, tStep);

    if( element->giveDomain()->hasXfemManager() ) {

		XfemManager *xMan = element->giveDomain()->giveXfemManager();

		int nEI = xMan->giveNumberOfEnrichmentItems();

		CrossSection *csInclusion = NULL;
		for ( int i = 1; i <= nEI; i++ ) {
			EnrichmentItem &ei = * ( xMan->giveEnrichmentItem(i) );
			if ( ei.isMaterialModified(* gp, * element, csInclusion) ) {
				StructuralCrossSection *structCSInclusion = dynamic_cast< StructuralCrossSection * >( csInclusion );

				if ( structCSInclusion != NULL ) {
					structCSInclusion->giveRealStresses(answer, gp, strain, tStep);
					return;
				} else {
					OOFEM_ERROR("PlaneStress2dXfem :: computeStressVector: failed to fetch StructuralCrossSection\n");
				}
			}
		}
    }
}

void XfemElementInterface :: computeCohesiveForces(FloatArray &answer, TimeStep *tStep)
{
    if ( hasCohesiveZone() ) {
        FloatArray solVec;
        element->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, solVec);

        size_t numSeg = mpCZIntegrationRules.size();
        for ( size_t segIndex = 0; segIndex < numSeg; segIndex++ ) {
            int numGP = mpCZIntegrationRules [ segIndex ]->giveNumberOfIntegrationPoints();

            for ( int gpIndex = 0; gpIndex < numGP; gpIndex++ ) {
                GaussPoint &gp = * ( mpCZIntegrationRules [ segIndex ]->getIntegrationPoint(gpIndex) );

                ////////////////////////////////////////////////////////
                // Compute a (slightly modified) N-matrix

                FloatMatrix NMatrix;
                computeNCohesive(NMatrix, gp, mCZEnrItemIndices [ segIndex ]);
                ////////////////////////////////////////////////////////


                // Traction
                FloatArray T2D;



                // Fetch material status and get normal
                StructuralInterfaceMaterialStatus *ms = dynamic_cast< StructuralInterfaceMaterialStatus * >( mpCZMat->giveStatus(& gp) );
                if ( ms == NULL ) {
                    OOFEM_ERROR("In XfemElementInterface :: computeCohesiveForces(): Failed to fetch material status.\n");
                }

                FloatArray crackNormal( ms->giveNormal() );

                // Compute jump vector
                FloatArray jump2D;
                computeDisplacementJump(gp, jump2D, solVec, NMatrix);

                computeGlobalCohesiveTractionVector(T2D, jump2D, crackNormal, NMatrix, gp, tStep);

                // Add to internal force
                FloatArray NTimesT;

                NTimesT.beTProductOf(NMatrix, T2D);
                CrossSection *cs  = element->giveCrossSection();
                double thickness = cs->give(CS_Thickness, & gp);
                double dA = thickness * gp.giveWeight();
                answer.add(dA, NTimesT);
            }
        }
    }
}

void XfemElementInterface :: computeGlobalCohesiveTractionVector(FloatArray &oT, const FloatArray &iJump, const FloatArray &iCrackNormal, const FloatMatrix &iNMatrix, GaussPoint &iGP, TimeStep *tStep)
{
    FloatMatrix F;
    F.resize(3, 3);
    F.beUnitMatrix();     // TODO: Compute properly


    FloatArray jump3D;
    jump3D.setValues(3, iJump.at(1), iJump.at(2), 0.0);


    FloatArray crackNormal3D;
    crackNormal3D.setValues(3, iCrackNormal.at(1), iCrackNormal.at(2), 0.0);

    FloatArray ez;
    ez.setValues(3, 0.0, 0.0, 1.0);
    FloatArray crackTangent3D;
    crackTangent3D.beVectorProductOf(crackNormal3D, ez);

    FloatMatrix locToGlob(3, 3);
    locToGlob.setColumn(crackTangent3D, 1);
    locToGlob.setColumn(crackNormal3D, 2);
    locToGlob.setColumn(ez, 3);

    FloatArray TLoc(3), jump3DLoc, TLocRenumbered(3);
    jump3DLoc.beTProductOf(locToGlob, jump3D);

    FloatArray jump3DLocRenumbered;
    jump3DLocRenumbered.setValues( 3, jump3DLoc.at(3), jump3DLoc.at(1), jump3DLoc.at(2) );

    mpCZMat->giveFirstPKTraction_3d(TLocRenumbered, & iGP, jump3DLocRenumbered, F, tStep);

    TLoc.setValues( 3, TLocRenumbered.at(2), TLocRenumbered.at(3), TLocRenumbered.at(1) );


    FloatArray T;
    T.beProductOf(locToGlob, TLoc);

    oT.setValues( 2, T.at(1), T.at(2) );
}

void XfemElementInterface :: computeCohesiveTangent(FloatMatrix &answer, TimeStep *tStep)
{
    if ( hasCohesiveZone() ) {
        FloatArray solVec;
        element->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, solVec);

        size_t numSeg = mpCZIntegrationRules.size();

        for ( size_t segIndex = 0; segIndex < numSeg; segIndex++ ) {
            int numGP = mpCZIntegrationRules [ segIndex ]->giveNumberOfIntegrationPoints();

            for ( int gpIndex = 0; gpIndex < numGP; gpIndex++ ) {
                GaussPoint &gp = * ( mpCZIntegrationRules [ segIndex ]->getIntegrationPoint(gpIndex) );

                ////////////////////////////////////////////////////////
                // Compute a (slightly modified) N-matrix

                FloatMatrix NMatrix;
                computeNCohesive(NMatrix, gp, mCZEnrItemIndices [ segIndex ]);

                ////////////////////////////////////////////////////////

                // Compute jump vector
                FloatArray jump2D;
                computeDisplacementJump(gp, jump2D, solVec, NMatrix);

                FloatArray jump3D;
                jump3D.setValues( 3, 0.0, jump2D.at(1), jump2D.at(2) );

                // Compute traction
                FloatMatrix F;
                F.resize(3, 3);
                F.beUnitMatrix();                     // TODO: Compute properly

                FloatMatrix K3DRenumbered, K3DGlob;


                FloatMatrix K2D;
                K2D.resize(2, 2);
                K2D.zero();

                if ( mpCZMat->hasAnalyticalTangentStiffness() ) {
                    ///////////////////////////////////////////////////
                    // Analytical tangent

                    FloatMatrix K3D;
                    mpCZMat->give3dStiffnessMatrix_dTdj(K3DRenumbered, TangentStiffness, & gp, tStep);

                    K3D.resize(3, 3);
                    K3D.zero();
                    K3D.at(1, 1) = K3DRenumbered.at(2, 2);
                    K3D.at(1, 2) = K3DRenumbered.at(2, 3);
                    K3D.at(1, 3) = K3DRenumbered.at(2, 1);

                    K3D.at(2, 1) = K3DRenumbered.at(3, 2);
                    K3D.at(2, 2) = K3DRenumbered.at(3, 3);
                    K3D.at(2, 3) = K3DRenumbered.at(3, 1);

                    K3D.at(3, 1) = K3DRenumbered.at(1, 2);
                    K3D.at(3, 2) = K3DRenumbered.at(1, 3);
                    K3D.at(3, 3) = K3DRenumbered.at(1, 1);


                    // Fetch material status and get normal
                    StructuralInterfaceMaterialStatus *ms = dynamic_cast< StructuralInterfaceMaterialStatus * >( mpCZMat->giveStatus(& gp) );
                    if ( ms == NULL ) {
                        OOFEM_ERROR("In XfemElementInterface :: computeCohesiveForces(): Failed to fetch material status.\n");
                    }

                    FloatArray crackNormal( ms->giveNormal() );

                    FloatArray crackNormal3D;
                    crackNormal3D.setValues(3, crackNormal.at(1), crackNormal.at(2), 0.0);

                    FloatArray ez;
                    ez.setValues(3, 0.0, 0.0, 1.0);
                    FloatArray crackTangent3D;
                    crackTangent3D.beVectorProductOf(crackNormal3D, ez);

                    FloatMatrix locToGlob(3, 3);
                    locToGlob.setColumn(crackTangent3D, 1);
                    locToGlob.setColumn(crackNormal3D, 2);
                    locToGlob.setColumn(ez, 3);


                    FloatMatrix tmp3(3, 3);
                    tmp3.beProductTOf(K3D, locToGlob);
                    K3DGlob.beProductOf(locToGlob, tmp3);

                    K2D.at(1, 1) = K3DGlob.at(1, 1);
                    K2D.at(1, 2) = K3DGlob.at(1, 2);
                    K2D.at(2, 1) = K3DGlob.at(2, 1);
                    K2D.at(2, 2) = K3DGlob.at(2, 2);
                } else {
                    ///////////////////////////////////////////////////
                    // Numerical tangent
                    double eps = 1.0e-9;

                    FloatArray T, TPert;

                    // Fetch material status and get normal
                    StructuralInterfaceMaterialStatus *ms = dynamic_cast< StructuralInterfaceMaterialStatus * >( mpCZMat->giveStatus(& gp) );
                    if ( ms == NULL ) {
                        OOFEM_ERROR("In XfemElementInterface :: computeCohesiveForces(): Failed to fetch material status.\n");
                    }

                    FloatArray crackNormal( ms->giveNormal() );

                    computeGlobalCohesiveTractionVector(T, jump2D, crackNormal, NMatrix, gp, tStep);


                    FloatArray jump2DPert;


                    jump2DPert = jump2D;
                    jump2DPert.at(1) += eps;
                    computeGlobalCohesiveTractionVector(TPert, jump2DPert, crackNormal, NMatrix, gp, tStep);

                    K2D.at(1, 1) = ( TPert.at(1) - T.at(1) ) / eps;
                    K2D.at(2, 1) = ( TPert.at(2) - T.at(2) ) / eps;

                    jump2DPert = jump2D;
                    jump2DPert.at(2) += eps;
                    computeGlobalCohesiveTractionVector(TPert, jump2DPert, crackNormal, NMatrix, gp, tStep);

                    K2D.at(1, 2) = ( TPert.at(1) - T.at(1) ) / eps;
                    K2D.at(2, 2) = ( TPert.at(2) - T.at(2) ) / eps;

                    computeGlobalCohesiveTractionVector(T, jump2D, crackNormal, NMatrix, gp, tStep);
                }


                FloatMatrix tmp, tmp2;
                tmp.beProductOf(K2D, NMatrix);
                tmp2.beTProductOf(NMatrix, tmp);

                CrossSection *cs  = element->giveCrossSection();
                double thickness = cs->give(CS_Thickness, & gp);
                double dA = thickness * gp.giveWeight();
                answer.add(dA, tmp2);
            }
        }
    }
}

void XfemElementInterface :: computeCohesiveTangentAt(FloatMatrix &answer, TimeStep *tStep)
{
    if ( hasCohesiveZone() ) {
        printf("Entering XfemElementInterface :: computeCohesiveTangentAt().\n");
    }
}

void XfemElementInterface :: XfemElementInterface_computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity)
{
    StructuralElement *structEl = dynamic_cast< StructuralElement * >(element);
    if ( structEl == NULL ) {
        OOFEM_ERROR("Error in XfemElementInterface :: XfemElementInterface_computeConsistentMassMatrix().\n");
    }

    int ndofs = structEl->computeNumberOfDofs();
    double density, dV;
    FloatMatrix n;
    IntegrationRule *iRule = element->giveIntegrationRule(0);
    IntArray mask;

    answer.resize(ndofs, ndofs);
    answer.zero();
    if ( !structEl->isActivated(tStep) ) {
        return;
    }

    structEl->giveMassMtrxIntegrationgMask(mask);

    mass = 0.;

    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        structEl->computeNmatrixAt(* ( gp->giveLocalCoordinates() ), n);
        density = structEl->giveMaterial()->give('d', gp);

        if ( ipDensity != NULL ) {
            // Override density if desired
            density = * ipDensity;
        }

        dV = structEl->computeVolumeAround(gp);
        mass += density * dV;

        if ( mask.isEmpty() ) {
            answer.plusProductSymmUpper(n, n, density * dV);
        } else {
            for ( int i = 1; i <= ndofs; i++ ) {
                for ( int j = i; j <= ndofs; j++ ) {
                    double summ = 0.;
                    for ( int k = 1; k <= n.giveNumberOfRows(); k++ ) {
                        if ( mask.at(k) == 0 ) {
                            continue;
                        }

                        summ += n.at(k, i) * n.at(k, j);
                    }

                    answer.at(i, j) += summ * density * dV;
                }
            }
        }
    }

    answer.symmetrized();
}

IRResultType
XfemElementInterface :: initializeCZFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    int material = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, material, _IFT_XfemElementInterface_CohesiveZoneMaterial);
    mCZMaterialNum = material;
    //    printf("In XfemElementInterface :: initializeCZFrom(): mCZMaterialNum: %d\n", mCZMaterialNum );


    // Number of Gauss points used when integrating the cohesive zone
    IR_GIVE_OPTIONAL_FIELD(ir, mCSNumGaussPoints, _IFT_XfemElementInterface_NumIntPointsCZ);
    //    printf("mCSNumGaussPoints: %d\n", mCSNumGaussPoints );

    int planeStrainFlag = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, planeStrainFlag, _IFT_XfemElementInterface_PlaneStrain);
    if ( planeStrainFlag == 1 ) {
        mUsePlaneStrain = true;
    }

    return IRRT_OK;
}

MaterialMode XfemElementInterface :: giveMaterialMode()
{
    if ( mUsePlaneStrain ) {
        return _PlaneStrain;
    } else {
        return _PlaneStress;
    }
}

void XfemElementInterface :: giveCZInputRecord(DynamicInputRecord &input)
{
    if ( mCZMaterialNum > 0 ) {
        input.setField(mCZMaterialNum, _IFT_XfemElementInterface_CohesiveZoneMaterial);
    }

    if ( mUsePlaneStrain ) {
        input.setField(1, _IFT_XfemElementInterface_PlaneStrain);
    }
}

void XfemElementInterface :: initializeCZMaterial()
{
    if ( mCZMaterialNum > 0 ) {
        mpCZMat = dynamic_cast< StructuralInterfaceMaterial * >( this->element->giveDomain()->giveMaterial(mCZMaterialNum) );

        if ( mpCZMat == NULL ) {
            OOFEM_ERROR("In XfemElementInterface :: initializeCZMaterial(): Failed to fetch pointer for mpCZMat.\n");
        }
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

void XfemElementInterface :: computeDisplacementJump(GaussPoint &iGP, FloatArray &oJump, const FloatArray &iSolVec, const FloatMatrix &iNMatrix)
{
    const int dim = 2;
    oJump.resize(dim);
    oJump.beProductOf(iNMatrix, iSolVec);
}

void XfemElementInterface :: computeNCohesive(FloatMatrix &oN, GaussPoint &iGP, int iEnrItemIndex)
{
    const int dim = 2;

    FloatArray Nc, globalCoord, localCoord;
    globalCoord = * ( iGP.giveCoordinates() );
    element->computeLocalCoordinates(localCoord, globalCoord);
    FEInterpolation *interp = element->giveInterpolation();
    interp->evalN( Nc, localCoord, FEIElementGeometryWrapper(element) );

    const int nDofMan = element->giveNumberOfDofManagers();

    // XFEM part of N-matrix
    XfemManager *xMan = element->giveDomain()->giveXfemManager();


    int counter = nDofMan * dim;

    std :: vector< std :: vector< double > > Nd(nDofMan);

    for ( int j = 1; j <= nDofMan; j++ ) {
        DofManager *dMan = element->giveDofManager(j);

        // Compute the total number of enrichments for node j
        int numEnrNode = 0;
        for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
            EnrichmentItem *ei = xMan->giveEnrichmentItem(i);
            if ( ei->isDofManEnriched(* dMan) ) {
                numEnrNode += ei->giveNumDofManEnrichments(* dMan);
            }
        }

        std :: vector< double > &NdNode = Nd [ j - 1 ];
        NdNode.assign(numEnrNode, 0.0);


        int globalNodeInd = dMan->giveGlobalNumber();


        int ndNodeInd = 0;
        for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
            EnrichmentItem *ei = xMan->giveEnrichmentItem(i);

            if ( ei->isDofManEnriched(* dMan) ) {
                int numEnr = ei->giveNumDofManEnrichments(* dMan);

                std :: vector< double >efJumps;
                ei->evaluateEnrFuncJumps(efJumps, globalNodeInd);

                for ( int k = 0; k < numEnr; k++ ) {
                    if ( i == iEnrItemIndex ) {
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
            NTot.at(column) = NdNode [ j - 1 ];
            column++;
        }
    }

    oN.beNMatrixOf(NTot, 2);
}
} // end namespace oofem
