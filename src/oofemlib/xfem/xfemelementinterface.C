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
//#include "patch.h"
#include "patchintegrationrule.h"
#include "delaunay.h"
#include "xfemmanager.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "enrichmentdomain.h"
#include "dynamicinputrecord.h"

#include "structuralinterfacematerial.h"

#include "XFEMDebugTools.h"
#include <string>
#include <sstream>
#include <math.h>

namespace oofem {

XfemElementInterface :: XfemElementInterface(Element *e) :
Interface(),
element(e),
mpCZMat(NULL),
mCZMaterialNum(-1),
mCSNumGaussPoints(4),
mpCZIntegrationRule(NULL),
mCrackLength(0.0),
mCZEnrItemIndex(-1)
{

}

XfemElementInterface :: ~XfemElementInterface()
{
	if( mpCZIntegrationRule != NULL ) {
		delete mpCZIntegrationRule;
		mpCZIntegrationRule = NULL;
	}

}

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


    std :: vector< FloatMatrix > Bd(nDofMan); // One Bd per node

    int counter = nDofMan * dim;

    std::vector< std::vector<double> > Nd(nDofMan);

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

        std::vector<double> &NdNode = Nd [ j - 1 ];
        NdNode.assign(numEnrNode, 0.0);


        int globalNodeInd = dMan->giveGlobalNumber();


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


				for(int k = 0; k < numEnr; k++) {
					NdNode[k] = ( efGP[k] - efNode[k] ) * Nc.at(j) ;
					counter ++;
				}
			}

		}
    }

    int numN = nDofMan;

    for ( int j = 1; j <= nDofMan; j++ ) {
    	numN += Nd[j-1].size();
    }

    FloatArray NTot;
    NTot.resize(numN);
    NTot.zero();
    int column = 1;

    for(int i = 1; i <= nDofMan; i++) {
        NTot.at(column) = Nc.at(i);
        column++;

        const std::vector<double> &NdNode = Nd[i-1];
        for(size_t j = 1; j <= NdNode.size(); j++) {
        	NTot.at(column) = NdNode[j-1];
        	column++;
        }
    }

    oAnswer.beNMatrixOf(NTot,2);
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
        FloatArray startPoint, endPoint;
        int enrItemInd = -1;
        this->XfemElementInterface_prepareNodesForDelaunay(pointPartitions, startPoint, endPoint, enrItemInd);
        mCrackLength = startPoint.distance(endPoint);
        mCZEnrItemIndex = enrItemInd;

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

        if(mpCZMat == NULL && mCZMaterialNum > 0) {
        	initializeCZMaterial();
        }

        if( mpCZMat != NULL ) {

        	if( mpCZIntegrationRule != NULL ) {
        		delete mpCZIntegrationRule;
        	}

            int czRuleNum = 1;
        	mpCZIntegrationRule = new GaussIntegrationRule(czRuleNum, element);
        	const FloatArray **coords = new const FloatArray*[2];
        	coords[0] = new FloatArray(startPoint);
        	coords[1] = new FloatArray(endPoint);
        	mpCZIntegrationRule->SetUpPointsOn2DEmbeddedLine(mCSNumGaussPoints, matMode, coords);

        	delete coords[0];
        	delete coords[1];
        	delete [] coords;


#if XFEM_DEBUG_VTK > 0
        	////////////////////////////////////////////////////////////////////////
        	// Write CZ GP to VTK

        	std :: vector< FloatArray >czGPCoord;

        	for(int i = 0; i < mpCZIntegrationRule->giveNumberOfIntegrationPoints(); i++) {
        		czGPCoord.push_back( *(mpCZIntegrationRule->getIntegrationPoint(i)->giveCoordinates()) );
        	}

            double time = 0.0;

            Element *el = element;

            Domain *dom = el->giveDomain();
            if(dom != NULL) {
            	EngngModel *em = dom->giveEngngModel();
            	if(em != NULL) {
            		TimeStep *ts = em->giveCurrentStep();
            		if(ts != NULL) {
            			time = ts->giveTargetTime();
            		}
            	}
            }

            int elIndex = el->giveGlobalNumber();
            std :: stringstream str;
            str << "CZGaussPointsTime" << time << "El" << elIndex << ".vtk";
            std :: string name = str.str();

            XFEMDebugTools :: WritePointsToVTK(name, czGPCoord);
        	////////////////////////////////////////////////////////////////////////
#endif

        }
    }
}

void XfemElementInterface :: XfemElementInterface_prepareNodesForDelaunay(std :: vector< std :: vector< FloatArray > > &oPointPartitions, FloatArray &oCrackStartPoint, FloatArray &oCrackEndPoint, int &oEnrItemIndex)
{
    XfemManager *xMan = this->element->giveDomain()->giveXfemManager();

    std :: vector< FloatArray >intersecPoints;
    std :: vector< int >intersecEdgeInd;

    // TODO:    Can we do this recursively to achieve proper splitting
    //			when several enrichment items interact with the
    //			same element?
    for ( int eiInd = 1; eiInd <= xMan->giveNumberOfEnrichmentItems(); eiInd++ ) {
        xMan->giveEnrichmentItem(eiInd)->computeIntersectionPoints(intersecPoints, intersecEdgeInd, element);



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


			// Export start and end points of
			// the intersection line.
			oCrackStartPoint = intersecPoints[0];
			oCrackEndPoint = intersecPoints[1];

			oEnrItemIndex = eiInd;

			// TODO: Handle multiple intersections
			return;
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

				// Export start and end points of
				// the intersection line.
				oCrackStartPoint = intersecPoints[0];
				oCrackEndPoint = tipCoord;

			} // If a tip was found
			else
			{
				oPointPartitions.resize(1);

				for ( int i = 1; i <= this->element->giveNumberOfDofManagers(); i++ )
				{
					const FloatArray &nodeCoord = *element->giveDofManager(i)->giveCoordinates();
					oPointPartitions[0].push_back(nodeCoord);
				}

				// Export start and end points of
				// the intersection line.
				oCrackStartPoint = intersecPoints[0];
				oCrackEndPoint = intersecPoints[0];
			}

			oEnrItemIndex = eiInd;

			return;
		}

    } // for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ )
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

void XfemElementInterface :: computeCohesiveForces(FloatArray &answer, TimeStep *tStep)
{

	if( hasCohesiveZone() ) {

		FloatArray solVec;
		element->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, solVec);


		int numGP = mpCZIntegrationRule->giveNumberOfIntegrationPoints();

		for(int gpIndex = 0; gpIndex < numGP; gpIndex++) {

			GaussPoint &gp = *(mpCZIntegrationRule->getIntegrationPoint(gpIndex));

		    ////////////////////////////////////////////////////////
		    // Compute a (slightly modified) N-matrix

			FloatMatrix NMatrix;
			computeNCohesive(NMatrix, gp);

		    ////////////////////////////////////////////////////////


			// Traction
			FloatArray T2D;


			FloatArray crackNormal;
			if( computeNormalInPoint( *(gp.giveCoordinates()), crackNormal ) ) {


				// Compute jump vector
				FloatArray jump2D;
				computeDisplacementJump(gp, jump2D, solVec, NMatrix);

				computeGlobalCohesiveTractionVector(T2D, jump2D, crackNormal, NMatrix, gp, tStep);

				// Add to internal force
				FloatArray NTimesT;

				NTimesT.beTProductOf(NMatrix, T2D);
				CrossSection *cs  = element->giveCrossSection();
				double thickness = cs->give(CS_Thickness);
				double dA = 0.5*mCrackLength*thickness*gp.giveWeight();
				answer.add(dA, NTimesT);
			}
			else {
				OOFEM_ERROR("In XfemElementInterface :: computeCohesiveForces: Failed to compute normal in Gauss point.\n");
			}
		}




	}
}

void XfemElementInterface :: computeGlobalCohesiveTractionVector(FloatArray &oT, const FloatArray &iJump, const FloatArray &iCrackNormal, const FloatMatrix &iNMatrix, GaussPoint &iGP, TimeStep *tStep)
{
	FloatMatrix F;
	F.resize(3,3);
	F.beUnitMatrix(); // TODO: Compute properly


	FloatArray jump3D;
	jump3D.setValues(3, iJump.at(1), iJump.at(2), 0.0);


	FloatArray crackNormal3D;
	crackNormal3D.setValues(3, iCrackNormal.at(1), iCrackNormal.at(2), 0.0);

	FloatArray ez;
	ez.setValues(3, 0.0, 0.0, 1.0);
	FloatArray crackTangent3D;
	crackTangent3D.beVectorProductOf(crackNormal3D, ez);

	FloatMatrix locToGlob(3,3);
	locToGlob.setColumn(crackTangent3D, 1);
	locToGlob.setColumn(crackNormal3D, 2);
	locToGlob.setColumn(ez, 3);

	FloatArray TLoc(3), jump3DLoc, TLocRenumbered(3);
	jump3DLoc.beTProductOf(locToGlob, jump3D);

	FloatArray jump3DLocRenumbered;
	jump3DLocRenumbered.setValues(3, jump3DLoc.at(3), jump3DLoc.at(1), jump3DLoc.at(2));

	mpCZMat->giveFirstPKTraction_3d(TLocRenumbered, &iGP, jump3DLocRenumbered, F, tStep);

	TLoc.setValues(3, TLocRenumbered.at(2), TLocRenumbered.at(3), TLocRenumbered.at(1));


	FloatArray T;
	T.beProductOf(locToGlob, TLoc);

	oT.setValues(2, T.at(1), T.at(2) );
}

void XfemElementInterface :: computeCohesiveTangent(FloatMatrix &answer, TimeStep *tStep)
{

	if( hasCohesiveZone() ) {

		FloatArray solVec;
		element->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, solVec);


		int numGP = mpCZIntegrationRule->giveNumberOfIntegrationPoints();

		for(int gpIndex = 0; gpIndex < numGP; gpIndex++) {

			GaussPoint &gp = *(mpCZIntegrationRule->getIntegrationPoint(gpIndex));

			////////////////////////////////////////////////////////
			// Compute a (slightly modified) N-matrix

			FloatMatrix NMatrix;
			computeNCohesive(NMatrix, gp);

			////////////////////////////////////////////////////////

			// Compute jump vector
			FloatArray jump2D;
			computeDisplacementJump(gp, jump2D, solVec, NMatrix);

			FloatArray jump3D;
			jump3D.setValues(3, 0.0, jump2D.at(1), jump2D.at(2));

			// Compute traction
			FloatMatrix F;
			F.resize(3,3);
			F.beUnitMatrix(); // TODO: Compute properly

			FloatMatrix K3DRenumbered, K3DGlob;


			FloatMatrix K2D;
			K2D.resize(2,2);
			K2D.zero();

			if( mpCZMat->hasAnalyticalTangentStiffness() ){
				///////////////////////////////////////////////////
				// Analytical tangent

				FloatMatrix K3D;
				mpCZMat->give3dStiffnessMatrix_dTdj(K3DRenumbered, TangentStiffness, &gp, tStep);

				K3D.resize(3,3);
				K3D.zero();
				K3D.at(1,1) = K3DRenumbered.at(2,2);
				K3D.at(1,2) = K3DRenumbered.at(2,3);
				K3D.at(1,3) = K3DRenumbered.at(2,1);

				K3D.at(2,1) = K3DRenumbered.at(3,2);
				K3D.at(2,2) = K3DRenumbered.at(3,3);
				K3D.at(2,3) = K3DRenumbered.at(3,1);

				K3D.at(3,1) = K3DRenumbered.at(1,2);
				K3D.at(3,2) = K3DRenumbered.at(1,3);
				K3D.at(3,3) = K3DRenumbered.at(1,1);


				FloatArray crackNormal;
				computeNormalInPoint( *(gp.giveCoordinates()), crackNormal );
				FloatArray crackNormal3D;
				crackNormal3D.setValues(3, crackNormal.at(1), crackNormal.at(2), 0.0);

				FloatArray ez;
				ez.setValues(3, 0.0, 0.0, 1.0);
				FloatArray crackTangent3D;
				crackTangent3D.beVectorProductOf(crackNormal3D, ez);

				FloatMatrix locToGlob(3,3);
				locToGlob.setColumn(crackTangent3D, 1);
				locToGlob.setColumn(crackNormal3D, 2);
				locToGlob.setColumn(ez, 3);


				FloatMatrix tmp3(3,3);
				tmp3.beProductTOf(K3D, locToGlob);
				K3DGlob.beProductOf(locToGlob, tmp3);

				K2D.at(1,1) = K3DGlob.at(1,1);
				K2D.at(1,2) = K3DGlob.at(1,2);
				K2D.at(2,1) = K3DGlob.at(2,1);
				K2D.at(2,2) = K3DGlob.at(2,2);
			}
			else {

				///////////////////////////////////////////////////
				// Numerical tangent
				double eps = 1.0e-9;

				FloatArray T, TPert;

				FloatArray crackNormal;
				if( computeNormalInPoint( *(gp.giveCoordinates()), crackNormal ) ) {

					computeGlobalCohesiveTractionVector(T, jump2D, crackNormal, NMatrix, gp, tStep);


					FloatArray jump2DPert;


					jump2DPert = jump2D;
					jump2DPert.at(1) += eps;
					computeGlobalCohesiveTractionVector(TPert, jump2DPert, crackNormal, NMatrix, gp, tStep);

					K2D.at(1,1) = (TPert.at(1) - T.at(1))/eps;
					K2D.at(2,1) = (TPert.at(2) - T.at(2))/eps;

					jump2DPert = jump2D;
					jump2DPert.at(2) += eps;
					computeGlobalCohesiveTractionVector(TPert, jump2DPert, crackNormal, NMatrix, gp, tStep);


					K2D.at(1,2) = (TPert.at(1) - T.at(1))/eps;
					K2D.at(2,2) = (TPert.at(2) - T.at(2))/eps;

					computeGlobalCohesiveTractionVector(T, jump2D, crackNormal, NMatrix, gp, tStep);

				}
			}


			FloatMatrix tmp, tmp2;
			tmp.beProductOf(K2D, NMatrix);
			tmp2.beTProductOf(NMatrix, tmp);

			CrossSection *cs  = element->giveCrossSection();
			double thickness = cs->give(CS_Thickness);
			double dA = 0.5*mCrackLength*thickness*gp.giveWeight();
			answer.add(dA, tmp2);
		}

	}
}

void XfemElementInterface :: computeCohesiveTangentAt(FloatMatrix &answer, TimeStep *tStep)
{
	if( hasCohesiveZone() ) {

		printf("Entering XfemElementInterface :: computeCohesiveTangentAt().\n");

	}
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

    return IRRT_OK;
}

void XfemElementInterface :: giveCZInputRecord(DynamicInputRecord &input)
{
	if( mCZMaterialNum > 0 ) {
		input.setField(mCZMaterialNum, _IFT_XfemElementInterface_CohesiveZoneMaterial);
	}
}

void XfemElementInterface :: initializeCZMaterial()
{
    if ( mCZMaterialNum > 0 ) {
    	mpCZMat = dynamic_cast<StructuralInterfaceMaterial*>(this->element->giveDomain()->giveMaterial(mCZMaterialNum) );

    	if( mpCZMat == NULL ) {
    		OOFEM_ERROR("In XfemElementInterface :: initializeCZMaterial(): Failed to fetch pointer for mpCZMat.\n");
    	}
    }
}

void XfemElementInterface :: updateYourselfCZ(TimeStep *tStep)
{
	if(mpCZIntegrationRule != NULL) {
		mpCZIntegrationRule->updateYourself(tStep);

		if(mpCZMat != NULL){
			int numGP = mpCZIntegrationRule->giveNumberOfIntegrationPoints();
			for(int i = 0; i < numGP; i++) {
				GaussPoint *gp = mpCZIntegrationRule->getIntegrationPoint(i);
				mpCZMat->updateYourself(gp, tStep);
			}
		}
	}


}

void XfemElementInterface :: computeDisplacementJump(GaussPoint &iGP, FloatArray &oJump, const FloatArray &iSolVec, const FloatMatrix &iNMatrix)
{
	const int dim = 2;
	oJump.resize(dim);
    oJump.beProductOf(iNMatrix, iSolVec);
}

void XfemElementInterface :: computeNCohesive(FloatMatrix &oN, GaussPoint &iGP)
{

	const int dim = 2;

    FloatArray Nc, globalCoord, localCoord;
    globalCoord = *(iGP.giveCoordinates());
    element->computeLocalCoordinates(localCoord, globalCoord);
    FEInterpolation *interp = element->giveInterpolation();
    interp->evalN( Nc, localCoord, FEIElementGeometryWrapper(element) );

    const int nDofMan = element->giveNumberOfDofManagers();

    // XFEM part of N-matrix
    XfemManager *xMan = element->giveDomain()->giveXfemManager();


    int counter = nDofMan * dim;

    std::vector< std::vector<double> > Nd(nDofMan);

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

        std::vector<double> &NdNode = Nd [ j - 1 ];
        NdNode.assign(numEnrNode, 0.0);


        int globalNodeInd = dMan->giveGlobalNumber();


        int ndNodeInd = 0;
        for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
        	EnrichmentItem *ei = xMan->giveEnrichmentItem(i);

        	if ( ei->isDofManEnriched(* dMan) ) {

                int numEnr = ei->giveNumDofManEnrichments(* dMan);

                std::vector<double> efJumps;
                ei->evaluateEnrFuncJumps(efJumps, globalNodeInd);

				for(int k = 0; k < numEnr; k++) {

					if( i == mCZEnrItemIndex ) {
						NdNode[ndNodeInd] = efJumps[k] * Nc.at(j);
					}
					else {
						NdNode[ndNodeInd] = 0.0;
					}
					counter ++;
					ndNodeInd++;
				}
			}

		}
    }

    int numN = nDofMan;

    for ( int j = 1; j <= nDofMan; j++ ) {
    	numN += Nd[j-1].size();
    }

    FloatArray NTot;
    NTot.resize(numN);
    NTot.zero();
    int column = 1;

    for(int i = 1; i <= nDofMan; i++) {

//        NTot.at(column) = Nc.at(i); // We do not want the continuous part.
        column++;

        const std::vector<double> &NdNode = Nd[i-1];
        for(size_t j = 1; j <= NdNode.size(); j++) {
        	NTot.at(column) = NdNode[j-1];
        	column++;
        }
    }

    oN.beNMatrixOf(NTot,2);
}

bool XfemElementInterface :: computeNormalInPoint(const FloatArray &iGlobalCoord, FloatArray &oNormal)
{
	bool foundNormal = false;

    XfemManager *xMan = element->giveDomain()->giveXfemManager();
    const IntArray &elNodes = element->giveDofManArray();

    FloatMatrix dNdx;
    FloatArray N;
    FEInterpolation *interp = element->giveInterpolation();
    FloatArray localCoord;
    interp->global2local(localCoord, iGlobalCoord, FEIElementGeometryWrapper(element) );
    interp->evaldNdx( dNdx, localCoord, FEIElementGeometryWrapper(element) );
    interp->evalN( N, localCoord, FEIElementGeometryWrapper(element) );

    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
    	EnrichmentItem *ei = xMan->giveEnrichmentItem(i);

    	double levelSet = 0.0;
    	ei->interpLevelSet(levelSet, N, elNodes);

    	double levelSetTang = 0.0;
    	ei->interpLevelSetTangential(levelSetTang, N, elNodes);

    	FloatArray gradLevelSet;
    	ei->interpGradLevelSet(gradLevelSet, dNdx, elNodes);


    	// If the normal level set is sufficiently small,
    	// and the tangential level set is positive,
    	// we have found the correct enrichment item.
    	double normalTol = 1.0e-2; // TODO: Should be related to the element size.
    	if( fabs(levelSet) < normalTol &&  levelSetTang > 0.0) {

    		// If so, take the normal as the gradient of the level set function

    		oNormal = gradLevelSet;
    		if( oNormal.computeNorm() > 1.0e-12 ) {
    			oNormal.normalize();
    			return true;
    		}
    	}
    }

	return foundNormal;
}

} // end namespace oofem
