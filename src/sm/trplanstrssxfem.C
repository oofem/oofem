/*
 * trplanstrssxfem.C
 *
 *  Created on: Jun 3, 2013
 *      Author: svennine
 */

#include "trplanstrssxfem.h"

#include "node.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "rcm2.h"
#endif


#include "structuralmaterial.h"
#include "xfemelementinterface.h"
#include "enrichmentfunction.h"
#include "enrichmentitem.h"
#include "enrichmentdomain.h"
#include "structuralcrosssection.h"
#include "vtkxmlexportmodule.h"
#//ifdef __OOFEG
#include "patchintegrationrule.h"
//#endif
#include "delaunay.h"

#include "XFEMDebugTools.h"
#include <string>
#include <sstream>


namespace oofem {

REGISTER_Element( TrPlaneStress2dXFEM );


TrPlaneStress2dXFEM::~TrPlaneStress2dXFEM() {

}

int TrPlaneStress2dXFEM::checkConsistency()
{
	TrPlaneStress2d :: checkConsistency();
    return 1;
}
#if 0
void TrPlaneStress2dXFEM :: XfemElementInterface_partitionElement(AList< Triangle > *answer, std :: vector< FloatArray > &together)
{
    Delaunay dl;
    dl.triangulate(together, answer);
    return;
/*
	// Two cases may occur when partitioning: we will get a subdomain
	// with either 3 or 4 nodes.
	if( together->giveSize() == 3 )
	{
		// 3 nodes

        FloatArray *p1 = new FloatArray();
        * p1 = * ( together->at(1) );
        FloatArray *p2 = new FloatArray();
        * p2 = * ( together->at(2) );
        FloatArray *p3 = new FloatArray();
        * p3 = * ( together->at(3) );

        Triangle *triangle = new Triangle(p1, p2, p3);
        if ( !triangle->isOrientedAnticlockwise() ) {
            triangle->changeToAnticlockwise();
        }

        answer->put(1, triangle);

	}
	else
	{
		if( together->giveSize() == 4 )
		{
			// 4 nodes

			// The array *together contains the four vertices of the area to be subdivided.
			// The intersection points are located in the first two entries in the array and
			// the vertex points are located in the last two entries in the array. We do not know
			// if the last two entries are numbered such that the four vertices form a proper quad.
			// We have to check if this is the case and if not, we swap the last two points.
			int nodeMap[4] = {1, 2, 3, 4};



			//////////////////////////////////////////////////////////////////////
			// We have two options for the subdivision.
			//
			// Which option we choose will probably not matter in the end.
			// However, trying to get triangles of good quality can't do
			// any harm. Start by checking the mesh quality we
			// will get for the two choices.


	        FloatArray *p1 = new FloatArray();
	        * p1 = * ( together->at(nodeMap[0]) );
	        FloatArray *p2 = new FloatArray();
	        * p2 = * ( together->at(nodeMap[1]) );
	        FloatArray *p3 = new FloatArray();
	        * p3 = * ( together->at(nodeMap[2]) );
	        FloatArray *p4 = new FloatArray();
	        * p4 = * ( together->at(nodeMap[3]) );


//	        bPoint2 b1_tmp( p1->at(1), p1->at(2) );
//	        bPoint2 b2_tmp( p2->at(1), p2->at(2) );
//	        bPoint2 b3_tmp( p3->at(1), p3->at(2) );
//	        bPoint2 b4_tmp( p4->at(1), p4->at(2) );

//	        bPoint2 cutLine( b2_tmp.x() - b1_tmp.x(), b2_tmp.y() - b1_tmp.y() );
	        FloatArray cutLine;
	        cutLine.beDifferenceOf(*p2, *p1);

//	        bPoint2  elLine( b4_tmp.x() - b3_tmp.x(), b4_tmp.y() - b3_tmp.y() );
	        FloatArray elLine;
	        elLine.beDifferenceOf(*p2, *p3);

//	        if( bDot(cutLine, elLine) > 0 )
	        if( cutLine.dotProduct(elLine) > 0.0 )
	        {
//	        	printf("Permuting node map.\n");
	        	nodeMap[2] = 4;
	        	nodeMap[3] = 3;
	        }

	        * p1 = * ( together->at(nodeMap[0]) );
	        * p2 = * ( together->at(nodeMap[1]) );
	        * p3 = * ( together->at(nodeMap[2]) );
	        * p4 = * ( together->at(nodeMap[3]) );

//	        bPoint2 b1( p1->at(1), p1->at(2) );
//	        bPoint2 b2( p2->at(1), p2->at(2) );
//	        bPoint2 b3( p3->at(1), p3->at(2) );
//	        bPoint2 b4( p4->at(1), p4->at(2) );


	        // Subdivision alternative 1
	        // 1 - 2 - 4
	        // and
	        // 2 - 3 - 4

	        // Triangle 1
//	        double a1_1 = bDist( b1, b2 );
	        double a1_1 = p1->distance(*p2);
//	        double b1_1 = bDist( b2, b4 );
	        double b1_1 = p2->distance(*p4);
//	        double c1_1 = bDist( b4, b1 );
	        double c1_1 = p4->distance(*p1);

//	        bSeg2 AB1_1(b1, b2);
//	        double h1_1 = bDist(b4, AB1_1);
	        double h1_1 = p4->distance(*p1,*p2);

	        double A1_1 = 0.5*a1_1*h1_1;
	        double R1_1 = 0.25*a1_1*b1_1*c1_1/A1_1;

	        double sin_a1_1 = 0.5*a1_1/R1_1;
	        double sin_b1_1 = 0.5*b1_1/R1_1;
	        double sin_c1_1 = 0.5*c1_1/R1_1;

	        double rho1_1 = ( sin_a1_1 + sin_b1_1 +  sin_c1_1 )/( 2.0*sin_a1_1*sin_b1_1*sin_c1_1 );
//	        printf("rho1_1: %e ", rho1_1);

	        // Triangle 2
//	        double a2_1 = bDist( b2, b3 );
	        double a2_1 = p2->distance(*p3);
//	        double b2_1 = bDist( b3, b4 );
	        double b2_1 = p3->distance(*p4);
//	        double c2_1 = bDist( b4, b2 );
	        double c2_1 = p4->distance(*p2);

//	        bSeg2 AB2_1(b2, b3);
//	        double h2_1 = bDist(b4, AB2_1);
	        double h2_1 = p4->distance(*p2,*p3);

	        double A2_1 = 0.5*a2_1*h2_1;
	        double R2_1 = 0.25*a2_1*b2_1*c2_1/A2_1;

	        double sin_a2_1 = 0.5*a2_1/R2_1;
	        double sin_b2_1 = 0.5*b2_1/R2_1;
	        double sin_c2_1 = 0.5*c2_1/R2_1;

	        double rho2_1 = ( sin_a2_1 + sin_b2_1 +  sin_c2_1 )/( 2.0*sin_a2_1*sin_b2_1*sin_c2_1 );
//	        printf("rho2_1: %e\n", rho2_1);

	        double worst1 = max( fabs(rho1_1-2.0), fabs(rho2_1-2.0) );
//	        printf("worst1: %e\n", worst1 );

	        // Subdivision alternative 2
	        // 1 - 2 - 3
	        // and
	        // 1 - 3 - 4

	        // Triangle 1
//	        double a1_2 = bDist( b1, b2 );
	        double a1_2 = p1->distance(p2);
//	        double b1_2 = bDist( b2, b3 );
	        double b1_2 = p2->distance(*p3);
//	        double c1_2 = bDist( b3, b1 );
	        double c1_2 = p3->distance(*p1);

//	        bSeg2 AB1_2(b1, b2);
//	        double h1_2 = bDist(b3, AB1_2);
	        double h1_2 = p3->distance(*p1,*p2);

	        double A1_2 = 0.5*a1_2*h1_2;
	        double R1_2 = 0.25*a1_2*b1_2*c1_2/A1_2;

	        double sin_a1_2 = 0.5*a1_2/R1_2;
	        double sin_b1_2 = 0.5*b1_2/R1_2;
	        double sin_c1_2 = 0.5*c1_2/R1_2;

	        double rho1_2 = ( sin_a1_2 + sin_b1_2 +  sin_c1_2 )/( 2.0*sin_a1_2*sin_b1_2*sin_c1_2 );
//	        printf("rho1_2: %e ", rho1_2);

	        // Triangle 2
//	        double a2_2 = bDist( b1, b3 );
	        double a2_2 = p1->distance(*p3);
//	        double b2_2 = bDist( b3, b4 );
	        double b2_2 = p3->distance(*p4);
//	        double c2_2 = bDist( b4, b1 );
	        double c2_2 = p4->distance(*p1);

//	        bSeg2 AB2_2(b1, b3);
//	        double h2_2 = bDist(b4, AB2_2);
	        double h2_2 = p4->distance(*p1,*p3);

	        double A2_2 = 0.5*a2_2*h2_2;
	        double R2_2 = 0.25*a2_2*b2_2*c2_2/A2_2;

	        double sin_a2_2 = 0.5*a2_2/R2_2;
	        double sin_b2_2 = 0.5*b2_2/R2_2;
	        double sin_c2_2 = 0.5*c2_2/R2_2;

	        double rho2_2 = ( sin_a2_2 + sin_b2_2 +  sin_c2_2 )/( 2.0*sin_a2_2*sin_b2_2*sin_c2_2 );
//	        printf("rho2_2: %e\n", rho2_2);

	        double worst2 = max( fabs(rho1_2-2.0), fabs(rho2_2-2.0) );
//	        printf("worst2: %e\n", worst2 );


			//////////////////////////////////////////////////////////////////////

	        if( worst1 < worst2 )
	        {
	        	// Choose subdivision 1
		        // 1 - 2 - 4
		        // and
		        // 2 - 3 - 4

//	        	printf("Choosing subdivision 1.\n");

		        FloatArray *p1A = new FloatArray();
		        * p1A = * ( together->at(nodeMap[0]) );
		        FloatArray *p2A = new FloatArray();
		        * p2A = * ( together->at(nodeMap[1]) );
		        FloatArray *p3A = new FloatArray();
		        * p3A = * ( together->at(nodeMap[3]) );

		        Triangle *triangle1 = new Triangle(p1A, p2A, p3A);

		        if ( !triangle1->isOrientedAnticlockwise() ) {
		            triangle1->changeToAnticlockwise();
		        }
		        answer->put(1, triangle1);

		        FloatArray *p1B = new FloatArray();
		        * p1B = * ( together->at(nodeMap[1]) );
		        FloatArray *p2B = new FloatArray();
		        * p2B = * ( together->at(nodeMap[2]) );
		        FloatArray *p3B = new FloatArray();
		        * p3B = * ( together->at(nodeMap[3]) );

		        Triangle *triangle2 = new Triangle(p1B, p2B, p3B);

		        if ( !triangle2->isOrientedAnticlockwise() ) {
		            triangle2->changeToAnticlockwise();
		        }
		        answer->put(2, triangle2);

	        }
	        else
	        {
	        	// Choose subdivision 2
		        // 1 - 2 - 3
		        // and
		        // 1 - 3 - 4

//	        	printf("Choosing subdivision 2.\n");

		        FloatArray *p1A = new FloatArray();
		        * p1A = * ( together->at(nodeMap[0]) );
		        FloatArray *p2A = new FloatArray();
		        * p2A = * ( together->at(nodeMap[1]) );
		        FloatArray *p3A = new FloatArray();
		        * p3A = * ( together->at(nodeMap[2]) );

		        Triangle *triangle1 = new Triangle(p1A, p2A, p3A);

		        if ( !triangle1->isOrientedAnticlockwise() ) {
		            triangle1->changeToAnticlockwise();
		        }
		        answer->put(1, triangle1);

		        FloatArray *p1B = new FloatArray();
		        * p1B = * ( together->at(nodeMap[0]) );
		        FloatArray *p2B = new FloatArray();
		        * p2B = * ( together->at(nodeMap[2]) );
		        FloatArray *p3B = new FloatArray();
		        * p3B = * ( together->at(nodeMap[3]) );

		        Triangle *triangle2 = new Triangle(p1B, p2B, p3B);

		        if ( !triangle2->isOrientedAnticlockwise() ) {
		            triangle2->changeToAnticlockwise();
		        }
		        answer->put(2, triangle2);

	        }

		}
	}
*/
}

void TrPlaneStress2dXFEM :: XfemElementInterface_updateIntegrationRule()
{

	XfemElementInterface :: XfemElementInterface_updateIntegrationRule();
	return;
/*
    XfemManager *xMan = this->element->giveDomain()->giveXfemManager();


//        IntArray activeEI;
//        xMan->giveActiveEIsFor(activeEI, element);


        std::vector< FloatArray > intersecPoints;
        std::vector< int > intersecEdgeInd;

//        printf("In TrPlaneStress2dXFEM :: XfemElementInterface_updateIntegrationRule(): Number of enrichment items: %d\n", xMan->giveNumberOfEnrichmentItems() );
        for( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++)
        {
        	EnrichmentItem *ei = xMan->giveEnrichmentItem(i);
        	ei->computeIntersectionPoints( intersecPoints, intersecEdgeInd, element);
        }


        // TODO: Right now, we only consider the case when we have exactly two
        //		 intersection points, i.e. one crack that passes through the
        //		 whole element.
        if( intersecPoints.size() == 2 )
        {
            AList< FloatArray > subDomain1, subDomain2;

            for ( int i = 1; i <= int(intersecPoints.size()); i++ ) {
                int sz = subDomain1.giveSize();
                subDomain1.put( sz + 1, new FloatArray(intersecPoints[i-1]) );

                FloatArray ip = intersecPoints[i-1];
                FloatArray *ipCopy = new FloatArray(ip);
                int sz2 = subDomain2.giveSize();
                subDomain2.put(sz2 + 1, ipCopy);
            }

            double x1 = intersecPoints[0].at(1);
            double x2 = intersecPoints[1].at(1);
            double y1 = intersecPoints[0].at(2);
            double y2 = intersecPoints[1].at(2);

            double tx = x2 - x1;
            double ty = y2 - y1;

            double nx = -ty;
            double ny =  tx;

            for ( int i = 1; i <= this->element->giveNumberOfDofManagers(); i++ )
            {
                double x = element->giveDofManager(i)->giveCoordinates()->at(1);
                double y = element->giveDofManager(i)->giveCoordinates()->at(2);

                double px = x1 - x;
                double py = y1 - y;

                FloatArray *node = element->giveDofManager(i)->giveCoordinates();
                FloatArray *nodesCopy = new FloatArray(*node);

                if( (px*nx + py*ny) > 0 )
                {
                    int sz = subDomain1.giveSize();
                    subDomain1.put(sz + 1, nodesCopy);
                }
                else
                {
                    int sz = subDomain2.giveSize();
                    subDomain2.put(sz + 1, nodesCopy);
                }

            }



            // Now we have two subdomains, one with three nodes and one with four nodes.

            AList< Triangle >triangles;
            AList< Triangle >triangles2;
            this->XfemElementInterface_partitionElement(& triangles, & subDomain1);
            this->XfemElementInterface_partitionElement(& triangles2, & subDomain2);

            int elIndex = element->giveGlobalNumber();


            ///////////////////////////////////////////
            // Write splitted elements to vtk.
            std::stringstream str1;
            str1 << "TriEl" << elIndex << "Side1.vtk";
            std::string name1 = str1.str();

            XFEMDebugTools::WriteTrianglesToVTK( name1, triangles );

            std::stringstream str2;
            str2 << "TriEl" << elIndex << "Side2.vtk";
            std::string name2 = str2.str();

            XFEMDebugTools::WriteTrianglesToVTK( name2, triangles2 );

            std::vector<Triangle> allTri;
//            for(int i = 1; i <= triangles2.giveSize(); i++)
//            {
//            	Triangle t(*(triangles2.at(i)));
//            	allTri.push_back( t );
//            }


            for ( int i = 1; i <= triangles2.giveSize(); i++ ) {
                int sz = triangles.giveSize();
//                triangles.put( sz + 1, triangles2.at(i) );
//                triangles2.unlink(i);

                triangles.put( sz + 1, new Triangle(*(triangles2.at(i)) ) );

            }


            for(int i = 1; i <= triangles.giveSize(); i++)
            {
            	Triangle t(*(triangles.at(i)));
            	allTri.push_back( t );
            }


            int ruleNum = 1;
//            int numGPPerTri = 3;
            int numGPPerTri = 19;
            AList< IntegrationRule >irlist;
            IntegrationRule *intRule = new PatchIntegrationRule(ruleNum, element, allTri);

            MaterialMode matMode = element->giveMaterialMode();
            intRule->SetUpPointsOnTriangle(numGPPerTri, matMode);
//            intRule->setUpIntegrationPoints(_Triangle, numGPPerTri, matMode);

            irlist.put(1, intRule);
            element->setIntegrationRules(&irlist);


        }
        else
        {
        	TrPlaneStress2d ::computeGaussPoints();

        }
*/
}

void TrPlaneStress2dXFEM :: XfemElementInterface_prepareNodesForDelaunay(std::vector< std::vector< FloatArray > > &oPointPartitions)
{
	XfemElementInterface :: XfemElementInterface_prepareNodesForDelaunay(oPointPartitions);
}
#endif

Interface *
TrPlaneStress2dXFEM :: giveInterface(InterfaceType it)
{
    if ( it == XfemElementInterfaceType ) {
        return ( XfemElementInterface * ) this;
    } else if ( it == VTKXMLExportModuleElementInterfaceType ) {
        return ( VTKXMLExportModuleElementInterface * ) this;
    } else {
        return TrPlaneStress2d :: giveInterface(it);
    }
}

int TrPlaneStress2dXFEM :: computeNumberOfDofs(EquationID ut)
{
    int ndofs = 0;

    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        ndofs += this->giveDofManager(i)->giveNumberOfDofs();
    }

    return ndofs;
}

void TrPlaneStress2dXFEM :: computeGaussPoints()
{
    XfemManager *xMan = this->giveDomain()->giveXfemManager();

    for(int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++)
    {
    	std::vector<FloatArray> intersecPoints;
    	EnrichmentItem *ei = xMan->giveEnrichmentItem(i);

        std::vector< int > intersecEdgeInd;
    	ei->computeIntersectionPoints(intersecPoints, intersecEdgeInd, this);
    	int numIntersecPoints = intersecPoints.size();

        if ( numIntersecPoints > 0 )
        {
            this->XfemElementInterface_updateIntegrationRule();
        } else
        {
            TrPlaneStress2d ::computeGaussPoints();
        }

    }

//	this->XfemElementInterface_updateIntegrationRule();
}


///@todo: Computation of N and B can be moved to the xfem interface and thus handle all continuum elements in the same way
void TrPlaneStress2dXFEM :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
	XfemElementInterface_createEnrBmatrixAt(answer, *gp, *this);
}


void TrPlaneStress2dXFEM :: computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
	XfemElementInterface_createEnrNmatrixAt(answer, *(gp->giveLocalCoordinates()), *this);
}



void
TrPlaneStress2dXFEM :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    // Returns the total id mask of the dof manager = regular id's + enriched id's

	XfemManager *xMan = this->domain->giveXfemManager();
	if( xMan != NULL )
	{
		this->giveDofManager(inode)->giveCompleteMasterDofIDArray(answer);
	}
	else
	{
	    // Continuous part
		TrPlaneStress2d ::giveDofManDofIDMask(inode, ut, answer);
	}



}
/*
void TrPlaneStress2dXFEM :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    this->computeStiffnessMatrix_withIRulesAsSubcells(answer, rMode, tStep);
}

void
TrPlaneStress2dXFEM :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    this->giveInternalForcesVector_withIRulesAsSubcells(answer, tStep, useUpdatedGpRecord);
}
*/

#ifdef __OOFEG
// TODO: FIX OOFEG implementation
void TrPlaneStress2dXFEM :: drawRawGeometry(oofegGraphicContext &context)
{
    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    XfemManager *xf = this->giveDomain()->giveXfemManager();
    if ( !xf->isElementEnriched(this) ) {
    	TrPlaneStress2d :: drawRawGeometry(context);
    } else {
        if ( numberOfIntegrationRules > 1 ) {
            int i;
//            PatchIntegrationRule *iRule;
            for ( i = 0; i < numberOfIntegrationRules; i++ ) {

            	// TODO: Implement visualization.
/*
            	iRule = dynamic_cast< PatchIntegrationRule * >( integrationRulesArray [ i ] );
                if ( iRule ) {
                    iRule->givePatch()->draw(context);
                }
*/
            }
        } else {
        	TrPlaneStress2d :: drawRawGeometry(context);
        }
    }
}

void TrPlaneStress2dXFEM :: drawScalar(oofegGraphicContext &context)
{
    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }
    XfemManager *xf = this->giveDomain()->giveXfemManager();
    if ( !xf->isElementEnriched(this) ) {
    	TrPlaneStress2d :: drawScalar(context);
    } else {
        if ( context.giveIntVarMode() == ISM_local ) {
            int indx;
            double val;
            FloatArray s(3), v;

            indx = context.giveIntVarIndx();

            TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
            PatchIntegrationRule *iRule;
            for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
                iRule = dynamic_cast< PatchIntegrationRule * >( integrationRulesArray [ i ] );

 #if 0
                val = iRule->giveMaterial();
 #else
                val = 0.0;
                for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
                    GaussPoint *gp = iRule->getIntegrationPoint(0);
                    giveIPValue(v, gp, context.giveIntVarType(), tStep);
                    val += v.at(indx);
                }

                val /= iRule->giveNumberOfIntegrationPoints();
 #endif
                s.at(1) = s.at(2) = s.at(3) = val;
            	// TODO: Implement visualization.
//                iRule->givePatch()->drawWD(context, s);

            }
        } else {
        	TrPlaneStress2d :: drawScalar(context);
        }
    }

}
#endif



int
TrPlaneStress2dXFEM :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode,
                                                                    TimeStep *stepN, const FloatArray &coords,
                                                                    FloatArray &answer)
{
	// TODO: Validate implementation.

    FloatArray lcoords, u;
    FloatMatrix n;
    int result;

    result = this->computeLocalCoordinates(lcoords, coords);

	XfemElementInterface_createEnrNmatrixAt(n, lcoords, *this);

    this->computeVectorOf(EID_MomentumBalance, mode, stepN, u);
    answer.beProductOf(n, u);

    return result;



/*
    FloatArray lcoords, u;
    FloatMatrix n;
    int result;

    result = this->computeLocalCoordinates(lcoords, coords);

    n.resize(2, 6);
    n.zero();

    n.at(1, 1) = lcoords.at(1);
    n.at(1, 3) = lcoords.at(2);
    n.at(1, 5) = lcoords.at(3);

    n.at(2, 2) = lcoords.at(1);
    n.at(2, 4) = lcoords.at(2);
    n.at(2, 6) = lcoords.at(3);

    this->computeVectorOf(EID_MomentumBalance, mode, stepN, u);
    answer.beProductOf(n, u);

    return result;
*/
}

void
TrPlaneStress2dXFEM :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
//    giveDofManDofIDMask(1, EID_MomentumBalance, answer);
	// TODO: For now, take only the continuous part
	int nodeInd = 1;
	TrPlaneStress2d ::giveDofManDofIDMask(nodeInd, EID_MomentumBalance, answer);
}


} /* namespace oofem */
