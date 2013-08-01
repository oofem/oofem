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

#include "mathfem.h"
#include "alist.h"
#include "enrichmentdomain.h"
#include "element.h"
#include "dofmanager.h"
#include "connectivitytable.h"
#include "classfactory.h"
#include "enrichmentfunction.h"
#include "xfemmanager.h"

#include <algorithm>

#include <cmath>

namespace oofem {

REGISTER_EnrichmentDomain( DofManList )
REGISTER_EnrichmentDomain( WholeDomain )
REGISTER_EnrichmentDomain( EDBGCircle )

#ifdef __BOOST_MODULE
REGISTER_EnrichmentDomain( EDCrack )
#endif
//REGISTER_EnrichmentDomain( BasicGeometryDomain<Line> )

// General 

EnrichmentDomain::EnrichmentDomain()
#ifdef __BOOST_MODULE
:
levelSetsNeedUpdate(true),
levelSetTol(1.0e-12), levelSetTol2(1.0e-12)
#endif
{

}


bool 
EnrichmentDomain :: isElementEnriched(Element *element) 
{
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        if ( this->isDofManagerEnriched( element->giveDofManager(i) ) ) {
            return true;
        }
    }
    return false;
}

bool EnrichmentDomain :: isAllElNodesEnriched(const Element *element)
{
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        if ( !this->isDofManagerEnriched( element->giveDofManager(i) ) ) {
            return false;
        }
    }
    return true;

}

#ifdef __BOOST_MODULE
void EnrichmentDomain :: updateLevelSets(XfemManager &ixFemMan)
{
	int nNodes = ixFemMan.giveDomain()->giveNumberOfDofManagers();

	levelSetPhi.resize(nNodes, 0.0);
	levelSetGamma.resize(nNodes, 0.0);


	int dim = ixFemMan.giveDomain()->giveNumberOfSpatialDimensions();

	if(dim == 2)
	{
		for ( int n = 1; n <= nNodes; n++ )
		{
			Node *node = ixFemMan.giveDomain()->giveNode(n);

			// Extract node coord
			const double &x = node->giveCoordinate(1);
			const double &y = node->giveCoordinate(2);

			// Calc normal sign dist
			double phi = 0.0;
			this->computeNormalSignDist(phi, x, y);
			levelSetPhi[n-1] = phi;

			// Calc tangential sign dist
			double gamma = 0.0;
			this->computeTangentialSignDist(gamma, x, y);
			levelSetGamma[n-1] = gamma;

		}
	}
	else
	{
		if(dim == 3)
		{
            OOFEM_ERROR2( "EnrichmentDomain :: updateLevelSets - Unsupported number of space dimensions: %d\n", dim );
		}
		else
		{
            OOFEM_ERROR2( "EnrichmentDomain :: updateLevelSets - Unsupported number of space dimensions: %d\n", dim );
		}
	}


	levelSetsNeedUpdate = false;


	updateNodeEnrMarker(ixFemMan);


}

void EnrichmentDomain :: updateNodeEnrMarker(XfemManager &ixFemMan)
{
	Domain *d = ixFemMan.giveDomain();
	int nEl = d->giveNumberOfElements();
	int nNodes = d->giveNumberOfDofManagers();

	nodeEnrichmentMarker.resize(nNodes, 0);

	// Loop over elements and use Phi and Gamma to determine if the element nodes are enriched.
	for(int elIndex = 1; elIndex <= nEl; elIndex++)
	{
		Element *el = d->giveElement(elIndex);
		int nElNodes = el->giveNumberOfNodes();

		int minSignPhi 	= 1, maxSignPhi 	= -1;
		int minSignGamma = 1, maxSignGamma = -1;

		for(int elNodeInd = 1; elNodeInd <= nElNodes; elNodeInd++)
		{
			int nGlob = el->giveNode(elNodeInd)->giveGlobalNumber();

			minSignPhi = min( bSign(minSignPhi), bSign(levelSetPhi[nGlob-1]) );
			maxSignPhi = max( bSign(maxSignPhi), bSign(levelSetPhi[nGlob-1]) );

			minSignGamma = min( bSign(minSignGamma), bSign( levelSetGamma[nGlob-1]) );
			maxSignGamma = max( bSign(maxSignGamma), bSign( levelSetGamma[nGlob-1]) );
		}


		if( minSignPhi*maxSignPhi < 0 && minSignGamma > 0 && maxSignGamma > 0 )
		{
			// Element completely cut by the crack
			// -> Apply step enrichment to all element nodes

			for(int elNodeInd = 1; elNodeInd <= nElNodes; elNodeInd++)
			{
				int nGlob = el->giveNode(elNodeInd)->giveGlobalNumber();

				if( nodeEnrichmentMarker[nGlob-1] == 0 )
				{
					nodeEnrichmentMarker[nGlob-1] = 1;
				}

			}


		}

		if( minSignPhi*maxSignPhi < 0 && minSignGamma*maxSignGamma < 0 )
		{
			// Element partly cut by the crack
			// -> Apply crack tip enrichment

		}


	}

}
#endif

void EnrichmentDomain_BG :: computeNormalSignDist(double &oDist, const double &iX, const double &iY)
{
	// TODO: Clean up interface
	FloatArray p; p.setValues(2, iX, iY);
	oDist = bg->computeDistanceTo(&p);
}

void EnrichmentDomain_BG :: computeTangentialSignDist(double &oDist, const double &iX, const double &iY)
{
	FloatArray p; p.setValues(2, iX, iY);
	oDist = bg->computeTangentialSignDist(&p);
}


// Node list

IRResultType DofManList :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro

    IntArray idList;
    IR_GIVE_FIELD(ir, idList, _IFT_DofManList_list);
    for ( int i = 1; i<=idList.giveSize(); i++) {
        this->dofManList.push_back( idList.at(i) );
    }
    return IRRT_OK;
    
}


bool DofManList :: isDofManagerEnriched(DofManager *dMan)
{
    int dManNumber = dMan->giveNumber();
    std::list< int > :: iterator p;
    p = std::find(this->dofManList.begin( ), this->dofManList.end( ), dManNumber);
    
    if ( p == this->dofManList.end( ) ) {
        return false;
    } else {
        return true;
    }
}




// Circle
bool 
EDBGCircle :: isDofManagerEnriched(DofManager *dMan)
{ 
#if 0
    // Only enrich the dofmans that are actually inside the domain
    FloatArray coords; 
    coords = *(dMan->giveCoordinates());
    return this->bg->isInside(coords);
#else
    // If any dofman of the neighboring elements is inside then the current dofman wil be enriched
    // => all dofmans of an element will be enriched if one dofman is inside.
    int node = dMan->giveGlobalNumber();
    Domain *d = dMan->giveDomain();
    Element *el;
    const IntArray *neighbours = d->giveConnectivityTable()->giveDofManConnectivityArray(node);
    for ( int i = 1; i <= neighbours->giveSize(); i++ ) {
        el = d->giveElement( neighbours->at(i) );
        for ( int j = 1; j <= el->giveNumberOfDofManagers(); j++ ) {
            if ( this->bg->isInside( * el->giveDofManager(j)->giveCoordinates() ) ) {
                return true;
            }
        }
    }

    return false;
#endif
};


bool
EDBGCircle :: isElementEnriched(Element *element) 
{
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        if ( this->isDofManagerEnriched( element->giveDofManager(i) ) ) {
            return true;
        }
    }
    return false;
};





#ifdef __BOOST_MODULE
// Crack
bool
EDCrack :: isDofManagerEnriched(DofManager *dMan)
{
	if( levelSetsNeedUpdate )
	{
		XfemManager *xMan = dMan->giveDomain()->giveXfemManager();
		updateLevelSets(*xMan);
	}

	// Use the the nodeEnrichmentMarker to check if the node is enriched.
	return ( nodeEnrichmentMarker[dMan->giveGlobalNumber()-1] != 0 );

};


bool
EDCrack :: isElementEnriched(Element *element)
{

    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        if ( this->isDofManagerEnriched( element->giveDofManager(i) ) ) {
            return true;
        }
    }

    return false;
};

void EDCrack :: computeIntersectionPoints(std::vector< FloatArray > &oIntersectionPoints, Element *element)
{

	if( isElementEnriched(element) )
	{
		// Use the level set functions to compute intersection points

		// Loop over element edges; an edge is intersected if the
		// node values of the level set functions have different signs

		int numEdges = element->giveNumberOfBoundarySides();

		for( int edgeIndex = 1; edgeIndex <= numEdges; edgeIndex++ )
		{
			IntArray bNodes;
			element->giveInterpolation()->boundaryGiveNodes(bNodes, edgeIndex);

			int nsLoc = bNodes.at( 1 );
			int nsGlob = element->giveNode(nsLoc)->giveGlobalNumber();
			int neLoc = bNodes.at( bNodes.giveSize() );
			int neGlob = element->giveNode(neLoc)->giveGlobalNumber();


			const double &phiS = levelSetPhi[nsGlob-1];
			const double &phiE = levelSetPhi[neGlob-1];

			if( phiS*phiE < levelSetTol2 )
			{
				// Intersection detected

				double xi = 0.0;

				if( fabs(phiS-phiE) > levelSetTol )
				{
					xi = (phiS+phiE)/(phiS-phiE);
				}

				if( xi < -1.0)
				{
					xi = -1.0;
				}

				if( xi > 1.0)
				{
					xi = 1.0;
				}


				FloatArray *ps = new FloatArray( *( element->giveDofManager(nsLoc)->giveCoordinates() ) );
				FloatArray *pe = new FloatArray( *( element->giveDofManager(neLoc)->giveCoordinates() ) );

				int nDim = ps->giveSize();
				FloatArray p;
				p.resize(nDim);

				for( int i = 1; i <= nDim; i++ )
				{
					(p.at(i)) = 0.5*(1.0-xi)*((ps->at(i))) + 0.5*(1.0+xi)*((pe->at(i)));
				}


				// Check that the intersection point has not already been identified.
				// This may happen if the crack intersects the element exactly at a node,
				// so that intersection is detected for both element edges in that node.

				bool alreadyFound = false;


				int numPointsOld = oIntersectionPoints.size();
				for(int k = 1; k <= numPointsOld; k++)
				{
					double dist = p.distance( oIntersectionPoints[k-1] );

					if( dist < levelSetTol )
					{
						alreadyFound = true;
						break;
					}

				}

				if(!alreadyFound)
				{
					oIntersectionPoints.push_back(p);
				}

			}

		}


	}

}

#endif // __BOOST_MODULE


} // end namespace oofem
