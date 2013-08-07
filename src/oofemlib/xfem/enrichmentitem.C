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

#include "xfemmanager.h"
#include "floatmatrix.h"
#include "enrichmentitem.h"
#include "element.h"
#include "enrichmentfunction.h"
#include "enrichmentdomain.h"
#include "cltypes.h"
#include "connectivitytable.h"
#include "oofem_limits.h"
#include "classfactory.h"
#include "mathfem.h"
#include "feinterpol.h"

#include <algorithm>

namespace oofem {

REGISTER_EnrichmentItem( CrackTip )
REGISTER_EnrichmentItem( CrackInterior )
REGISTER_EnrichmentItem( Inclusion )
REGISTER_EnrichmentItem( Delamination )

REGISTER_EnrichmentItem( Crack )

EnrichmentItem :: EnrichmentItem(int n, XfemManager *xMan, Domain *aDomain) : FEMComponent(n, aDomain),
mpEnrichmentDomain(NULL),
mEnrDomainIndex(0),
mpEnrichmentFunc(NULL),
mEnrFuncIndex(0),
mpEnrichmentFront(NULL),
mEnrFrontIndex(0),
mLevelSetsNeedUpdate(true),
mLevelSetTol(1.0e-12), mLevelSetTol2(1.0e-12)
{
    this->xMan = xMan;
    this->enrichmentFunctionList = new AList< EnrichmentFunction >(0);
    this->enrichmentDomainList = new AList< EnrichmentDomain >(0); 
    this->numberOfEnrichmentFunctions = 1;
    this->numberOfEnrichmentDomains = 1;
    this->startOfDofIdPool = -1;
    this->enrichesDofsWithIdArray = new IntArray;
}

EnrichmentItem :: ~EnrichmentItem()
{
    delete this->enrichmentFunctionList;
    delete this->enrichmentDomainList;

	delete this->enrichesDofsWithIdArray;

	if( mpEnrichmentDomain != NULL )
	{
		delete mpEnrichmentDomain;
		mpEnrichmentDomain = NULL;
	}

	if( mpEnrichmentFunc != NULL )
	{
		delete mpEnrichmentFunc;
		mpEnrichmentFunc = NULL;
	}

	if( mpEnrichmentFront != NULL)
	{
		delete mpEnrichmentFront;
		mpEnrichmentFront = NULL;
	}
}

IRResultType EnrichmentItem :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->enrichmentDomainNumbers, _IFT_EnrichmentItem_domains);
    this->numberOfEnrichmentDomains = this->enrichmentDomainNumbers.giveSize();

    IR_GIVE_OPTIONAL_FIELD(ir, enrichmentFunction, _IFT_EnrichmentItem_function);
    IR_GIVE_FIELD(ir, mEnrDomainIndex, _IFT_EnrichmentItem_domain);

    IR_GIVE_FIELD(ir, mEnrFrontIndex, _IFT_EnrichmentItem_front);

    mEnrFuncIndex = enrichmentFunction;

    
    return IRRT_OK;
}

int EnrichmentItem :: instanciateYourself(DataReader *dr)
{
    const char *__proc = "instanciateYourself"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro
    std :: string name;

    // Instanciate enrichment functions
    this->enrichmentFunctionList->growTo(numberOfEnrichmentFunctions);
    for ( int i = 1; i <= this->numberOfEnrichmentFunctions; i++ ) {
        InputRecord *mir = dr->giveInputRecord(DataReader :: IR_enrichFuncRec, i);
        result = mir->giveRecordKeywordField(name);

        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, "", mir, result);
        }
        
        EnrichmentFunction *ef = classFactory.createEnrichmentFunction( name.c_str(), i, this->xMan->giveDomain() );
        mpEnrichmentFunc = classFactory.createEnrichmentFunction( name.c_str(), i, this->xMan->giveDomain() );
        mpEnrichmentFunc->initializeFrom(mir);

        if ( ef == NULL ) {
            OOFEM_ERROR2( "EnrichmentItem::instanciateYourself: unknown enrichment function (%s)", name.c_str() );
        }
        enrichmentFunctionList->put(i, ef);
        ef->initializeFrom(mir);
    }
/*
    // Create new EnrichmentFunction
    printf("mEnrFuncIndex: %d\n", mEnrFuncIndex);
    InputRecord *irEnrFunc = dr->giveInputRecord(DataReader :: IR_enrichFuncRec, mEnrFuncIndex);
    std::string enrFuncName;
    result = irEnrFunc->giveRecordKeywordField(enrFuncName);
    mpEnrichmentFunc = classFactory.createEnrichmentFunction( enrFuncName.c_str(), mEnrFuncIndex, this->xMan->giveDomain() );

    if ( mpEnrichmentFunc == NULL ) {
        OOFEM_ERROR2( "EnrichmentItem::instanciateYourself: unknown enrichment function (%s)", enrFuncName.c_str() );
    }
    mpEnrichmentFunc->initializeFrom(irEnrFunc);
*/

    // Instanciate enrichment domains
    enrichmentDomainList->growTo(numberOfEnrichmentDomains);

    for ( int i = 1; i <= numberOfEnrichmentDomains; i++ ) {
        InputRecord *mir = dr->giveInputRecord(DataReader :: IR_geoRec, i);
        result = mir->giveRecordKeywordField(name);
        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, "", mir, result);
        }

        EnrichmentDomain *ed = classFactory.createEnrichmentDomain( name.c_str() ); 
        mpEnrichmentDomain = classFactory.createEnrichmentDomain( name.c_str() );
        if ( ed == NULL ) {
            OOFEM_ERROR2( "EnrichmentItem::instanciateYourself: unknown enrichment domain (%s)", name.c_str() );
        }
        this->enrichmentDomainList->put(i, ed);
        ed->initializeFrom(mir);
        mpEnrichmentDomain->initializeFrom(mir);

    }

    // TODO: Switch to factory pattern
    // Create EnrichmentFront
    printf("In EnrichmentItem :: instanciateYourself: mEnrFrontIndex: %d\n", mEnrFrontIndex );
    if(mEnrFrontIndex == 0)
    {
    	mpEnrichmentFront = new EnrFrontDoNothing();
    }

    if(mEnrFrontIndex == 1)
    {
    	mpEnrichmentFront = new EnrFrontExtend();
    }

/*
    // Create new EnrichmentDomain
    InputRecord *irEnrDom = dr->giveInputRecord(DataReader :: IR_geoRec, mEnrDomainIndex);
    result = irEnrDom->giveRecordKeywordField(name);
    mpEnrichmentDomain = classFactory.createEnrichmentDomain( name.c_str() );
    if ( mpEnrichmentDomain == NULL ) {
        OOFEM_ERROR2( "EnrichmentItem::instanciateYourself: unknown enrichment domain (%s)", name.c_str() );
    }
*/



    // Set start of the enrichment dof pool for the given EI
    int xDofPoolAllocSize = this->giveEnrichesDofsWithIdArray()->giveSize() * this->giveNumberOfEnrDofs() * this->giveNumberOfEnrichmentDomains(); 
    this->startOfDofIdPool = this->giveDomain()->giveNextFreeDofID(xDofPoolAllocSize);



    mpEnrichmentDomain->CallNodeEnrMarkerUpdate(*this, *xMan);

    return 1;
}

int
EnrichmentItem :: giveNumberOfEnrDofs()
{
    // returns the array of dofs a particular EI s
    int temp=0;
    for ( int i = 1; i <= this->giveNumberOfEnrichmentfunctions(); i++ ) {
        EnrichmentFunction *ef = this->giveEnrichmentFunction(i);
        // This is per regular dof
        temp += ef->giveNumberOfDofs(); // = number of functions associated with a particular enrichment function, e.g. 4 for branch function.

    }
    return temp;
}

EnrichmentFunction *EnrichmentItem :: giveEnrichmentFunction(int n)
{
    // Returns the n-th geometry.
    if ( enrichmentFunctionList->includes(n) ) {
        return enrichmentFunctionList->at(n);
    } else {
        OOFEM_ERROR2("giveEnrichmentFunction: undefined enrichment function (%d)", n);
    }

    return NULL;
}


bool EnrichmentItem :: isDofManEnriched(DofManager *dMan)
{
	return isDofManEnriched(*dMan);
}

bool EnrichmentItem :: isDofManEnrichedByEnrichmentDomain(DofManager *dMan, int edNumber)
{
	return isDofManEnriched(*dMan);
}

bool EnrichmentItem :: isElementEnriched(const Element *element)
{
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
    	if ( this->isDofManEnriched( *(element->giveDofManager(i)) ) ) {
            return true;
        }
    }
    return false;
}

bool EnrichmentItem :: isElementEnrichedByEnrichmentDomain(const Element *element, int edNumber)
{
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        DofManager *dMan = element->giveDofManager(i);
        if ( isDofManEnrichedByEnrichmentDomain(dMan, edNumber) ){
            return true;
        }
    }
    return false;
}

bool EnrichmentItem :: isDofManEnriched(const DofManager &iDMan) const
{
#if defined(ENABLE_XFEM_CPP11)
	auto begin 	= mEnrNodeIndices.begin();
	auto end 	= mEnrNodeIndices.end();
	int nodeInd	= iDMan.giveGlobalNumber();

	return std::binary_search(begin, end, nodeInd);
#else
	int nodeInd	= iDMan.giveGlobalNumber();
	return std::binary_search(mEnrNodeIndices.begin(), mEnrNodeIndices.end(), nodeInd);
#endif
}

int  EnrichmentItem :: giveNumDofManEnrichments(const DofManager &iDMan) const
{
#if defined(ENABLE_XFEM_CPP11)
	auto begin 	= mEnrNodeIndices.begin();
	auto end 	= mEnrNodeIndices.end();
	int nodeInd	= iDMan.giveGlobalNumber();

	auto it 	= std::find(begin, end, nodeInd);

	if( it != end )
	{
		// TODO: Check the actual number of enrichments.
		return 1;
	}
#else
	std::vector<int>::const_iterator begin = mEnrNodeIndices.begin();
	std::vector<int>::const_iterator end 	= mEnrNodeIndices.end();
	int nodeInd	= iDMan.giveGlobalNumber();

	std::vector<int>::const_iterator it = std::find(begin, end, nodeInd);

	if( it != end )
	{
		// TODO: Check the actual number of enrichments.
		return 1;
	}

#endif

	return 0;
}

bool EnrichmentItem :: isMaterialModified(GaussPoint &iGP, Element &iEl, StructuralMaterial *&opSM) const
{
	return false;
}

void EnrichmentItem :: updateGeometry()
{
    mpEnrichmentDomain->CallNodeEnrMarkerUpdate(*this, *xMan);
}


void
EnrichmentItem :: computeDofManDofIdArray(IntArray &answer, DofManager *dMan, int enrichmentDomainNumber)
{
    // Gives an array containing the dofId's that should be created as new dofs (what dofs to enrich).
    IntArray *enrichesDofsWithIdArray = this->giveEnrichesDofsWithIdArray();
    int eiEnrSize = enrichesDofsWithIdArray->giveSize();

    // Go through the list of dofs that the EI supports and compare with the available dofs in the dofMan.
    // Store matches in dofMask
    IntArray dofMask(eiEnrSize); dofMask.zero();
    int count = 0;
    for ( int i = 1; i <= eiEnrSize; i++ ) {
        if ( dMan->hasDofID( (DofIDItem) enrichesDofsWithIdArray->at(i) ) ) {
            count++;
            dofMask.at(count) = dMan->giveDofWithID( enrichesDofsWithIdArray->at(i) )->giveNumber();
        }
    }

    answer.resize(count);
    int xDofAllocSize = eiEnrSize * this->giveNumberOfEnrDofs(); // number of new dof id's the ei will allocate
    for ( int i = 1; i <= count; i++ ) {
        answer.at(i) = this->giveStartOfDofIdPool() + (enrichmentDomainNumber-1)*xDofAllocSize + dofMask.at(i)-1 ;
    }

}

void
EnrichmentItem :: giveEIDofIdArray(IntArray &answer, int enrichmentDomainNumber)
{
    // Returns an array containing the dof Id's of the new enrichment dofs pertinent to the ei.
    // Note: the dof managers may not support these dofsall potential dof id's
    IntArray *enrichesDofsWithIdArray = this->giveEnrichesDofsWithIdArray();
    int eiEnrSize = enrichesDofsWithIdArray->giveSize();

    answer.resize(eiEnrSize);
    int xDofAllocSize = eiEnrSize * this->giveNumberOfEnrDofs(); // number of new dof id's the ei will allocate
    for ( int i = 1; i <= eiEnrSize; i++ ) {
        answer.at(i) = this->giveStartOfDofIdPool() + (enrichmentDomainNumber-1)*xDofAllocSize + (i-1);
    }
}

void EnrichmentItem :: evaluateEnrFuncAt(double &oEnrFunc, const FloatArray &iPos, const double &iLevelSet) const
{
	mpEnrichmentFunc->evaluateEnrFuncAt(oEnrFunc, iPos, iLevelSet, mpEnrichmentDomain);
}

void EnrichmentItem :: evaluateEnrFuncDerivAt(FloatArray &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet) const
{
	mpEnrichmentFunc->evaluateEnrFuncDerivAt(oEnrFuncDeriv, iPos, iLevelSet, iGradLevelSet, mpEnrichmentDomain);
}

void EnrichmentItem :: updateLevelSets(XfemManager &ixFemMan)
{

	int nNodes = ixFemMan.giveDomain()->giveNumberOfDofManagers();

	mLevelSetNormalDir.resize(nNodes, 0.0);
	mLevelSetTangDir.resize(nNodes, 0.0);


	for ( int n = 1; n <= nNodes; n++ )
	{
		Node *node = ixFemMan.giveDomain()->giveNode(n);

		// Extract node coord
		const FloatArray &pos( *node->giveCoordinates() );

		// Calc normal sign dist
		double phi = 0.0;
		mpEnrichmentDomain->computeNormalSignDist(phi, pos);
		mLevelSetNormalDir[n-1] = phi;

		// Calc tangential sign dist
		double gamma = 0.0;
		mpEnrichmentDomain->computeTangentialSignDist(gamma, pos);
		mLevelSetTangDir[n-1] = gamma;

	}

	mLevelSetsNeedUpdate = false;

}

void EnrichmentItem :: updateNodeEnrMarker(XfemManager &ixFemMan, const EnrichmentDomain_BG &iEnrichmentDomain_BG)
{
	updateLevelSets(ixFemMan);

	Domain *d = ixFemMan.giveDomain();
	int nEl = d->giveNumberOfElements();
	int nNodes = d->giveNumberOfDofManagers();

	mNodeEnrMarker.resize(nNodes, 0);


	// Loop over elements and use the level sets to mark nodes belonging to completely cut elements.
	for(int elIndex = 1; elIndex <= nEl; elIndex++)
	{
		Element *el = d->giveElement(elIndex);
		int nElNodes = el->giveNumberOfNodes();

		int minSignPhi 	= 1, maxSignPhi 	= -1;
		int minSignGamma = 1, maxSignGamma = -1;

		for(int elNodeInd = 1; elNodeInd <= nElNodes; elNodeInd++)
		{
			int nGlob = el->giveNode(elNodeInd)->giveGlobalNumber();

			minSignPhi = std::min( sgn(minSignPhi), sgn(mLevelSetNormalDir[nGlob-1]) );
			maxSignPhi = std::max( sgn(maxSignPhi), sgn(mLevelSetNormalDir[nGlob-1]) );

			minSignGamma = std::min( sgn(minSignGamma), sgn( mLevelSetTangDir[nGlob-1]) );
			maxSignGamma = std::max( sgn(maxSignGamma), sgn( mLevelSetTangDir[nGlob-1]) );
		}


		if( minSignPhi*maxSignPhi < 0 && minSignGamma > 0 && maxSignGamma > 0 )
		{
			// Element completely cut by the crack
			// -> Apply step enrichment to all element nodes

			for(int elNodeInd = 1; elNodeInd <= nElNodes; elNodeInd++)
			{
				int nGlob = el->giveNode(elNodeInd)->giveGlobalNumber();

				if( mNodeEnrMarker[nGlob-1] == 0 )
				{
					mNodeEnrMarker[nGlob-1] = 1;
				}

			}


		}
/*
		if( minSignPhi*maxSignPhi < 0 && minSignGamma*maxSignGamma < 0 )
		{
			if(mMarkTipNodes)
			{
				// Check if the element is intersected by the crack
				// (This check is actually needed!)

				bool edgeIntersected = false;
				for(int elNodeInd = 1; elNodeInd <= nElNodes-1; elNodeInd++)
				{
					int niGlob = el->giveNode(elNodeInd  )->giveGlobalNumber();
					int njGlob = el->giveNode(elNodeInd+1)->giveGlobalNumber();

					if( levelSetPhi[niGlob-1]*levelSetPhi[njGlob-1] < 0.0 )
					{
						if( levelSetGamma[niGlob-1] > 0.0 && levelSetGamma[njGlob-1] > 0.0 )
						{
							edgeIntersected = true;
						}
					}

				}

				int niGlob = el->giveNode(1  		)->giveGlobalNumber();
				int njGlob = el->giveNode(nElNodes	)->giveGlobalNumber();

				if( levelSetPhi[niGlob-1]*levelSetPhi[njGlob-1] < 0.0 )
				{
					if( levelSetGamma[niGlob-1] > 0.0 && levelSetGamma[njGlob-1] > 0.0 )
					{
						edgeIntersected = true;
					}
				}

				if( edgeIntersected )
				{

						printf("Inserting index %d\n", elIndex );
						tipEls.push_back(elIndex);
				}
			}
		}
*/

/*
		if( minSignPhi*maxSignPhi < 0 && minSignGamma*maxSignGamma < 0 )
		{
			// Element partly cut by the crack
			// -> Apply crack tip enrichment

			if(mMarkTipNodes && false)
			{
				// Check if the element is intersected by the crack
				// (This check is actually needed!)

				bool edgeIntersected = false;
				for(int elNodeInd = 1; elNodeInd <= nElNodes-1; elNodeInd++)
				{
					int niGlob = el->giveNode(elNodeInd  )->giveGlobalNumber();
					int njGlob = el->giveNode(elNodeInd+1)->giveGlobalNumber();

					if( levelSetPhi[niGlob-1]*levelSetPhi[njGlob-1] < 0.0 )
					{
						if( levelSetGamma[niGlob-1] > 0.0 && levelSetGamma[njGlob-1] > 0.0 )
						{
							edgeIntersected = true;
						}
					}

				}

				int niGlob = el->giveNode(1  		)->giveGlobalNumber();
				int njGlob = el->giveNode(nElNodes	)->giveGlobalNumber();

				if( levelSetPhi[niGlob-1]*levelSetPhi[njGlob-1] < 0.0 )
				{
					if( levelSetGamma[niGlob-1] > 0.0 && levelSetGamma[njGlob-1] > 0.0 )
					{
						edgeIntersected = true;
					}
				}

				if( edgeIntersected )
				{
					for(int elNodeInd = 1; elNodeInd <= nElNodes; elNodeInd++)
					{
						int nGlob = el->giveNode(elNodeInd)->giveGlobalNumber();

//						if( nodeEnrichmentMarker[nGlob-1] == 0 )
//						{
							nodeEnrichmentMarker[nGlob-1] = 2;
//						}

					}
				}

			}

		}
*/


	}

	// Mark tip nodes for special treatment.
	mpEnrichmentFront->MarkNodesAsFront(mNodeEnrMarker, *xMan, mLevelSetNormalDir, mLevelSetTangDir);

/*

	// Loop over the tip elements and determine which nodes to enrich.
	for(int i = 0; i < tipEls.size(); i++)
	{
		int elIndex = tipEls[i];
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
*/

/*
		// Element partly cut by the crack
		// -> Apply crack tip enrichment

		if(mMarkTipNodes)
		{

			// Decide where to remove enrichment
			bool foundFrontNodes = false;
			for(int elNodeInd = 1; elNodeInd <= nElNodes-1; elNodeInd++)
			{
				int niGlob = el->giveNode(elNodeInd  )->giveGlobalNumber();
				int njGlob = el->giveNode(elNodeInd+1)->giveGlobalNumber();

				if( levelSetPhi[niGlob-1]*levelSetPhi[njGlob-1] < 0.0 )
				{

					if( (nodeEnrichmentMarker[niGlob-1] == 0 && nodeEnrichmentMarker[njGlob-1] == 1) || (nodeEnrichmentMarker[niGlob-1] == 1 && nodeEnrichmentMarker[njGlob-1] == 0 ) )
					{
						nodeEnrichmentMarker[niGlob-1] = 0;
						nodeEnrichmentMarker[njGlob-1] = 0;

						foundFrontNodes = true;
						break;
					}
				}

			}

			if(!foundFrontNodes)
			{
				int niGlob = el->giveNode(1  		)->giveGlobalNumber();
				int njGlob = el->giveNode(nElNodes	)->giveGlobalNumber();

				if( levelSetPhi[niGlob-1]*levelSetPhi[njGlob-1] < 0.0 )
				{

					if( (nodeEnrichmentMarker[niGlob-1] == 0 && nodeEnrichmentMarker[njGlob-1] == 1) || (nodeEnrichmentMarker[niGlob-1] == 1 && nodeEnrichmentMarker[njGlob-1] == 0 ) )
					{
						nodeEnrichmentMarker[niGlob-1] = 0;
						nodeEnrichmentMarker[njGlob-1] = 0;

						foundFrontNodes = true;
					}
				}
			}




		}



	}
*/

/*
	// Loop over elements once more, this time to mark crack tip nodes.
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



		if( minSignPhi*maxSignPhi < 0 && minSignGamma*maxSignGamma < 0 )
		{
			// Element partly cut by the crack
			// -> Apply crack tip enrichment

			if(mMarkTipNodes)
			{

				// Decide where to remove enrichment
				bool foundFrontNodes = false;
				for(int elNodeInd = 1; elNodeInd <= nElNodes-1; elNodeInd++)
				{
					int niGlob = el->giveNode(elNodeInd  )->giveGlobalNumber();
					int njGlob = el->giveNode(elNodeInd+1)->giveGlobalNumber();

					if( levelSetPhi[niGlob-1]*levelSetPhi[njGlob-1] < 0.0 )
					{

						if( (nodeEnrichmentMarker[niGlob-1] == 0 && nodeEnrichmentMarker[njGlob-1] == 1) || (nodeEnrichmentMarker[niGlob-1] == 1 && nodeEnrichmentMarker[njGlob-1] == 0 ) )
						{
							nodeEnrichmentMarker[niGlob-1] = 0;
							nodeEnrichmentMarker[njGlob-1] = 0;

							foundFrontNodes = true;
							break;
						}
					}

				}

				if(!foundFrontNodes)
				{
					int niGlob = el->giveNode(1  		)->giveGlobalNumber();
					int njGlob = el->giveNode(nElNodes	)->giveGlobalNumber();

					if( levelSetPhi[niGlob-1]*levelSetPhi[njGlob-1] < 0.0 )
					{

						if( (nodeEnrichmentMarker[niGlob-1] == 0 && nodeEnrichmentMarker[njGlob-1] == 1) || (nodeEnrichmentMarker[niGlob-1] == 1 && nodeEnrichmentMarker[njGlob-1] == 0 ) )
						{
							nodeEnrichmentMarker[niGlob-1] = 0;
							nodeEnrichmentMarker[njGlob-1] = 0;

							foundFrontNodes = true;
						}
					}
				}




			}

		}


	}
*/


	mEnrNodeIndices.clear();

	// Loop over nodes and add the indices of enriched nodes.
	// Since we loop over the nodes in order from 1 to nNodes,
	// mEnrNodeIndices will automatically be sorted.
	for( int i = 1; i <= nNodes; i++)
	{
		if( mNodeEnrMarker[i-1] > 0 )
		{
			mEnrNodeIndices.push_back(i);
		}
	}

}

void EnrichmentItem :: updateNodeEnrMarker(XfemManager &ixFemMan, const DofManList &iEnrichmentDomain_BG)
{

}

void EnrichmentItem :: updateNodeEnrMarker(XfemManager &ixFemMan, const WholeDomain &iEnrichmentDomain_BG)
{

}

void EnrichmentItem :: computeIntersectionPoints(std::vector< FloatArray > &oIntersectionPoints, Element *element)
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


			const double &phiS = mLevelSetNormalDir[nsGlob-1];
			const double &phiE = mLevelSetNormalDir[neGlob-1];

			if( phiS*phiE < mLevelSetTol2 )
			{
				// Intersection detected

				double xi = 0.0;

				if( fabs(phiS-phiE) > mLevelSetTol )
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


				FloatArray ps( *( element->giveDofManager(nsLoc)->giveCoordinates() ) );
				FloatArray pe( *( element->giveDofManager(neLoc)->giveCoordinates() ) );

				int nDim = ps.giveSize();
				FloatArray p;
				p.resize(nDim);

				for( int i = 1; i <= nDim; i++ )
				{
					(p.at(i)) = 0.5*(1.0-xi)*((ps.at(i))) + 0.5*(1.0+xi)*((pe.at(i)));
				}


				// Check that the intersection point has not already been identified.
				// This may happen if the crack intersects the element exactly at a node,
				// so that intersection is detected for both element edges in that node.

				bool alreadyFound = false;


				int numPointsOld = oIntersectionPoints.size();
				for(int k = 1; k <= numPointsOld; k++)
				{
					double dist = p.distance( oIntersectionPoints[k-1] );

					if( dist < mLevelSetTol )
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


Inclusion :: Inclusion(int n, XfemManager *xm, Domain *aDomain) :
EnrichmentItem(n, xm, aDomain),
mat(NULL)
{ 
    this->enrichesDofsWithIdArray->setValues(3, D_u, D_v, D_w);
}

Inclusion :: ~Inclusion()
{
	if(mat != NULL)
	{
		mat = NULL;
	}
}

bool Inclusion :: isMaterialModified(GaussPoint &iGP, Element &iEl, StructuralMaterial *&opSM) const
{
	// Check if the point is located inside the inclusion

    FloatArray N;
    FEInterpolation *interp = iEl.giveInterpolation();
    interp->evalN(     N  , *iGP.giveCoordinates(), FEIElementGeometryWrapper(&iEl) );

    const IntArray &elNodes = iEl.giveDofManArray();

	double levelSetGP = 0.0;
	this->interpLevelSet(levelSetGP, N, elNodes);

	if(levelSetGP < 0.0)
	{
		opSM = static_cast< StructuralMaterial * > (mat);
		return true;
	}

	return false;
}

IRResultType Inclusion :: initializeFrom(InputRecord *ir)
{
    EnrichmentItem :: initializeFrom(ir);
    const char *__proc = "initializeFrom";
    IRResultType result;
    int material = 0;
    IR_GIVE_FIELD(ir, material, _IFT_Inclusion_material);
    this->mat = this->giveDomain()->giveMaterial(material);

    this->numberOfEnrichmentFunctions = 1;
    // Not sure this should be input at but instead be determined by the ei which describes the physical model /JB
    //IR_GIVE_OPTIONAL_FIELD(ir, numberOfEnrichmentFunctions, _IFT_XfemManager_numberOfEnrichmentFunctions, "numberofenrichmentfunctions");
    return IRRT_OK;
}



// DELAMINATION

Delamination :: Delamination(int n, XfemManager *xm, Domain *aDomain) : EnrichmentItem(n, xm, aDomain)
{ 
    //this->enrichesDofsWithIdArray->setValues(7, D_u, D_v, D_w, W_u, W_v, W_w, Gamma);
    this->enrichesDofsWithIdArray->setValues(6, D_u, D_v, D_w, W_u, W_v, W_w);
}


IRResultType Delamination :: initializeFrom(InputRecord *ir)
{
    this->numberOfEnrichmentFunctions = 1; // must be set before EnrichmentItem :: initializeFrom(ir) is called
    EnrichmentItem :: initializeFrom(ir);
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->enrichmentDomainXiCoords, _IFT_Delamination_xiCoords);
    if ( this->numberOfEnrichmentDomains != this->enrichmentDomainXiCoords.giveSize() ) {
        OOFEM_ERROR3( "EnrichmentItem :: initializeFrom: size of enrichmentDomainXiCoords (%i) differs from numberOfEnrichmentDomains (%i)", 
                       this->enrichmentDomainXiCoords.giveSize(), this->numberOfEnrichmentDomains );
    }

    //write an instanciate method

    return IRRT_OK;
}




double 
Delamination :: giveDelaminationZCoord(int n, Element *element) 
{
    //AList<double> *xiCoordList;
    int nDelam = this->giveNumberOfEnrichmentDomains(); // max possible number
    int pos = 1;
    for ( int i = 1; i <= nDelam; i++ ) {
        if( this->isElementEnriched(element) ) {
            //xiCoordList.put( pos, this->delaminationZCoords.at(i) );
            pos++;
        } 
    }
    return 0.;
};

// Remove!
int
Delamination :: giveDelaminationGroupAt(double zeta) 
{
    //double zRef = Shell7Base :: giveLocalZetaCoord(gp);
    int nDelam = this->giveNumberOfEnrichmentDomains();
    for ( int j = 1; j <= nDelam; j++ ) {
        //double zDelam = this->giveDelaminationZCoord(j);
        double zDelam = 0.;
        if ( zeta  < zDelam ) { //belong to the delamination group just below delamination #j. How to deal with poins that lie onthe boundary?
            return j;
        }

    }
    return nDelam + 1;
}

double 
Delamination :: heaviside(double xi, double xi0)
{
    if( xi < xi0 ) {
        return 0.0;
    } else {
        return 1.0;
    }
}


double 
Delamination :: giveDelaminationGroupMidZ(int dGroup, Element *e)
{
    double zTop=0., zBottom=0.;
    this->giveDelaminationGroupZLimits(dGroup, zTop, zBottom, e);
    return 0.5 * ( zTop + zBottom );
}


double 
Delamination :: giveDelaminationGroupThickness(int dGroup, Element *e)
{
    double zTop, zBottom;
    this->giveDelaminationGroupZLimits(dGroup, zTop, zBottom, e);
    return zTop - zBottom;
}

void 
Delamination :: giveDelaminationGroupZLimits(int &dGroup, double &zTop, double &zBottom, Element *e)
{
    int nDelam = this->giveNumberOfEnrichmentDomains();
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * > (e->giveCrossSection());
    
    if ( dGroup == 1 ) {
        zBottom = - layeredCS->giveMidSurfaceZcoordFromBottom();
        zTop    =  0.;//this->giveDelaminationZCoord(dGroup);
    } else if (dGroup == nDelam + 1) {
        zBottom =  0.;//this->giveDelaminationZCoord(dGroup-1);
        zTop    = -layeredCS->giveMidSurfaceZcoordFromBottom() + layeredCS->computeIntegralThick();
    } else {
        zBottom =  0.;//this->giveDelaminationZCoord(dGroup-1);
        zTop    =  0.;//this->giveDelaminationZCoord(dGroup);
    }
#ifdef DEBUG
    if ( zBottom > zTop ) {
        OOFEM_ERROR2("giveDelaminationGroupZLimits: Bottom z-coord is larger than top z-coord in dGroup. (%i)", dGroup);
    }
#endif
}


Crack :: Crack(int n, XfemManager *xm, Domain *aDomain) : EnrichmentItem(n, xm, aDomain)
{
    this->enrichesDofsWithIdArray->setValues(3, D_u, D_v, D_w);
}

IRResultType Crack :: initializeFrom(InputRecord *ir)
{
    this->numberOfEnrichmentFunctions = 1;
    EnrichmentItem :: initializeFrom(ir);

    return IRRT_OK;
}

void EnrFrontExtend :: MarkNodesAsFront(std::vector<int> &ioNodeEnrMarker, XfemManager &ixFemMan, const std::vector<double> &iLevelSetNormalDir, const std::vector<double> &iLevelSetTangDir)
{

	// Extend the set of enriched nodes as follows:
	// If any node of the neighboring elements is enriched, the current node is also enriched.

	Domain &d = *(ixFemMan.giveDomain());

	// Loop over all nodes
	int nNodes = d.giveNumberOfDofManagers();

	std::vector<int> newEnrNodes;
	for(int i = 1; i <= nNodes; i++)
	{
		// Check if the node is already enriched
		bool alreadyEnr = (ioNodeEnrMarker[i-1] > 0);

#if defined(ENABLE_XFEM_CPP11)
		auto begin 	= newEnrNodes.begin();
		auto end 	= newEnrNodes.end();
#else
		std::vector<int>::const_iterator begin 	= newEnrNodes.begin();
		std::vector<int>::const_iterator end 	= newEnrNodes.end();
#endif
		if( std::binary_search(begin, end, i) )
		{
			alreadyEnr = true;
		}


		if( !alreadyEnr )
		{
			bool goOn = true;

			// Loop over neighbors
			const IntArray &neigh = *(d.giveConnectivityTable()->giveDofManConnectivityArray(i) );
			for ( int j = 1; j <= neigh.giveSize(); j++ )
			{
				if(!goOn)
				{
					break;
				}

				Element &el = *(d.giveElement( neigh.at(j) ));

				// Loop over neighbor element nodes
				for ( int k = 1; k <= el.giveNumberOfDofManagers(); k++ ) {

					int kGlob = el.giveDofManager(k)->giveGlobalNumber();
					if( iLevelSetNormalDir[kGlob-1] < 0 )
					{
						newEnrNodes.push_back(i);
						goOn = false;
						break;
					}

				}

			}

		}

	}


	// Mark the new nodes to be enriched
	for(int i = 0; i < newEnrNodes.size(); i++)
	{
		ioNodeEnrMarker[ newEnrNodes[i]-1 ] = 1;
	}
}

} // end namespace oofem
