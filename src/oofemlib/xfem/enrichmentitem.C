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
#include "fracturemanager.h"
#include "mathfem.h"
#include "feinterpol.h"
#include "masterdof.h"
#include "propagationlaw.h"
#include <algorithm>
#include <limits>

namespace oofem {
REGISTER_EnrichmentItem(Inclusion)
REGISTER_EnrichmentItem(Delamination)

REGISTER_EnrichmentItem(Crack)

EnrichmentItem :: EnrichmentItem(int n, XfemManager *xMan, Domain *aDomain) : FEMComponent(n, aDomain),
    mpEnrichmentDomain(NULL),
    mEnrDomainIndex(0),
    mpEnrichmentFunc(NULL),
    mEnrFuncIndex(0),
    mpEnrichmentFront(NULL),
    mEnrFrontIndex(0),
    mpPropagationLaw(NULL),
    mPropLawIndex(0),
    mLevelSetsNeedUpdate(true),
    mLevelSetTol(1.0e-12), mLevelSetTol2(1.0e-12)
{
    this->xMan = xMan;
    this->enrichmentFunctionList = new AList< EnrichmentFunction >(0);
    this->enrichmentDomainList = new AList< EnrichmentDomain >(0);
    this->numberOfEnrichmentFunctions = 1;
    this->numberOfEnrichmentDomains = 1;
    this->startOfDofIdPool = -1;
    this->endOfDofIdPool = -1;
    this->mpEnrichesDofsWithIdArray = new IntArray;
}

EnrichmentItem :: ~EnrichmentItem()
{
    delete this->enrichmentFunctionList;
    delete this->enrichmentDomainList;

    delete this->mpEnrichesDofsWithIdArray;

    if ( mpEnrichmentDomain != NULL ) {
        delete mpEnrichmentDomain;
        mpEnrichmentDomain = NULL;
    }

    if ( mpEnrichmentFunc != NULL ) {
        delete mpEnrichmentFunc;
        mpEnrichmentFunc = NULL;
    }

    if ( mpEnrichmentFront != NULL ) {
        delete mpEnrichmentFront;
        mpEnrichmentFront = NULL;
    }

    if ( mpPropagationLaw != NULL ) {
        delete mpPropagationLaw;
        mpPropagationLaw = NULL;


    }

}

bool EnrichmentItem :: isElementFullyEnrichedByEnrichmentDomain(const Element *element, int edNumber) 
{
    // Checks if all of the dofmanagers are enriched
    int count = 0;
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        DofManager *dMan = element->giveDofManager(i);
        if ( isDofManEnrichedByEnrichmentDomain(dMan, edNumber) ){
            count++;
        }
    }
    return count == element->giveNumberOfDofManagers();
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


    IR_GIVE_OPTIONAL_FIELD(ir, mPropLawIndex, _IFT_EnrichmentItem_propagationlaw);


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
        if(mpEnrichmentFunc != NULL) {
        	mpEnrichmentFunc->initializeFrom(mir);
        }
        else {
            OOFEM_ERROR2( "EnrichmentItem::instanciateYourself: failed to create enrichment function (%s)", name.c_str() );
        }

        if ( ef == NULL ) {
            OOFEM_ERROR2( "EnrichmentItem::instanciateYourself: unknown enrichment function (%s)", name.c_str() );
        }

        enrichmentFunctionList->put(i, ef);
        ef->initializeFrom(mir);
    }


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


    // Instantiate EnrichmentFront
    std::string enrFrontName;

    InputRecord *enrFrontir = dr->giveInputRecord(DataReader :: IR_enrichFrontRec, mEnrFrontIndex);
    result = enrFrontir->giveRecordKeywordField(enrFrontName);

    mpEnrichmentFront = classFactory.createEnrichmentFront( enrFrontName.c_str() );
    if(mpEnrichmentFront != NULL) {
    	mpEnrichmentFront->initializeFrom(enrFrontir);
    }
    else {
        OOFEM_ERROR2( "EnrichmentItem::instanciateYourself: Failed to create enrichment front (%s)", enrFrontName.c_str() );
    }



    // Instantiate PropagationLaw
    std::string propLawName;
    
    InputRecord *propLawir = dr->giveInputRecord(DataReader :: IR_propagationLawRec, mPropLawIndex);
    result = propLawir->giveRecordKeywordField(propLawName);

    if( mPropLawIndex == 0 ) {
    	// Dummy propagation law
    	printf("Creating dummy propagation law.\n");
    	mpPropagationLaw = new PLDoNothing();
    }
    else {
    	// Propagation law from input record
    	printf("Creating propagation law from input record. propLawName.c_str(): %s \n", propLawName.c_str() );
		mpPropagationLaw = classFactory.createPropagationLaw( propLawName.c_str() );
		if(mpPropagationLaw != NULL) {
			mpPropagationLaw->initializeFrom(propLawir);
		}
		else {
			OOFEM_ERROR2( "EnrichmentItem::instanciateYourself: Failed to create propagation law (%s)", propLawName.c_str() );
		}

    }

    // Set start of the enrichment dof pool for the given EI
//    int xDofPoolAllocSize = this->giveEnrichesDofsWithIdArray()->giveSize() * this->giveNumberOfEnrDofs() * this->giveNumberOfEnrichmentDomains(); 
    int xDofPoolAllocSize = this->giveEnrichesDofsWithIdArray()->giveSize() * this->giveNumberOfEnrDofs() * 1; 
    this->startOfDofIdPool = this->giveDomain()->giveNextFreeDofID(xDofPoolAllocSize);
    this->endOfDofIdPool = this->startOfDofIdPool + xDofPoolAllocSize - 1;


    mpEnrichmentDomain->CallNodeEnrMarkerUpdate(* this, * xMan);

    return 1;
}

int
EnrichmentItem :: giveNumberOfEnrDofs() const
{
	// TODO: Take branch functions into account when computing the total number of dofs.
    // returns the array of dofs a particular EI s

	int numEnrDofs = mpEnrichmentFunc->giveNumberOfDofs();

	if(mpEnrichmentFront != NULL) {
		numEnrDofs = max(numEnrDofs, mpEnrichmentFront->giveMaxNumEnrichments() );
	}

	return numEnrDofs;
}

bool EnrichmentItem :: isDofManEnrichedByEnrichmentDomain(DofManager *dMan, int edNumber) const
{
    EnrichmentDomain *ed = this->enrichmentDomainList->at(edNumber);
    //return ed->isDofManagerEnriched(dMan);
    return this->isDofManEnriched(*dMan);

}

//{
//    return isDofManEnriched(* dMan);
//}
    // Note: the dof managers may not support all/any of these new potential dof id's. Eg. a 
    // beam will not be able to account for a pressure dof. 

bool EnrichmentItem :: isElementEnriched(const Element *element) const
{
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        if ( this->isDofManEnriched( * ( element->giveDofManager(i) ) ) ) {
            return true;
        }
    }

    return false;
}

/*
bool EnrichmentItem :: isElementEnrichedByEnrichmentDomain(const Element *element, int edNumber) const
{
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        DofManager *dMan = element->giveDofManager(i);
        if ( isDofManEnrichedByEnrichmentDomain(dMan, edNumber) ) {
            return true;
        }
    }

    return false;
}
*/


bool EnrichmentItem :: isDofManEnriched(const DofManager &iDMan) const
{
#if defined( ENABLE_XFEM_CPP11 )
    auto begin      = mEnrNodeIndices.begin();
    auto end        = mEnrNodeIndices.end();
    int nodeInd     = iDMan.giveGlobalNumber();

    return std :: binary_search(begin, end, nodeInd);

#else
    int nodeInd     = iDMan.giveGlobalNumber();
    return std :: binary_search(mEnrNodeIndices.begin(), mEnrNodeIndices.end(), nodeInd);

#endif
}

int EnrichmentItem :: giveNumDofManEnrichments(const DofManager &iDMan) const
{
#if defined( ENABLE_XFEM_CPP11 )
    auto begin      = mEnrNodeIndices.begin();
    auto end        = mEnrNodeIndices.end();
    int nodeInd     = iDMan.giveGlobalNumber();

    auto it         = std :: find(begin, end, nodeInd);

    if ( it != end ) {
        int enrichmentType = mNodeEnrMarker [ * it - 1 ];

        if ( enrichmentType == 1 ) {
            // Bulk enrichment
            return 1;
        } else   {
            // Front enrichment
            return mpEnrichmentFront->giveNumEnrichments();
        }
    }

#else
    std :: vector< int > :: const_iterator begin = mEnrNodeIndices.begin();
    std :: vector< int > :: const_iterator end    = mEnrNodeIndices.end();
    int nodeInd     = iDMan.giveGlobalNumber();

    std :: vector< int > :: const_iterator it = std :: find(begin, end, nodeInd);

    if ( it != end ) {
        int enrichmentType = mNodeEnrMarker [ * it - 1 ];

        if ( enrichmentType == 1 ) {
            // Bulk enrichment
            return 1;
        } else   {
            // Front enrichment
            return mpEnrichmentFront->giveNumEnrichments(iDMan);
        }
    }

#endif

    return 0;
}

bool EnrichmentItem :: isMaterialModified(GaussPoint &iGP, Element &iEl, StructuralMaterial * &opSM) const
{
    return false;
}

void EnrichmentItem :: updateGeometry()
{
    // Update enrichments ...
    mpEnrichmentDomain->CallNodeEnrMarkerUpdate(* this, * xMan);

    // ... and create new dofs if necessary.
    createEnrichedDofs();
}

void EnrichmentItem :: propagateFronts()
{
    // Propagate interfaces
    mpPropagationLaw->propagateInterfaces(*mpEnrichmentDomain);

    updateGeometry();
}

void 
EnrichmentItem :: addEnrichmentDomain( EnrichmentDomain *ed )
{
    // Appends the enrichment domain ed to the list
    // Does not check if there is a duplicate in the list
    this->numberOfEnrichmentDomains++;
    enrichmentDomainList->growTo(numberOfEnrichmentDomains);

    ed->setNumber(this->numberOfEnrichmentDomains);
    this->enrichmentDomainList->put(this->numberOfEnrichmentDomains, ed);
}

void
EnrichmentItem :: computeDofManDofIdArray(IntArray &answer, DofManager *dMan)
{
    // Gives an array containing the dofId's that should be created as new dofs (which dofs to enrich).
	const IntArray *enrichesDofsWithIdArray = this->giveEnrichesDofsWithIdArray();

    // Number of new dofs for one enrichment function
    int eiEnrSize = enrichesDofsWithIdArray->giveSize();

    // Number of active enrichment functions
    int numEnrFunc = this->giveNumDofManEnrichments(* dMan);

    // Go through the list of dofs that the EI supports and compare with the available dofs in the dofMan.
    // Store matches in dofMask
    IntArray dofMask(eiEnrSize *numEnrFunc);
    dofMask.zero();
    int count = 0;

    for ( int i = 1; i <= numEnrFunc; i++ ) {
        for ( int j = 1; j <= eiEnrSize; j++ ) {
            if ( dMan->hasDofID( ( DofIDItem ) enrichesDofsWithIdArray->at(j) ) ) {
                count++;
                dofMask.at(count) = dMan->giveDofWithID( enrichesDofsWithIdArray->at(j) )->giveNumber();
            }
        }
    }

    answer.resize(count);
    for ( int i = 1; i <= count; i++ ) {
        answer.at(i) = this->giveStartOfDofIdPool() + i - 1;
    }

}


// Old - JB
void
EnrichmentItem :: computeDofManDofIdArray(IntArray &answer, DofManager *dMan, int enrichmentDomainNumber)
{

	// Gives an array containing the dofId's that should be created as new dofs (what dofs to enrich).
    const IntArray *enrichesDofsWithIdArray = this->giveEnrichesDofsWithIdArray();
    int eiEnrSize = enrichesDofsWithIdArray->giveSize();

    // Go through the list of dofs that the EI supports and compare with the available dofs in the dofMan.
    // If the dofMan has support for the particular dof add it to the list.
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
EnrichmentItem :: giveEIDofIdArray(IntArray &answer, int enrichmentDomainNumber) const
{
    // Returns an array containing the dof Id's of the new enrichment dofs pertinent to the ei.
    // Note: the dof managers may not support these dofs/all potential dof id's
	const IntArray *enrichesDofsWithIdArray = this->giveEnrichesDofsWithIdArray();
    int eiEnrSize = enrichesDofsWithIdArray->giveSize();

    answer.resize(eiEnrSize);
    int xDofAllocSize = eiEnrSize * this->giveNumberOfEnrDofs(); // number of new dof id's the ei will allocate
    for ( int i = 1; i <= eiEnrSize; i++ ) {
        answer.at(i) = this->giveStartOfDofIdPool() + ( enrichmentDomainNumber - 1 ) * xDofAllocSize + ( i - 1 );
    }
}


void
EnrichmentItem :: giveEIDofIdArray(IntArray &answer) const
{
    // Returns an array containing the dof Id's of the new enrichment dofs pertinent to the ei.
    // Note: the dof managers may not support these dofs/all potential dof id's
	const IntArray *enrichesDofsWithIdArray = this->giveEnrichesDofsWithIdArray();
    int eiEnrSize = enrichesDofsWithIdArray->giveSize();

    answer.resize(eiEnrSize);
    int xDofAllocSize = eiEnrSize * this->giveNumberOfEnrDofs(); // number of new dof id's the ei will allocate
    for ( int i = 1; i <= eiEnrSize; i++ ) {
        answer.at(i) = this->giveStartOfDofIdPool() + i - 1;
    }
}

void EnrichmentItem :: evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, int iNodeInd) const
{

    if(iNodeInd == -1) {
		// Bulk enrichment
		oEnrFunc.resize(1, 0.0);
		mpEnrichmentFunc->evaluateEnrFuncAt(oEnrFunc[0], iPos, iLevelSet, mpEnrichmentDomain);
    } else {
	    if(mNodeEnrMarker[iNodeInd-1] == 1)
	    {
		    // Bulk enrichment
		    oEnrFunc.resize(1, 0.0);
		    mpEnrichmentFunc->evaluateEnrFuncAt(oEnrFunc[0], iPos, iLevelSet, mpEnrichmentDomain);
	    }
	    else
	    {
		    // Front enrichment
		    mpEnrichmentFront->evaluateEnrFuncAt(oEnrFunc, iPos, iLevelSet, iNodeInd);
	    }
    }
}


//void EnrichmentItem :: evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iPos, const double &iLevelSet) const
//{
//	// Bulk enrichment
//	oEnrFunc.resize(1, 0.0);
//	mpEnrichmentFunc->evaluateEnrFuncAt(oEnrFunc[0], iPos, iLevelSet, mpEnrichmentDomain);
//}


void EnrichmentItem :: evaluateEnrFuncDerivAt(std::vector<FloatArray> &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, int iNodeInd) const
{
	if(mNodeEnrMarker[iNodeInd-1] == 1)
	{
		// Bulk enrichment
		oEnrFuncDeriv.resize(1);
		mpEnrichmentFunc->evaluateEnrFuncDerivAt(oEnrFuncDeriv[0], iPos, iLevelSet, iGradLevelSet, mpEnrichmentDomain);
	}
	else
	{
		// Front enrichment
		mpEnrichmentFront->evaluateEnrFuncDerivAt(oEnrFuncDeriv, iPos, iLevelSet, iGradLevelSet, iNodeInd);
	}
}

void 
EnrichmentItem :: updateGeometry(TimeStep *tStep, FractureManager *fMan)
{

    //this->needsUpdate = false;
    Domain *domain= this->giveDomain();   

    for ( int i = 1; i <= domain->giveNumberOfElements(); i++ ) { 
        printf( "\n -------------------------------\n");
        Element *el = domain->giveElement(i);

        for ( int j = 1; j <= fMan->failureCriterias->giveSize(); j++ ) {
            FailureCriteria *fc = fMan->failureCriterias->at(j);
            fMan->evaluateFailureCriteria(fc, el, tStep);

            if ( Delamination *dei = dynamic_cast< Delamination * > (this) )  {
                dei->updateGeometry(tStep, fMan, el, fc); //not an overloaded function, change the name
            }
        }
    }

}

void EnrichmentItem :: updateLevelSets(XfemManager &ixFemMan)
{
    int nNodes = ixFemMan.giveDomain()->giveNumberOfDofManagers();


    mLevelSetNormalDir.resize(nNodes, 0.0);
    mLevelSetTangDir.resize(nNodes, 0.0);
    mLevelSetSurfaceNormalDir.resize(nNodes, 0.0);

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


void Delamination :: updateLevelSets(XfemManager &ixFemMan)
{
    EnrichmentItem ::updateLevelSets(ixFemMan);

    int nNodes = ixFemMan.giveDomain()->giveNumberOfDofManagers();

    mLevelSetSurfaceNormalDir.resize(nNodes, 0.0);

	for ( int n = 1; n <= nNodes; n++ )
	{
		Node *node = ixFemMan.giveDomain()->giveNode(n);

		// Extract node coord
		const FloatArray &pos( *node->giveCoordinates() );

		// Calc sign dist to normal surface // New -JB
		double delta = 0.0;

        //mpEnrichmentDomain->computeSurfaceNormalSignDist(delta, pos);
        //mLevelSetSurfaceNormalDir[n-1] = delta;
        mLevelSetSurfaceNormalDir[n-1] = this->delamXiCoord;
	}

	mLevelSetsNeedUpdate = false;

}



void EnrichmentItem :: updateNodeEnrMarker(XfemManager &ixFemMan, const EnrichmentDomain_BG &iEnrichmentDomain_BG)
{
    updateLevelSets(ixFemMan);

    Domain *d = ixFemMan.giveDomain();
    int nEl = d->giveNumberOfElements();
    int nNodes = d->giveNumberOfDofManagers();

	mNodeEnrMarker.assign(nNodes, 0);
    std::vector<TipInfo> tipInfoArray;

    // Loop over elements and use the level sets to mark nodes belonging to completely cut elements.
    for ( int elIndex = 1; elIndex <= nEl; elIndex++ ) {
        Element *el = d->giveElement(elIndex);
        int nElNodes = el->giveNumberOfNodes();

        double minSignPhi  = 1, maxSignPhi         = -1;
        double minPhi = std::numeric_limits<double>::max();
        double maxPhi = std::numeric_limits<double>::min();

        FloatArray elCenter;
        elCenter.setValues(2, 0.0, 0.0);

        for ( int elNodeInd = 1; elNodeInd <= nElNodes; elNodeInd++ ) {
            int nGlob = el->giveNode(elNodeInd)->giveGlobalNumber();

            minSignPhi = std :: min( sgn(minSignPhi), sgn(mLevelSetNormalDir [ nGlob - 1 ]) );
            maxSignPhi = std :: max( sgn(maxSignPhi), sgn(mLevelSetNormalDir [ nGlob - 1 ]) );

            minPhi = std :: min( minPhi, mLevelSetNormalDir [ nGlob - 1 ] );
            maxPhi = std :: max( maxPhi, mLevelSetNormalDir [ nGlob - 1 ] );

            elCenter.at(1) += el->giveDofManager(elNodeInd)->giveCoordinate(1) / double( nElNodes );
            elCenter.at(2) += el->giveDofManager(elNodeInd)->giveCoordinate(2) / double( nElNodes );
        }


        int numEdgeIntersec = 0;

        if( minPhi*maxPhi < mLevelSetTol ) // If the level set function changes sign within the element.
        {
        	// Count the number of element edges intersected by the interface
        	int numEdges = nElNodes; // TODO: Is this assumption always true?

        	for( int edgeIndex = 1; edgeIndex <= numEdges; edgeIndex++ )
        	{
        		IntArray bNodes;
        		el->giveInterpolation()->boundaryGiveNodes(bNodes, edgeIndex);

        		int niLoc = bNodes.at( 1 );
        		int niGlob = el->giveNode(niLoc)->giveGlobalNumber();
        		int njLoc = bNodes.at( bNodes.giveSize() );
        		int njGlob = el->giveNode(njLoc)->giveGlobalNumber();

        		if( mLevelSetNormalDir[niGlob-1]*mLevelSetNormalDir[njGlob-1] < mLevelSetTol )
        		{
        			double xi = calcXiZeroLevel(mLevelSetNormalDir[niGlob-1], mLevelSetNormalDir[njGlob-1]);

        			const double &gammaS = mLevelSetTangDir[niGlob-1];
        			const double &gammaE = mLevelSetTangDir[njGlob-1];
        			double gamma = 0.5*(1.0-xi)*gammaS + 0.5*(1.0+xi)*gammaE;

        			if( gamma > 0.0 )
        			{
        				numEdgeIntersec++;
        			}
        		}
        	}


        	if(numEdgeIntersec >= 2) {
        		// If we captured a completely cut element.
        		for(int elNodeInd = 1; elNodeInd <= nElNodes; elNodeInd++)
        		{
        			int nGlob = el->giveNode(elNodeInd)->giveGlobalNumber();

        			if( mNodeEnrMarker[nGlob-1] == 0 ) {
        				mNodeEnrMarker[nGlob-1] = 1;
        			}

        		}
        	}
        	else {

        		// Store indices of elements containing an interface tip.
        		if( numEdgeIntersec == 1 )
        		{
        			TipInfo tipInfo;
        			if( mpEnrichmentDomain->giveClosestTipInfo(elCenter, tipInfo) )
        			{
        				// Prevent storage of duplicates
        				const double tol2 = 1.0e-20;
        				bool alreadyAdded = false;

        				for(size_t i = 0; i < tipInfoArray.size(); i++) {
        					if( tipInfoArray[i].mGlobalCoord.distance_square( tipInfo.mGlobalCoord ) < tol2 ) {
        						alreadyAdded = true;
        						break;
        					}
        				}

        				if(!alreadyAdded) {
        					tipInfo.mElIndex = elIndex;
        					tipInfoArray.push_back(tipInfo);
        				}
        			}
        		}
        	}
        }
	}

	// Mark tip nodes for special treatment.
	mpEnrichmentFront->MarkNodesAsFront(mNodeEnrMarker, *xMan, mLevelSetNormalDir, mLevelSetTangDir, tipInfoArray);


	// Loop over nodes and add the indices of enriched nodes.
	// Since we loop over the nodes in order from 1 to nNodes,
	// mEnrNodeIndices will automatically be sorted.
	mEnrNodeIndices.clear();
	for( int i = 1; i <= nNodes; i++)
	{
		if( mNodeEnrMarker[i-1] > 0 )
		{
			mEnrNodeIndices.push_back(i);
		}
	}

}

void EnrichmentItem :: updateNodeEnrMarker(XfemManager &ixFemMan, const DofManList &iDofManList)
{
    updateLevelSets(ixFemMan);

    Domain *d = ixFemMan.giveDomain();
    int nNodes = d->giveNumberOfDofManagers();
    mNodeEnrMarker.resize(nNodes, 0);

    // Loop over nodes in the DofManList and mark nodes as enriched.
    const std :: vector< int > &dofList = iDofManList.giveDofManList();
    for ( int i = 0; i < int( dofList.size() ); i++ ) {
        mNodeEnrMarker [ dofList [ i ] - 1 ] = 1;
        mEnrNodeIndices.push_back(dofList[i]); // new /JB
    }

    // Set level set fields to zero
    mLevelSetNormalDir.resize(nNodes, 0.0);
    mLevelSetTangDir.resize(nNodes, 0.0);
    mLevelSetSurfaceNormalDir.resize(nNodes, 0.0); // New /JB

}

void EnrichmentItem :: updateNodeEnrMarker(XfemManager &ixFemMan, const WholeDomain &iWholeDomain)
{
    // Mark all nodes for enrichment
    Domain *d = ixFemMan.giveDomain();
    int nNodes = d->giveNumberOfDofManagers();
    mNodeEnrMarker.resize(nNodes, 1);

    // Set level set fields to zero
    mLevelSetNormalDir.resize(nNodes, 0.0);
    mLevelSetTangDir.resize(nNodes, 0.0);
}

void EnrichmentItem :: createEnrichedDofs()
{
    // Creates new dofs due to the enrichment and appends them to the dof managers

    int nrDofMan = this->giveDomain()->giveNumberOfDofManagers();
    IntArray dofIdArray;

    // Create new dofs
    for ( int i = 1; i <= nrDofMan; i++ ) {
    	DofManager *dMan = this->giveDomain()->giveDofManager(i);

    	if ( isDofManEnriched(* dMan) ) {
    		computeDofManDofIdArray(dofIdArray, dMan);
    		int nDofs = dMan->giveNumberOfDofs();
    		for ( int m = 1; m <= dofIdArray.giveSize(); m++ ) {

    			if ( !dMan->hasDofID( ( DofIDItem ) ( dofIdArray.at(m) ) ) ) {
    				dMan->appendDof( new MasterDof( nDofs + m, dMan, ( DofIDItem ) ( dofIdArray.at(m) ) ) );
    			}
    		}
    	}
    }

    // TODO: Map values from old to new dofs



    // Remove old dofs
    int poolStart 	= giveStartOfDofIdPool();
    int poolEnd 	= giveEndOfDofIdPool();

    for ( int i = 1; i <= nrDofMan; i++ ) {
    	DofManager *dMan = this->giveDomain()->giveDofManager(i);

    	computeDofManDofIdArray(dofIdArray, dMan);
    	std::vector<DofIDItem> dofsToRemove;
    	int numNodeDofs = dMan->giveNumberOfDofs();
    	for(int j = 1; j <= numNodeDofs; j++) {

    		Dof *dof = dMan->giveDof(j);
    		DofIDItem dofID = dof->giveDofID();

    		if( dofID >= DofIDItem(poolStart) && dofID <= DofIDItem(poolEnd) ) {

    			bool dofIsInIdArray = false;
    			for(int k = 1; k <= dofIdArray.giveSize(); k++) {
    				if( dofID == DofIDItem(dofIdArray.at(k)) ) {
    					dofIsInIdArray = true;
    					break;
    				}
    			}

    			if(!dofIsInIdArray) {
					dofsToRemove.push_back(dofID);
    			}

    		}

    	}

    	for(size_t j = 0; j < dofsToRemove.size(); j++) {
    		dMan->removeDof(dofsToRemove[j]);
    	}

//    	if(dofsToRemove.size() > 0) {
//    		printf("Node: %d Number of dofs: %d\n", i, dMan->giveNumberOfDofs() );
//    	}

/*
    	if( dMan->giveNumberOfDofs() > 2 ) {
    		printf("dMan->giveNumberOfDofs(): %d dofs: ", dMan->giveNumberOfDofs() );
        	computeDofManDofIdArray(dofIdArray, dMan);
        	dofIdArray.printYourself();
    	}
*/
    }


}

void EnrichmentItem :: computeIntersectionPoints(std :: vector< FloatArray > &oIntersectionPoints, std :: vector< int > &oIntersectedEdgeInd, Element *element)
{
	if( isElementEnriched(element) ) {
		// Use the level set functions to compute intersection points

		// Loop over element edges; an edge is intersected if the
		// node values of the level set functions have different signs

//		int numEdges = element->giveNumberOfBoundarySides();
		int numEdges = element->giveNumberOfNodes(); // TODO: Is this assumption always true?

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


			const double &gammaS = mLevelSetTangDir[nsGlob-1];
			const double &gammaE = mLevelSetTangDir[neGlob-1];

			if( phiS*phiE < mLevelSetTol2 )
			{
				// Intersection detected

				double xi = calcXiZeroLevel(phiS, phiE);
				double gamma = 0.5*(1.0-xi)*gammaS + 0.5*(1.0+xi)*gammaE;

				// If we are inside in tangential direction
				if(gamma > 0.0)
				{
					if( fabs( phiS - phiE ) < mLevelSetTol ) {
						// If the crack is parallel to the edge.

						FloatArray ps( *( element->giveDofManager(nsLoc)->giveCoordinates() ) );
						FloatArray pe( *( element->giveDofManager(neLoc)->giveCoordinates() ) );

						// Check that the intersection points have not already been identified.
						// This may happen if the crack intersects the element exactly at a node,
						// so that intersection is detected for both element edges in that node.

						bool alreadyFound = false;

						int numPointsOld = oIntersectionPoints.size();
						for(int k = 1; k <= numPointsOld; k++)
						{
							double dist = ps.distance( oIntersectionPoints[k-1] );

							if( dist < mLevelSetTol )
							{
								alreadyFound = true;
								break;
							}

						}

						if(!alreadyFound)
						{
							oIntersectionPoints.push_back(ps);
							oIntersectedEdgeInd.push_back(edgeIndex);
						}

						alreadyFound = false;

						numPointsOld = oIntersectionPoints.size();
						for(int k = 1; k <= numPointsOld; k++)
						{
							double dist = pe.distance( oIntersectionPoints[k-1] );

							if( dist < mLevelSetTol )
							{
								alreadyFound = true;
								break;
							}

						}

						if(!alreadyFound)
						{
							oIntersectionPoints.push_back(pe);
							oIntersectedEdgeInd.push_back(edgeIndex);
						}

					}
					else {


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
							oIntersectedEdgeInd.push_back(edgeIndex);
						}
					}
				}
			}

		}


	}


}

bool EnrichmentItem :: giveElementTipCoord(FloatArray &oCoord, int iElIndex) const
{
	if( mpEnrichmentFront != NULL ) {
		return mpEnrichmentFront->giveElementTipCoord(oCoord, iElIndex);
	}
	else {
		return false;
	}
}

double EnrichmentItem :: calcXiZeroLevel(const double &iQ1, const double &iQ2) const
{
	double xi = 0.0;

	if( fabs(iQ1-iQ2) > mLevelSetTol )
	{
		xi = (iQ1+iQ2)/(iQ1-iQ2);
	}

	if( xi < -1.0)
	{
		xi = -1.0;
	}

	if( xi > 1.0)
	{
		xi = 1.0;
	}

	return xi;
}

void EnrichmentItem :: calcPolarCoord(double &oR, double &oTheta, const FloatArray &iOrigin, const FloatArray &iPos, const FloatArray &iN, const FloatArray &iT)
{
    FloatArray q;
    q.beDifferenceOf(iPos, iOrigin);
    q.normalize();

    // Compute polar coordinates
    oR = iOrigin.distance(iPos);

    if( q.dotProduct(iN) > 0.0 ) {
    	oTheta =  acos( q.dotProduct(iT) );
    }
    else {
    	oTheta = -acos( q.dotProduct(iT) );
    }
}


Inclusion :: Inclusion(int n, XfemManager *xm, Domain *aDomain) :
    EnrichmentItem(n, xm, aDomain),
    mat(NULL)
{
    mpEnrichesDofsWithIdArray->setValues(3, D_u, D_v, D_w);
}

Inclusion :: ~Inclusion()
{
    if ( mat != NULL ) {
        mat = NULL;
    }
}

bool Inclusion :: isMaterialModified(GaussPoint &iGP, Element &iEl, StructuralMaterial * &opSM) const
{
    // Check if the point is located inside the inclusion

    FloatArray N;
    FEInterpolation *interp = iEl.giveInterpolation();
    interp->evalN( N, * iGP.giveCoordinates(), FEIElementGeometryWrapper(& iEl) );

    const IntArray &elNodes = iEl.giveDofManArray();

    double levelSetGP = 0.0;
    this->interpLevelSet(levelSetGP, N, elNodes);

    if ( levelSetGP < 0.0 ) {
        opSM = static_cast< StructuralMaterial * >(mat);
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


//------------------
// DELAMINATION
//------------------

void 
Delamination :: updateGeometry(TimeStep *tStep, FractureManager *fMan, Element *el, FailureCriteria *fc)
{
    // This method needs to be updated wrt new implementation!!

    for ( int i = 1; i <= fc->quantities.size(); i++ ) {
        if ( fc->hasFailed(i) ) { // interface has failed
            printf( "Element %d fails in interface %d \n", el->giveNumber(), i );
            fMan->setUpdateFlag(true);

            // should create a new *ed at the level given by interface i iff it does not already exist
            // add coord
            FloatArray xiCoords;
            dynamic_cast <LayeredCrossSection * > ( el->giveCrossSection() )->giveInterfaceXiCoords(xiCoords);
            double xi = xiCoords.at(i); // current xi-coord
            // find if xi is in this->enrichmentDomainXiCoords
            bool flag=false;
            int num = 0;
            /*for ( int j = 1; j <= this->enrichmentDomainXiCoords.giveSize(); j++ ) {
                if ( abs(xi-this->enrichmentDomainXiCoords.at(j)) < 1.0e-6 ) {
                    flag = true;
                    num = j;
                }
            }*/


            IntArray dofManNumbers, elDofMans;
            elDofMans = el->giveDofManArray();
            for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {  
                // ugly piece of code that will skip enrichment of dofmans that have any bc's
                // which is not generally what you want
                #if 1
                bool hasBc= false;
                for ( int j = 1; j <= el->giveDofManager(i)->giveNumberOfDofs(); j++ ) {
                    if ( el->giveDofManager(i)->giveDof(j)->hasBc(tStep) ) {
                        hasBc = true;
                        continue;
                    }
                }
                #endif
                if ( !hasBc) {
                    dofManNumbers.followedBy(elDofMans.at(i));
                }
            }
            

            //dofManNumbers.printYourself();

            if ( flag ) { //in list only add dofmans
                dynamic_cast< DofManList * > ( this->giveEnrichmentDomain(num) )->addDofManagers( dofManNumbers );
                //dynamic_cast< DofManList * > ( this->giveEnrichmentDomain(num) )->updateEnrichmentDomain(dofManNumbers);
            } else { //create ed
                int numED = this->giveNumberOfEnrichmentDomains();
                EnrichmentDomain *ed = classFactory.createEnrichmentDomain( "DofManList" ); 
                DofManList *dml = dynamic_cast< DofManList * > ( ed );
                dml->addDofManagers( dofManNumbers ); // add the dofmans of the el to the list
                this->addEnrichmentDomain(ed);
                //this->enrichmentDomainXiCoords.resizeWithValues(numED+1);
                //this->enrichmentDomainXiCoords.at(numED+1) = xiCoords.at(i);

            }                        

        }
    }


}

Delamination :: Delamination(int n, XfemManager *xm, Domain *aDomain) : EnrichmentItem(n, xm, aDomain)
{
    mpEnrichesDofsWithIdArray->setValues(6, D_u, D_v, D_w, W_u, W_v, W_w);
}


IRResultType Delamination :: initializeFrom(InputRecord *ir)
{
    this->numberOfEnrichmentFunctions = 1; // must be set before EnrichmentItem :: initializeFrom(ir) is called
    EnrichmentItem :: initializeFrom(ir);
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->delamXiCoord, _IFT_Delamination_xiCoord);

    int material = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, material, _IFT_Delamination_CohesiveZoneMaterial);
    if ( material > 0 ) {
        this->mat = this->giveDomain()->giveMaterial(material);
    }
    //write an instanciate method

    //enrichmentDomainInterfaceList

    return IRRT_OK;
}


void 
Delamination :: giveActiveDelaminationXiCoords(FloatArray &xiCoords, Element *element) 
{
    // Goes through the list of delaminations and checks which are active for a given element
    int nDelam = this->giveNumberOfEnrichmentDomains(); // max possible number
    int pos = 1;
    xiCoords.resize(0);
    for ( int i = 1; i <= nDelam; i++ ) {
        if( this->isElementFullyEnrichedByEnrichmentDomain(element, i) ) {
            xiCoords.resizeWithValues(pos);
            xiCoords.at(pos) = this->giveDelaminationXiCoord(i);
            pos++;
        } 
    }
};


// remove???
double
Delamination :: giveDelaminationZCoord(int n, Element *element)
{
    //AList<double> *xiCoordList;
    int nDelam = this->giveNumberOfEnrichmentDomains(); // max possible number
    int pos = 1;
    for ( int i = 1; i <= nDelam; i++ ) {
        if ( this->isElementEnriched(element) ) {
            //xiCoordList.put( pos, this->delaminationZCoords.at(i) );
            pos++;
        }
    }

    return 0.;
}

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
    if ( xi < xi0 ) {
        return 0.0;
    } else {
        return 1.0;
    }
}


double
Delamination :: giveDelaminationGroupMidZ(int dGroup, Element *e)
{
    double zTop = 0., zBottom = 0.;
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

// unneccesary??
void
Delamination :: giveDelaminationGroupZLimits(int &dGroup, double &zTop, double &zBottom, Element *e)
{
    int nDelam = this->giveNumberOfEnrichmentDomains();
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( e->giveCrossSection() );

    if ( dGroup == 1 ) {
        zBottom = -layeredCS->giveMidSurfaceZcoordFromBottom();
        zTop    =  0.; //this->giveDelaminationZCoord(dGroup);
    } else if ( dGroup == nDelam + 1 ) {
        zBottom =  0.; //this->giveDelaminationZCoord(dGroup-1);
        zTop    = -layeredCS->giveMidSurfaceZcoordFromBottom() + layeredCS->computeIntegralThick();
    } else {
        zBottom =  0.; //this->giveDelaminationZCoord(dGroup-1);
        zTop    =  0.; //this->giveDelaminationZCoord(dGroup);
    }

#ifdef DEBUG
    if ( zBottom > zTop ) {
        OOFEM_ERROR2("giveDelaminationGroupZLimits: Bottom z-coord is larger than top z-coord in dGroup. (%i)", dGroup);
    }

#endif
}


Crack :: Crack(int n, XfemManager *xm, Domain *aDomain) : EnrichmentItem(n, xm, aDomain)
{
    mpEnrichesDofsWithIdArray->setValues(3, D_u, D_v, D_w);
}

IRResultType Crack :: initializeFrom(InputRecord *ir)
{
    this->numberOfEnrichmentFunctions = 1;
    EnrichmentItem :: initializeFrom(ir);

    return IRRT_OK;
}

REGISTER_EnrichmentFront(EnrFrontDoNothing)
REGISTER_EnrichmentFront(EnrFrontExtend)
REGISTER_EnrichmentFront(EnrFrontLinearBranchFuncRadius)

bool EnrichmentFront :: giveElementTipCoord(FloatArray &oCoord, int iElIndex) const
{
	for(size_t i = 0; i < mTipInfo.size(); i++)
	{
		if(mTipInfo[i].mElIndex == iElIndex)
		{
			oCoord = mTipInfo[i].mGlobalCoord;
			return true;
		}

	}

	return false;
}

void EnrichmentFront :: addTipIndexToNode(int iNodeInd, int iTipInd)
{
	// If the node is already enriched by the tip,
	// append the new index to the list.
	for(size_t i = 0; i < mNodeTipIndices.size(); i++) {
		if(mNodeTipIndices[i].first == iNodeInd) {
			for(size_t j = 0; j < mNodeTipIndices[i].second.size(); j++) {
				if(mNodeTipIndices[i].second[j] == iTipInd) {
					// If the index is already present, we do not
					// need to do anything
					return;
				}
			}
			mNodeTipIndices[i].second.push_back(iTipInd);
			return;
		}
	}

	// If not, create a new pair
	std::vector<int> tipIndices;
	tipIndices.push_back(iTipInd);
	std::pair<int, std::vector<int> > nodeTipInd = make_pair(iNodeInd, tipIndices);
	mNodeTipIndices.push_back(nodeTipInd);
}

void EnrichmentFront :: giveNodeTipIndices(int iNodeInd, std::vector<int> &oTipIndices) const
{
	for(size_t i = 0; i < mNodeTipIndices.size(); i++) {
		if(mNodeTipIndices[i].first == iNodeInd) {
			oTipIndices = mNodeTipIndices[i].second;
			return;
		}
	}
}

void EnrFrontExtend :: MarkNodesAsFront(std::vector<int> &ioNodeEnrMarker, XfemManager &ixFemMan, const std::vector<double> &iLevelSetNormalDir, const std::vector<double> &iLevelSetTangDir, const std::vector<TipInfo> &iTipInfo)
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
	for(int i = 0; i < int(newEnrNodes.size()); i++)
	{
		ioNodeEnrMarker[ newEnrNodes[i]-1 ] = 1;
	}
}

EnrFrontLinearBranchFuncRadius :: EnrFrontLinearBranchFuncRadius():
mEnrichmentRadius(0.0)
{
	mpBranchFunc = new LinElBranchFunction();
}

EnrFrontLinearBranchFuncRadius :: ~EnrFrontLinearBranchFuncRadius()
{
	if(mpBranchFunc != NULL)
	{
		delete mpBranchFunc;
		mpBranchFunc = NULL;
	}
}


void EnrFrontLinearBranchFuncRadius :: MarkNodesAsFront(std::vector<int> &ioNodeEnrMarker, XfemManager &ixFemMan, const std::vector<double> &iLevelSetNormalDir, const std::vector<double> &iLevelSetTangDir, const std::vector<TipInfo> &iTipInfo)
{
	// Enrich all nodes within a prescribed radius around the crack tips.
	// TODO: If performance turns out to be an issue, we may wish
	// to put the nodes in a Kd tree (or similar) to speed up searching.
	// For now, loop over all nodes.

	mTipInfo = iTipInfo;
	mNodeTipIndices.clear();

	Domain *d = ixFemMan.giveDomain();
	int nNodes = d->giveNumberOfDofManagers();

	for(int i = 1; i <= nNodes; i++)
	{
		DofManager *dMan = d->giveDofManager(i);
		const FloatArray &nodePos = *(dMan->giveCoordinates());

		for(int j = 0; j < int(iTipInfo.size()); j++)
		{
			double radius2 = iTipInfo[j].mGlobalCoord.distance_square(nodePos);

			if( radius2 < mEnrichmentRadius* mEnrichmentRadius )
			{
				ioNodeEnrMarker[i-1] = 2;
				addTipIndexToNode(i, j);
			}
		}
	}

}

int  EnrFrontLinearBranchFuncRadius :: giveNumEnrichments(const DofManager &iDMan) const
{
	std::vector<int> tipIndices;
	int nodeInd = iDMan.giveGlobalNumber();
	giveNodeTipIndices(nodeInd, tipIndices);

	return 4*tipIndices.size();
}

void EnrFrontLinearBranchFuncRadius :: evaluateEnrFuncAt(std::vector<double> &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, int iNodeInd) const
{
	oEnrFunc.clear();

	std::vector<int> tipIndices;
	giveNodeTipIndices(iNodeInd, tipIndices);

	for(size_t i = 0; i < tipIndices.size(); i++) {
		FloatArray xTip;
		int tipInd = tipIndices[i];
		xTip.setValues(2, mTipInfo[tipInd].mGlobalCoord.at(1), mTipInfo[tipInd].mGlobalCoord.at(2));

		FloatArray pos;
		pos.setValues(2, iPos.at(1), iPos.at(2));

	    // Crack tip tangent and normal
		const FloatArray &t = mTipInfo[tipInd].mTangDir;
		const FloatArray &n = mTipInfo[tipInd].mNormalDir;

		double r = 0.0, theta = 0.0;
		EnrichmentItem::calcPolarCoord(r, theta, xTip, pos, n, t);

	    mpBranchFunc->evaluateEnrFuncAt(oEnrFunc, r, theta);
	}
}

void EnrFrontLinearBranchFuncRadius :: evaluateEnrFuncDerivAt(std::vector<FloatArray> &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, int iNodeInd) const
{

	oEnrFuncDeriv.clear();

	std::vector<int> tipIndices;
	giveNodeTipIndices(iNodeInd, tipIndices);

	for(size_t i = 0; i < tipIndices.size(); i++) {

		int tipInd = tipIndices[i];
		const FloatArray &xTip = mTipInfo[tipInd].mGlobalCoord;

        // Crack tip tangent and normal
		const FloatArray &t = mTipInfo[tipInd].mTangDir;
		const FloatArray &n = mTipInfo[tipInd].mNormalDir;

		double r = 0.0, theta = 0.0;
		EnrichmentItem::calcPolarCoord(r, theta, xTip, iPos, n, t);


		size_t sizeStart = oEnrFuncDeriv.size();
        mpBranchFunc->evaluateEnrFuncDerivAt(oEnrFuncDeriv, r, theta);

        /**
         * Transform to global coordinates.
         */
        FloatMatrix E;
        E.resize(2,2);
        E.setColumn(t,1);
        E.setColumn(n,2);


        for(size_t j = sizeStart; j < oEnrFuncDeriv.size(); j++) {
        	FloatArray enrFuncDerivGlob;
        	enrFuncDerivGlob.beProductOf(E, oEnrFuncDeriv[j]);
        	oEnrFuncDeriv[j] = enrFuncDerivGlob;
        }
	}
}

IRResultType EnrFrontLinearBranchFuncRadius :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;

    IR_GIVE_FIELD(ir, mEnrichmentRadius, _IFT_EnrFrontLinearBranchFuncRadius_Radius);

    printf("In EnrFrontLinearBranchFuncRadius :: initializeFrom(): mEnrichmentRadius: %e\n", mEnrichmentRadius );

	return IRRT_OK;
}





} // end namespace oofem
