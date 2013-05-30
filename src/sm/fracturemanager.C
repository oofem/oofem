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

#include "fracturemanager.h"
#include "inputrecord.h"
#include "intarray.h"
#include "connectivitytable.h"
#include "floatarray.h"
#include "alist.h"
#include "domain.h"
#include "element.h"
#include "dofmanager.h"
#include "cltypes.h"
//#include "usrdefsub.h"

#include "masterdof.h"
#include "datareader.h"
#include "datastream.h"
#include "contextioerr.h"

#include "layeredcrosssection.h"
#include "xfemmanager.h"
#include "enrichmentdomain.h"
#include "classfactory.h"
namespace oofem {
FractureManager :: FractureManager(Domain *domain)
{
    this->domain = domain;
}

FractureManager :: ~FractureManager()
{
}

void
FractureManager :: clear()
{
}



IRResultType FractureManager :: initializeFrom(InputRecord *ir)
{
    /*
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, numberOfEnrichmentItems, _IFT_XfemManager_numberOfEnrichmentItems);
    */
    return IRRT_OK;
}


int FractureManager :: instanciateYourself(DataReader *dr)
{
    /*
    const char *__proc = "instanciateYourself"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro
    std :: string name;

    enrichmentItemList->growTo(numberOfEnrichmentItems);
    for ( int i = 1; i <= numberOfEnrichmentItems; i++ ) {
        InputRecord *mir = dr->giveInputRecord(DataReader :: IR_enrichItemRec, i);
        result = mir->giveRecordKeywordField(name);

        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, "", mir, result);
        }

        EnrichmentItem *ei = CreateUsrDefEnrichmentItem( name.c_str(), i, this, this->giveDomain() );
        if ( ei == NULL ) {
            OOFEM_ERROR2( "XfemManager::instanciateYourself: unknown enrichment item (%s)", name.c_str() );
        }

        ei->initializeFrom(mir);
        ei->instanciateYourself(dr);
        this->enrichmentItemList->put(i, ei);
        this->createEnrichedDofs();
    }
    */
    return 1;
}


// remove?
void 
FractureManager :: evaluateFailureCriterias(TimeStep *tStep)
{
    // For element wise evaluation of failure criteria
    /* Go through all elements and evaluate all failure criterias that are applicable to the given element.
       Should maybe be applied to element sets instead  

       How to know if a certain fc should be evaluated? what if the element has already failed? 
    */
    Domain *domain= this->giveDomain();   

    for ( int i = 1; i <= domain->giveNumberOfElements(); i++ ) { 
        printf( "\n -------------------------------\n");
        Element *el = domain->giveElement(i);
        for ( int j = 1; j <= this->failureCriterias->giveSize(); j++ ) {
            FailureCriteria *fc = this->failureCriterias->at(j);
            this->evaluateFailureCriteria(fc, el, tStep);

            if( fc->hasFailed() ) {
                printf( "Element %d fails \n", el->giveNumber() );
                this->setUpdateFlag(true);
            }
        }

    }

}


void 
FractureManager :: evaluateFailureCriteria(FailureCriteria *fc, Element *el, TimeStep *tStep) 
{
    // 1. compute quantities necessary for evaluation of fc
    // 2. evaluate fc
    // 3. if failure, update necessary geometry

    // If we have a layered cross section it needs special treatment
    // or only check CS when a special quantity is asked for, e.g interlam eval.
    // Probably must separate FC into 2 groups: regular and layered CS

    // reset fc
    fc->quantities.clear();

    // Compute fc quantities
    if ( ! fc->evaluateFCQuantities() ) { // cannot evaluate on its own, ask element for implementation through an interface
        FailureModuleElementInterface *fmInterface =
            dynamic_cast< FailureModuleElementInterface * >( el->giveInterface(FailureModuleElementInterfaceType) );
        if ( fmInterface ) { // if element supports the failure module interface
            fmInterface->computeFailureCriteriaQuantities(fc, tStep); // compute quantities
        }
    }

    // Compare quantity with thresholds
    fc->evaluateFailureCriteria();

     
}

void 
FractureManager :: updateXFEM(TimeStep *tStep)
{ 
    XfemManager *xMan = this->giveDomain()->giveXfemManager();
    EnrichmentItem *ei;
    for ( int k = 1; k <= xMan->giveNumberOfEnrichmentItems(); k++ ) {    
        ei = xMan->giveEnrichmentItem(k);
        ei->updateGeometry(tStep, this);
    }
    //if ( this->requiresUnknownsDictionaryUpdate() ) {
        xMan->createEnrichedDofs();
    //}
}


void
FractureManager :: update(TimeStep *tStep)
{
    this->setUpdateFlag(false);
    this->updateXFEM(tStep);

}


//---------------------------------------------------

bool 
FailureCriteria :: evaluateFailureCriteria() 
{
    // Compare quantity with threshold    
    failedFlags.resize(this->quantities.size());
    for ( int i = 1; i <= this->quantities.size(); i++ ) { // if there are several quantities like interfaces

        failedFlags.at(i-1) = false;
        for ( int j = 1; j <= this->quantities[i-1].size(); j++ ) { // all the evaluation points (often the integration points)
            FloatArray &values = this->quantities[i-1][j-1];        // quantities in each evaluation point (e.g. max stress will check stress components in different directions)
            for ( int k = 1; k <= values.giveSize(); k++ ) {
                // assumes there is only one value to compare against which is generally 
                // not true, e.g. tension/compression thresholds
                if ( values.at(k) >= this->thresholds.at(k) ) {
                    //this->failedFlag = true;
                    failedFlags.at(i-1) = true;
                }
            }
        }
    }

    return false;
}


} // end namespace oofem
