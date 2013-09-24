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



#include "shell7basexfem.h"


//#define VERBOSE

namespace oofem {
FractureManager :: FractureManager(Domain *domain)
{
    this->domain = domain;
    this->updateFlag = false;
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
    ///@todo Write proper method
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro

    // list of the failure criterias to evaluate
    IntArray criteriaList; 
    IR_GIVE_FIELD(ir, criteriaList, _IFT_fracManager_criteriaList);

    
    this->criteriaManagers.resize( criteriaList.giveSize() );
    

    //IR_GIVE_OPTIONAL_FIELD(ir, mNumGpPerTri, _IFT_XfemManager_numberOfGpPerTri);

    return IRRT_OK;
}


int FractureManager :: instanciateYourself(DataReader *dr)
{
    ///@todo Write proper method
    
    // initialize failure criteria managers
    int numEl = this->domain->giveNumberOfElements();
    for ( int i = 1; i <= this->criteriaManagers.size(); i++ ) {
        // if local
        this->criteriaManagers.at(i-1) = new FailureCriteriaManager(Local, this);
        FailureCriteriaManager *cMan = this->criteriaManagers.at(i-1);
        cMan->list.resize(numEl);
        for ( int j = 1; j <= numEl; j++ ) { 
            cMan->list.at(j - 1) = new DamagedNeighborLayered(Local, this);
            cMan->list.at(j - 1)->thresholds.resize(1);
            cMan->list.at(j - 1)->thresholds.at(1) = -10.0;
            cMan->list.at(j - 1)->el = domain->giveElement(j);           
        }
    }

    return 1;
}


bool
FailureCriteria :: evaluateFCQuantities(Element *el, TimeStep *tStep)
{
    // Should contain calls to default implementations for evaluation of failure criteria quantities
    return false;
}






void
FractureManager :: evaluateYourself(TimeStep *tStep)
{
    
    this->setUpdateFlag(false);
    this->evaluateFailureCriterias(tStep);

}



void 
FractureManager :: evaluateFailureCriterias(TimeStep *tStep) 
{
    // Go through all the failure criteria managers. These in turn keep track of the failure criterias 
    // which are responsible for their own evaluation

    for ( int i = 1; i <= this->criteriaManagers.size(); i++ ) {

#ifdef VERBOSE
        printf( "\n  Evaluating failure criteria %i \n", i);
#endif
        FailureCriteriaManager *cMan = this->criteriaManagers.at(i-1);
        if ( cMan->giveType() == Local ) { 

            for ( int j = 1; j <= cMan->list.size(); j++ ) {
#ifdef VERBOSE
                printf( "\n    Evaluating for element %i \n", j);
#endif
                FailureCriteria *fc = cMan->list.at(j-1);
                fc->computeFailureCriteriaQuantities(tStep);

                // temporary
                //if ( fc->giveClassName() == "DamagedNeighborLayered") {
                if ( DamagedNeighborLayered *dfc = dynamic_cast<DamagedNeighborLayered *>(fc) ) {
                    this->setUpdateFlag( dfc->evaluateFailureCriteria() );
                }

            }

        } else if (cMan->giveType() == Nonlocal ) {
            OOFEM_ERROR1("FractureManager :: evaluateFailureCriterias - Nonlocal criteria not supported yet");
        } else {
            OOFEM_ERROR1("FractureManager :: evaluateFailureCriterias - Unknown failure criteria");
        }
    }
}

//---------------------------------------------------

bool 
FailureCriteria :: computeFailureCriteriaQuantities(TimeStep *tStep) 
{
    
    Element *el = this->el;

    // If the quantity cannot be evaluated ask element for implementation through an interface
    if ( ! this->evaluateFCQuantities(el, tStep) ) { 
        FailureModuleElementInterface *fmInterface =
            dynamic_cast< FailureModuleElementInterface * >( el->giveInterface(FailureModuleElementInterfaceType) );
    
        if ( fmInterface ) { // if element supports the failure module interface
            fmInterface->computeFailureCriteriaQuantities(this, tStep); // compute quantities
 
        }
    }

    return true;
}



void 
FractureManager :: updateXFEM(TimeStep *tStep)
{ 
    if ( this->giveUpdateFlag() ) {
        XfemManager *xMan = this->giveDomain()->giveXfemManager();
        EnrichmentItem *ei;
    
        for ( int k = 1; k <= xMan->giveNumberOfEnrichmentItems(); k++ ) {    
#ifdef VERBOSE
            printf( "\n Updating geometry of enrichment item %i ", k);
#endif
            ei = xMan->giveEnrichmentItem(k);

            for ( int i = 1; i <= this->criteriaManagers.size(); i++ ) {
#ifdef VERBOSE
                printf( "based on failure criteria %i \n", i);
#endif
                FailureCriteriaManager *cMan = this->criteriaManagers.at(i-1);

                for ( int j = 1; j <= cMan->list.size(); j++ ) { // each criteria (often each element)
#ifdef VERBOSE
                    printf( "\n  Element %i ", j);
#endif
                    FailureCriteria *fc = cMan->list.at(j-1);
                    ei->updateGeometry(fc, tStep);
                }
            }    
        }
    
        //xMan->createEnrichedDofs();
        //xMan->updateYourself();
    }

}




//=======================
// DamagedNeighborLayered
//=======================

bool
DamagedNeighborLayered :: evaluateFailureCriteria() 
{
    // Go through all the layers and compare against threshold value
    bool criteriaFulfilled = false;
    failedFlags.resize(this->layerDamageValues.giveSize());
    for ( int i = 1; i <= this->failedFlags.size(); i++ ) { // if there are several quantities like interfaces

        failedFlags.at(i-1) = false;
        if ( this->layerDamageValues.at(i) > this->thresholds.at(1) ) {
            this->failedFlags.at(i-1) = true;
            criteriaFulfilled = true;
        }
    }
    return criteriaFulfilled;
};

} // end namespace oofem
