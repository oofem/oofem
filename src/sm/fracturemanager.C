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
    ///@todo Write proper method
    
    return IRRT_OK;
}


int FractureManager :: instanciateYourself(DataReader *dr)
{
    ///@todo Write proper method
    return 1;
}








bool
FailureCriteria :: evaluateFCQuantities(Element *el, TimeStep *tStep)
{
    // Should be default implementation
    switch ( this->giveType() ) {

    default:
        return false;
    }
}






void
FractureManager :: update(TimeStep *tStep)
{


    this->setUpdateFlag(false);

    this->evaluateFailureCriterias(tStep);

    this->updateXFEM(tStep);

}



void 
FractureManager :: evaluateFailureCriterias(TimeStep *tStep) 
{
    // Go through all the failure criterias. They are responsible for their own evaluation
    this->setUpdateFlag(false);

    for ( int i = 1; i <= this->criteriaManagers.size(); i++ ) {

        printf( "\n  Evaluating failure criteria %i \n", i);
        FailureCriteriaManager *cMan = this->criteriaManagers.at(i-1);
        //if ( this->giveType() == Local ) { 

            for ( int j = 1; j <= cMan->list.size(); j++ ) {

                printf( "\n    Evaluating for element %i \n", j);

                FailureCriteria *fc = cMan->list.at(j-1);
                fc->computeFailureCriteriaQuantities(tStep);

                // temporary
                DamagedNeighborLayered *dfc = dynamic_cast<DamagedNeighborLayered *>(fc);
                dfc->evaluateFailureCriteria();
            }
        //} else if (this->giveType() == Nonlocal ) {
        //OOFEM_ERROR1("FractureManager :: evaluateFailureCriterias -Nonlocal criteria not supported yet");
        //} else {
        //OOFEM_ERROR1("FractureManager :: evaluateFailureCriterias - Unknown failure criteria");
        //}
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



bool 
FailureCriteria :: evaluateFailureCriteria() 
{
    // Should be standard implementation to check threshold values
    // Or should a criterias overload this method?

    failedFlags.resize(this->quantities.size());
    for ( int i = 1; i <= this->quantities.size(); i++ ) { // if there are several quantities like interfaces

        failedFlags.at(i-1) = false;
        for ( int j = 1; j <= this->quantities[i-1].size(); j++ ) { // all the evaluation points (often the integration points)
            FloatArray &values = this->quantities[i-1][j-1];        // quantities in each evaluation point (e.g. max stress will check stress components in different directions)
            for ( int k = 1; k <= values.giveSize(); k++ ) {
                // assumes there is only one value to compare against which is generally 
                // not true, e.g. tension/compression thresholds
                //if ( values.at(k) >= this->thresholds.at(k) ) {
                if ( values.at(k) > this->thresholds.at(1) ) {
                    failedFlags.at(i-1) = true;
                }
            }
        }
    }

    return true;
}


void 
FractureManager :: updateXFEM(FailureCriteria *fc, TimeStep *tStep)
{ 
    // Update enrichment domains based on fc's
    XfemManager *xMan = this->giveDomain()->giveXfemManager();
    EnrichmentItem *ei;
    printf( "\n Updating enrichment item geometries... \n");
    for ( int k = 1; k <= xMan->giveNumberOfEnrichmentItems(); k++ ) {    
        ei = xMan->giveEnrichmentItem(k);
        ei->updateGeometry(fc, tStep);
    }
    //if ( this->requiresUnknownsDictionaryUpdate() ) {
        xMan->createEnrichedDofs();
    //}
}


void 
FractureManager :: updateXFEM(TimeStep *tStep)
{ 
    //// Update enrichment domains based on fc's
    //XfemManager *xMan = this->giveDomain()->giveXfemManager();
    //EnrichmentItem *ei;
    //printf( "\n Updating enrichment item geometries... \n");
    //for ( int k = 1; k <= xMan->giveNumberOfEnrichmentItems(); k++ ) {    
    //    ei = xMan->giveEnrichmentItem(k);
    //    ei->updateGeometry(tStep, this); 
    //}
    ////if ( this->requiresUnknownsDictionaryUpdate() ) {
    //    xMan->createEnrichedDofs();
    ////}


    XfemManager *xMan = this->giveDomain()->giveXfemManager();
    EnrichmentItem *ei;
    
    for ( int k = 1; k <= xMan->giveNumberOfEnrichmentItems(); k++ ) {    
        printf( "\n Updating geometry of enrichment item %i ", k);
        ei = xMan->giveEnrichmentItem(k);

        for ( int i = 1; i <= this->criteriaManagers.size(); i++ ) {

            printf( "based on failure criteria %i \n", i);
            FailureCriteriaManager *cMan = this->criteriaManagers.at(i-1);

            for ( int j = 1; j <= cMan->list.size(); j++ ) { // each element

                printf( "\n  Element %i ", j);

                FailureCriteria *fc = cMan->list.at(j-1);
                ei->updateGeometry(fc, tStep);
            }
        }    
    }
    
    //if ( this->requiresUnknownsDictionaryUpdate() ) {
        xMan->createEnrichedDofs();
    //}

}


DamagedNeighborLayered :: DamagedNeighborLayered(FailureCriteriaType type, FractureManager *fMan) 
    : FailureCriteria(type,  fMan)
{
//
}


bool
DamagedNeighborLayered :: evaluateFailureCriteria() 
{
    
    // Compare quantities with thresholds    
    // go through all the layers and compare against threshold value
    failedFlags.resize(this->layerDamageValues.giveSize());
    for ( int i = 1; i <= this->failedFlags.size(); i++ ) { // if there are several quantities like interfaces

        failedFlags.at(i-1) = false;
        if ( this->layerDamageValues.at(i) > this->thresholds.at(1) ) {
            failedFlags.at(i-1) = true;
            //fMan->setUpdateFlag(true);
        }
    }
    return true;
};

} // end namespace oofem
