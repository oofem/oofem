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

#include "fracturemanager.h"
#include "inputrecord.h"
#include "intarray.h"
#include "floatarray.h"
#include "domain.h"
#include "cltypes.h"
#include "datareader.h"
#include "datastream.h"
#include "contextioerr.h"
#include "xfem/xfemmanager.h"
#include "classfactory.h"
#include "element.h"


namespace oofem {
REGISTER_FailureCriteria(DamagedNeighborLayered)


//===================================================
//  Fracture Manager
//===================================================
FractureManager :: FractureManager(Domain *domain)
{
    this->domain = domain;
    this->updateFlag = false;
}

FractureManager :: ~FractureManager() { }

void
FractureManager :: clear() { }



IRResultType FractureManager :: initializeFrom(InputRecord *ir)
{
    // Read number of failure criterias to evaluate
    IRResultType result; // Required by IR_GIVE_FIELD macro

    int numCriterias;
    IR_GIVE_FIELD(ir, numCriterias, _IFT_FracManager_numcriterias);
    this->criteriaList.resize(numCriterias);
    bool verbose = false;
    IR_GIVE_OPTIONAL_FIELD(ir, verbose, _IFT_FracManager_verbose);

#define VERBOSE


    return IRRT_OK;
}


int FractureManager :: instanciateYourself(DataReader *dr)
{
    IRResultType result; // Required by IR_GIVE_FIELD macro
    std :: string name;

    // Create and initialize all failure criterias
    for ( int i = 1; i <= ( int ) this->criteriaList.size(); i++ ) {
        InputRecord *mir = dr->giveInputRecord(DataReader :: IR_failCritRec, i);
        result = mir->giveRecordKeywordField(name);

        if ( result != IRRT_OK ) { ///@todo Make so that this can't fail.
            IR_IOERR("", mir, result);
        }

        FailureCriteria *failCriteria = classFactory.createFailureCriteria(name.c_str(), i, this);
        if ( !failCriteria ) {
            OOFEM_ERROR( "unknown failure criteria (%s)", name.c_str() );
            return 0;
        }
        failCriteria->initializeFrom(mir);

        // Case: One status per element
        if ( failCriteria->giveType() == ELLocal ) { // if ELLocal, allocate one failure criteria for each element - should be per gp
            int numEl = this->domain->giveNumberOfElements();
            failCriteria->list.resize(numEl);
            for ( int j = 1; j <= numEl; j++ ) {
                Element *el = domain->giveElement(j);
                FailureCriteriaStatus *fcs = failCriteria->CreateStatus(el, failCriteria);
                failCriteria->list.at(j - 1) = fcs;
            }
        } else if ( failCriteria->giveType() == IPLocal ) {
            OOFEM_ERROR("IPLocal criteria not supported yet");
        } else if ( failCriteria->giveType() == Nonlocal ) {
            OOFEM_ERROR("Nonlocal criteria not supported yet");
        } else {
            OOFEM_ERROR("Unknown failure criteria");
        }

        this->criteriaList.at(i - 1) = failCriteria;
    }

    return 1;
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
    // Go through all the failure criterias. These in turn keep track of the failure criteria statuses
    // which are responsible for their own evaluation

    for ( int i = 1; i <= ( int ) this->criteriaList.size(); i++ ) {
#ifdef VERBOSE
        printf("\n  Evaluating failure criteria %i \n", i);
#endif
        FailureCriteria *failCrit = this->criteriaList.at(i - 1);
        if ( failCrit->giveType() == ELLocal ) {
            for ( int j = 1; j <= ( int ) failCrit->list.size(); j++ ) {
#ifdef VERBOSE
                printf("\n    Evaluating for element %i \n", j);
#endif
                FailureCriteriaStatus *fcStatus = failCrit->list.at(j - 1);
                failCrit->computeFailureCriteriaQuantities(fcStatus, tStep);

                this->setUpdateFlag( failCrit->evaluateFailureCriteria(fcStatus) );
            }
        } else if ( failCrit->giveType() == Nonlocal ) {
            OOFEM_ERROR("Nonlocal criteria not supported yet");
        } else {
            OOFEM_ERROR("Unknown failure criteria");
        }
    }
}


void
FractureManager :: updateXFEM(TimeStep *tStep)
{
    if ( this->giveUpdateFlag() ) {
        XfemManager *xMan = this->giveDomain()->giveXfemManager();
        EnrichmentItem *ei;

        for ( int k = 1; k <= xMan->giveNumberOfEnrichmentItems(); k++ ) {
#ifdef VERBOSE
            printf("\n Updating geometry of enrichment item %i ", k);
#endif
            ei = xMan->giveEnrichmentItem(k);

            for ( int i = 1; i <= ( int ) this->criteriaList.size(); i++ ) {
#ifdef VERBOSE
                printf("based on failure criteria %i \n", i);
#endif
                FailureCriteria *failCrit = this->criteriaList.at(i - 1);

                for ( int j = 1; j <= ( int ) failCrit->list.size(); j++ ) { // each criteria (often each element)
#ifdef VERBOSE
                    printf("\n  Element %i ", j);
#endif
                    FailureCriteriaStatus *fcStatus = failCrit->list.at(j - 1);
                    ei->updateGeometry(fcStatus, tStep);
                }
            }
        }
    }
}





//===================================================
// Failure Criterias
//===================================================
#if 1



bool
FailureCriteria :: computeFailureCriteriaQuantities(FailureCriteriaStatus *fcStatus, TimeStep *tStep)
{
    Element *el = fcStatus->el;

    // If the quantity cannot be evaluated ask element for implementation through an interface
    if ( !this->evaluateFCQuantities(el, tStep) ) {
        FailureModuleElementInterface *fmInterface =
            dynamic_cast< FailureModuleElementInterface * >( el->giveInterface(FailureModuleElementInterfaceType) );

        if ( fmInterface ) { // if element supports the failure module interface
            fmInterface->computeFailureCriteriaQuantities(fcStatus, tStep); // compute quantities
        }
    }

    return true;
}


// DamagedNeighborLayered
bool
DamagedNeighborLayered :: evaluateFailureCriteria(FailureCriteriaStatus *fcStatus)
{
    // Go through all the layers and compare against threshold value
    DamagedNeighborLayeredStatus *status = dynamic_cast< DamagedNeighborLayeredStatus * >(fcStatus);
    bool criteriaFulfilled = false;
    status->failedFlags.resize( status->layerDamageValues.giveSize() );
    for ( int i = 1; i <= ( int ) status->failedFlags.size(); i++ ) { // if there are several quantities like interfaces
        status->failedFlags.at(i - 1) = false;
        if ( status->layerDamageValues.at(i) > this->DamageThreshold ) {
            status->failedFlags.at(i - 1) = true;
            criteriaFulfilled = true;
        }
    }
    return criteriaFulfilled;
};




IRResultType FailureCriteria :: initializeFrom(InputRecord *ir)
{
    //IRResultType result; // Required by IR_GIVE_FIELD macro

    return IRRT_OK;
}


IRResultType DamagedNeighborLayered :: initializeFrom(InputRecord *ir)
{
    IRResultType result; // Required by IR_GIVE_FIELD macro

    // Read damage threshold value
    IR_GIVE_FIELD(ir, this->DamageThreshold, _IFT_DamagedNeighborLayered_DamageThreshold);

    this->setType(ELLocal);

    return IRRT_OK;
}

#endif



//===================================================
// Failure Criteria Status
//===================================================
IRResultType FailureCriteriaStatus :: initializeFrom(InputRecord *ir)
{
    //IRResultType result; // Required by IR_GIVE_FIELD macro

    return IRRT_OK;
}
} // end namespace oofem
