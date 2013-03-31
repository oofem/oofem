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
#include "conTable.h"
#include "flotarry.h"
#include "alist.h"
#include "domain.h"
#include "element.h"
#include "dofmanager.h"
#include "cltypes.h"
#include "usrdefsub.h"
#include "masterdof.h"
#include "datareader.h"
#include "datastream.h"
#include "contextioerr.h"

#include "layeredcrosssection.h"
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



void 
FractureManager :: evaluateFailureCriterias(TimeStep *tStep)
{
    // For element wise evaluation of failure criteria
    Domain *domain;
    domain = this->giveDomain();   
    Element *el;
    FailureCriteria *fc;
    FailureCriteriaType type = FC_MaxShearStress;
    for ( int i = 1; i <= domain->giveNumberOfElements(); i++ ) {
        el = domain->giveElement(i);
        for ( int j = 1; j <= this->failureCriterias->giveSize(); j++ ) {
            fc = this->failureCriterias->at(j);
            this->evaluateFailureCriteria(fc, type, el, tStep);
        }

    }
     
    //EnrichmentItem *ei;
    //for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        //domain = this->giveDomain(idomain);        
        //xMan = domain->giveXfemManager(1);
        

            //for ( int j = 1; j <= xMan->giveNumberOfEnrichmentItems(); j++ ) {    
              //  ei = xMan->giveEnrichmentItem(j);
                // Different actions depending on ei
                // Delamination 
               // if ( Delamination *dei = dynamic_cast< Delamination * > (ei) )  {
               //     for ( int k = 1; k <= ei->giveNumberOfEnrichmentDomains(); k++ ) {                
               //         EnrichmentDomain *ed; 
               //         ed = ei->giveEnrichmentDomain(k);
               //         this->evaluatePropagationLawForDelamination(el, ed, tStep); 
               //     }
               // }
           // }
       // }

    //}

    // Create additional dofs if any Enrichment Domain is updated
    /*
    if ( this->requiresUnknownsDictionaryUpdate() ) {
        xMan->createEnrichedDofs();
    }
    */
}


void 
FractureManager :: evaluateFailureCriteria(FailureCriteria *fc, FailureCriteriaType type, Element *el, TimeStep *tStep) 
{
    // If we have a layered cross section it needs special treatment
    if ( dynamic_cast< LayeredCrossSection *>( el->giveCrossSection() ) ) {
        FailureModuleElementInterface *fmInterface =
            dynamic_cast< FailureModuleElementInterface * >( el->giveInterface(FailureModuleElementInterfaceType) );
        if ( fmInterface ) { // if element supports the failure interface
            std::vector < FloatArray > quantities;
            fmInterface->evaluateFailureCriteriaQuantities(fc, tStep); // computeQuantities
        }
    }


    // Compare quantity with threshold
    printf( "\n -------------------------------\n");
    fc->evaluateFailureCriteria();
    /*
    for ( int i = 1; i <= fc->quantities.size(); i++ ) {
        for ( int j = 1; j <= fc->quantities[i-1].size(); j++ ) {
            FloatArray &values = fc->quantities[i-1][j-1];
            for ( int k = 1; k <= values.giveSize(); k++ ) {
                // assumeds that there is only one value to compare against which is generally 
                // not true, e.g. tension/compression thresholds
                if ( values.at(k) >= fc->thresholds.at(k) ) {
                    printf( "Eval. point %d in interface %d in element %d fails (%e >= %e) \n", 
                        j, i, el->giveNumber(), values.at(k), fc->thresholds.at(k) );
                }
            }
        }
    }
    */

    // Store which layers, interfaces that has failed and send it back 
}


bool 
FailureCriteria :: evaluateFailureCriteria() 
{
    // Compare quantity with threshold
    // should save bool for all values
    printf( "\n -------------------------------\n");
    for ( int i = 1; i <= this->quantities.size(); i++ ) {
        for ( int j = 1; j <= this->quantities[i-1].size(); j++ ) {
            FloatArray &values = this->quantities[i-1][j-1];
            for ( int k = 1; k <= values.giveSize(); k++ ) {
                // assumeds that there is only one value to compare against which is generally 
                // not true, e.g. tension/compression thresholds
                if ( values.at(k) >= this->thresholds.at(k) ) {
                    printf( "Eval. point %d in interface %d fails (%e >= %e) \n", 
                        j, i, values.at(k), this->thresholds.at(k) );
                    return true;
                }
            }
        }
    }
    return false;
}


} // end namespace oofem
