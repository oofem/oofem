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

#include "delamination.h"
#include "classfactory.h"
#include "fracturemanager.h"
#include "element.h"
#include "../sm/CrossSections/layeredcrosssection.h"
#include "dynamicinputrecord.h"
#include "dynamicdatareader.h"
#include "xfem/enrichmentfunction.h"
#include "xfem/enrichmentdomain.h"
#include "xfem/propagationlaw.h"

namespace oofem {
REGISTER_EnrichmentItem(Delamination)

//------------------
// DELAMINATION
//------------------

void
Delamination :: updateGeometry(FailureCriteriaStatus *fc, TimeStep *tStep)
{
    if ( fc->hasFailed( this->giveNumber() ) ) { // interface has failed
        //printf( "...fails in interface %d \n", this->giveNumber() );
        IntArray dofManNumbers, elDofMans;
        Element *el = fc->el;
        elDofMans = el->giveDofManArray();

        for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {
            // ugly piece of code that will skip enrichment of dofmans that have any bc's
            // which is not generally what you want

#if 1
            bool hasBc = false;
            for ( Dof *dof: *el->giveDofManager(i) ) {
                if ( dof->hasBc(tStep) ) {
                    hasBc = true;
                    continue;
                }
            }

#endif
            if ( !hasBc ) {
                dofManNumbers.followedBy( elDofMans.at(i) );
            }
        }

        dynamic_cast< DofManList * >( this->mpEnrichmentDomain )->addDofManagers(dofManNumbers);     // fix JB
    }
}




Delamination :: Delamination(int n, XfemManager *xm, Domain *aDomain) : EnrichmentItem(n, xm, aDomain)
{
    mpEnrichesDofsWithIdArray = {
        D_u, D_v, D_w, W_u, W_v, W_w
    };
    this->interfaceNum = -1;
    this->crossSectionNum = -1;
    this->matNum = 0;
}


IRResultType Delamination :: initializeFrom(InputRecord *ir)
{
    EnrichmentItem :: initializeFrom(ir);
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    // Compute the delamination xi-coord
    IR_GIVE_FIELD(ir, this->interfaceNum, _IFT_Delamination_interfacenum); // interface number from the bottom
    IR_GIVE_FIELD(ir, this->crossSectionNum, _IFT_Delamination_csnum);

    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveDomain()->giveCrossSection(this->crossSectionNum) );
    if ( layeredCS == NULL ) {
        OOFEM_ERROR("requires a layered cross section reference as input");
    }

    this->delamXiCoord = -1.0;
    double totalThickness = layeredCS->give(CS_Thickness, NULL, NULL, false); // no position available
    for ( int i = 1; i <= this->interfaceNum; i++ ) {
        this->delamXiCoord += layeredCS->giveLayerThickness(i) / totalThickness * 2.0;
    }


    IR_GIVE_OPTIONAL_FIELD(ir, this->matNum, _IFT_Delamination_CohesiveZoneMaterial);
    if ( this->matNum > 0 ) {
        this->mat = this->giveDomain()->giveMaterial(this->matNum);
    }

    
    this->xiBottom = -2.0;
    this->xiTop = -2.0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->xiBottom, _IFT_Delamination_xiBottom);
    IR_GIVE_OPTIONAL_FIELD(ir, this->xiTop, _IFT_Delamination_xiTop);

    return IRRT_OK;
}



void
Delamination :: appendInputRecords(DynamicDataReader &oDR)
{
    ///@todo almost everything is copied from EnrichmentItem :: giveInputRecord, should be written in a better way
    DynamicInputRecord *eiRec = new DynamicInputRecord();
    FEMComponent :: giveInputRecord(* eiRec);

    eiRec->setField(mEnrFrontIndex,                     _IFT_EnrichmentItem_front);
    eiRec->setField(mPropLawIndex,                      _IFT_EnrichmentItem_propagationlaw);

    // Delamination specific records
    eiRec->setField(this->interfaceNum, _IFT_Delamination_interfacenum);
    eiRec->setField(this->crossSectionNum, _IFT_Delamination_csnum);
    eiRec->setField(this->matNum, _IFT_Delamination_CohesiveZoneMaterial);


    oDR.insertInputRecord(DataReader :: IR_enrichItemRec, eiRec);

    // Enrichment function
    DynamicInputRecord *efRec = new DynamicInputRecord();
    mpEnrichmentFunc->giveInputRecord(* efRec);
    oDR.insertInputRecord(DataReader :: IR_enrichFuncRec, efRec);


    // Enrichment domain
    DynamicInputRecord *edRec = new DynamicInputRecord();
    mpEnrichmentDomain->giveInputRecord(* edRec);
    oDR.insertInputRecord(DataReader :: IR_geoRec, edRec);

    // Enrichment front
    if ( mEnrFrontIndex != 0 ) {
        DynamicInputRecord *efrRecStart = new DynamicInputRecord();
        mpEnrichmentFrontStart->giveInputRecord(* efrRecStart);
        oDR.insertInputRecord(DataReader :: IR_enrichFrontRec, efrRecStart);

        DynamicInputRecord *efrRecEnd = new DynamicInputRecord();
        mpEnrichmentFrontEnd->giveInputRecord(* efrRecEnd);
        oDR.insertInputRecord(DataReader :: IR_enrichFrontRec, efrRecEnd);
    }

    if ( mPropLawIndex != 0 ) {
        // Propagation law
        DynamicInputRecord *plRec = new DynamicInputRecord();
        this->mpPropagationLaw->giveInputRecord(* plRec);
        oDR.insertInputRecord(DataReader :: IR_propagationLawRec, plRec);
    }
}
} // end namespace oofem
