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

#include "nucleationcriterion.h"
#include "error.h"
#include "enrichmentitem.h"
#include "dynamicinputrecord.h"
#include "dynamicdatareader.h"
#include "xfem/enrichmentfunction.h"
#include "classfactory.h"

#include <memory>

namespace oofem {

NucleationCriterion::NucleationCriterion(Domain *ipDomain):
mpDomain(ipDomain),
mpEnrichmentFunc(NULL)
{

}

NucleationCriterion::~NucleationCriterion()
{
	delete mpEnrichmentFunc;
}

std::vector<std::unique_ptr<EnrichmentItem>> NucleationCriterion::nucleateEnrichmentItems() {
	OOFEM_ERROR("Not implemented.")

	std::vector<std::unique_ptr<EnrichmentItem>> eiList;
	return std::move( eiList );
}

IRResultType NucleationCriterion::initializeFrom(InputRecord *ir) {

    return IRRT_OK;
}

int NucleationCriterion::instanciateYourself(DataReader *dr) {

//	printf("Entering NucleationCriterion::instanciateYourself(DataReader *dr)\n");

    IRResultType result; // Required by IR_GIVE_FIELD macro
    std :: string name;

    // Instantiate enrichment function
    InputRecord *mir = dr->giveInputRecord(DataReader :: IR_crackNucleationRec, 1);
    result = mir->giveRecordKeywordField(name);

    if ( result != IRRT_OK ) {
        mir->report_error(this->giveClassName(), __func__, "", result, __FILE__, __LINE__);
    }

    mpEnrichmentFunc = classFactory.createEnrichmentFunction( name.c_str(), 1, mpDomain );
    if ( mpEnrichmentFunc != NULL ) {
        mpEnrichmentFunc->initializeFrom(mir);
    } else {
        OOFEM_ERROR( "failed to create enrichment function (%s)", name.c_str() );
    }

    return IRRT_OK;
}

void NucleationCriterion :: appendInputRecords(DynamicDataReader &oDR)
{
    DynamicInputRecord *ir = new DynamicInputRecord();

    ir->setRecordKeywordField( this->giveInputRecordName(), 1 );


//    eiRec->setField(mEnrFrontIndex,                     _IFT_EnrichmentItem_front);
//    eiRec->setField(mPropLawIndex,                      _IFT_EnrichmentItem_propagationlaw);
//
//    if ( mInheritBoundaryConditions ) {
//        eiRec->setField(_IFT_EnrichmentItem_inheritbc);
//    }

    oDR.insertInputRecord(DataReader :: IR_crackNucleationRec, ir);


//    // Enrichment function
//    DynamicInputRecord *efRec = new DynamicInputRecord();
//    mpEnrichmentFunc->giveInputRecord(* efRec);
//    oDR.insertInputRecord(DataReader :: IR_enrichFuncRec, efRec);
//
//    // Geometry
//    DynamicInputRecord *geoRec = new DynamicInputRecord();
//    mpBasicGeometry->giveInputRecord(* geoRec);
//    oDR.insertInputRecord(DataReader :: IR_geoRec, geoRec);
//
//
//    // Enrichment front
//    if ( mEnrFrontIndex != 0 ) {
//        DynamicInputRecord *efrRecStart = new DynamicInputRecord();
//        mpEnrichmentFrontStart->giveInputRecord(* efrRecStart);
//        oDR.insertInputRecord(DataReader :: IR_enrichFrontRec, efrRecStart);
//
//        DynamicInputRecord *efrRecEnd = new DynamicInputRecord();
//        mpEnrichmentFrontEnd->giveInputRecord(* efrRecEnd);
//        oDR.insertInputRecord(DataReader :: IR_enrichFrontRec, efrRecEnd);
//    }
//
//    // Propagation law
//    if ( mPropLawIndex != 0 ) {
//        DynamicInputRecord *plRec = new DynamicInputRecord();
//        this->mpPropagationLaw->giveInputRecord(* plRec);
//        oDR.insertInputRecord(DataReader :: IR_propagationLawRec, plRec);
//    }
}

} /* namespace oofem */
