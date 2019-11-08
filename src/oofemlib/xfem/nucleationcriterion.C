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
    mpDomain(ipDomain)
{ }

NucleationCriterion::~NucleationCriterion() { }

std::vector<std::unique_ptr<EnrichmentItem>> NucleationCriterion::nucleateEnrichmentItems()
{
    OOFEM_ERROR("Not implemented.")

    std::vector<std::unique_ptr<EnrichmentItem>> eiList;
    return eiList;
}

void NucleationCriterion::initializeFrom(InputRecord &ir) 
{
}

int NucleationCriterion::instanciateYourself(DataReader &dr)
{
    std :: string name;

    // Instantiate enrichment function
    auto &mir = dr.giveInputRecord(DataReader :: IR_crackNucleationRec, 1);
    mir.giveRecordKeywordField(name);

    mpEnrichmentFunc = classFactory.createEnrichmentFunction( name.c_str(), 1, mpDomain );
    if ( mpEnrichmentFunc ) {
        mpEnrichmentFunc->initializeFrom(mir);
    } else {
        OOFEM_ERROR( "failed to create enrichment function (%s)", name.c_str() );
    }
    return 1;
}

void NucleationCriterion :: appendInputRecords(DynamicDataReader &oDR)
{
    auto ir = std::make_unique<DynamicInputRecord>();

    ir->setRecordKeywordField( this->giveInputRecordName(), 1 );


//    eiRec->setField(mEnrFrontIndex,                     _IFT_EnrichmentItem_front);
//    eiRec->setField(mPropLawIndex,                      _IFT_EnrichmentItem_propagationlaw);
//
//    if ( mInheritBoundaryConditions ) {
//        eiRec->setField(_IFT_EnrichmentItem_inheritbc);
//    }

    oDR.insertInputRecord(DataReader :: IR_crackNucleationRec, std::move(ir));


//    // Enrichment function
//    auto efRec = std::make_unique<DynamicInputRecord>();
//    mpEnrichmentFunc->giveInputRecord(* efRec);
//    oDR.insertInputRecord(DataReader :: IR_enrichFuncRec, std::move(efRec));
//
//    // Geometry
//    auto geoRec = std::make_unique<DynamicInputRecord>();
//    mpBasicGeometry->giveInputRecord(* geoRec);
//    oDR.insertInputRecord(DataReader :: IR_geoRec, std::move(geoRec));
//
//
//    // Enrichment front
//    if ( mEnrFrontIndex != 0 ) {
//        auto efrRecStart = std::make_unique<DynamicInputRecord>();
//        mpEnrichmentFrontStart->giveInputRecord(* efrRecStart);
//        oDR.insertInputRecord(DataReader :: IR_enrichFrontRec, std::move(efrRecStart));
//
//        auto efrRecEnd = std::make_unique<DynamicInputRecord>();
//        mpEnrichmentFrontEnd->giveInputRecord(* efrRecEnd);
//        oDR.insertInputRecord(DataReader :: IR_enrichFrontRec, std::move(efrRecEnd));
//    }
//
//    // Propagation law
//    if ( mPropLawIndex != 0 ) {
//        auto plRec = std::make_unique<DynamicInputRecord>();
//        this->mpPropagationLaw->giveInputRecord(* plRec);
//        oDR.insertInputRecord(DataReader :: IR_propagationLawRec, std::move(plRec));
//    }
}

} /* namespace oofem */
