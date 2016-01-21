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
#include "dof.h"
#include "../sm/CrossSections/layeredcrosssection.h"
#include "dynamicinputrecord.h"
#include "dynamicdatareader.h"
#include "xfem/enrichmentfunction.h"
#include "xfem/propagationlaw.h"
#include "xfem/xfemmanager.h"
#include "xfem/enrichmentfronts/enrichmentfrontdonothing.h"

namespace oofem {
REGISTER_EnrichmentItem(Delamination)

//------------------
// DELAMINATION
//------------------

Delamination :: Delamination(int n, XfemManager *xm, Domain *aDomain) : ListBasedEI(n, xm, aDomain)
{
    mpEnrichesDofsWithIdArray = {
        D_u, D_v, D_w, W_u, W_v, W_w
    };
    this->interfaceNum.clear();
    this->crossSectionNum = -1;
    this->matNum = 0;
    this->xiBottom = -1.0;
    this->xiTop = -1.0;
}

int Delamination :: instanciateYourself(DataReader *dr)
{
    IRResultType result; // Required by IR_GIVE_FIELD macro
    std :: string name;

    // Instantiate enrichment function
    InputRecord *mir = dr->giveInputRecord(DataReader :: IR_enrichFuncRec, 1);
    result = mir->giveRecordKeywordField(name);

    if ( result != IRRT_OK ) {
        mir->report_error(this->giveClassName(), __func__, "", result, __FILE__, __LINE__);
    }

    mpEnrichmentFunc = classFactory.createEnrichmentFunction( name.c_str(), 1, this->giveDomain() );
    if ( mpEnrichmentFunc != NULL ) {
        mpEnrichmentFunc->initializeFrom(mir);
    } else {
        OOFEM_ERROR( "failed to create enrichment function (%s)", name.c_str() );
    }


    // Instantiate enrichment domain
    mir = dr->giveInputRecord(DataReader :: IR_geoRec, 1);
    result = mir->giveRecordKeywordField(name);
    if ( result != IRRT_OK ) {
        mir->report_error(this->giveClassName(), __func__, "", result, __FILE__, __LINE__);
    }

    IntArray idList;
    IR_GIVE_FIELD(mir, idList, _IFT_ListBasedEI_list);
    for ( int i = 1; i <= idList.giveSize(); i++ ) {
        this->dofManList.push_back( idList.at(i) );
    }

    std :: sort( dofManList.begin(), this->dofManList.end() );
    //IR_GIVE_FIELD(ir, this->xi, _IFT_DofManList_DelaminationLevel);

    // Instantiate EnrichmentFront
    if ( mEnrFrontIndex == 0 ) {
        mpEnrichmentFrontStart = new EnrFrontDoNothing();
        mpEnrichmentFrontEnd = new EnrFrontDoNothing();
    } else {
        std :: string enrFrontNameStart, enrFrontNameEnd;

        InputRecord *enrFrontStartIr = dr->giveInputRecord(DataReader :: IR_enrichFrontRec, mEnrFrontIndex);
        result = enrFrontStartIr->giveRecordKeywordField(enrFrontNameStart);

        mpEnrichmentFrontStart = classFactory.createEnrichmentFront( enrFrontNameStart.c_str() );
        if ( mpEnrichmentFrontStart != NULL ) {
            mpEnrichmentFrontStart->initializeFrom(enrFrontStartIr);
        } else {
            OOFEM_ERROR( "Failed to create enrichment front (%s)", enrFrontNameStart.c_str() );
        }

        InputRecord *enrFrontEndIr = dr->giveInputRecord(DataReader :: IR_enrichFrontRec, mEnrFrontIndex);
        result = enrFrontEndIr->giveRecordKeywordField(enrFrontNameEnd);

        mpEnrichmentFrontEnd = classFactory.createEnrichmentFront( enrFrontNameEnd.c_str() );
        if ( mpEnrichmentFrontEnd != NULL ) {
            mpEnrichmentFrontEnd->initializeFrom(enrFrontEndIr);
        } else {
            OOFEM_ERROR( "Failed to create enrichment front (%s)", enrFrontNameEnd.c_str() );
        }
    }


    // Instantiate PropagationLaw
    if ( mPropLawIndex == 0 ) {
        mpPropagationLaw = new PLDoNothing();
    } else {
        std :: string propLawName;

        InputRecord *propLawir = dr->giveInputRecord(DataReader :: IR_propagationLawRec, mPropLawIndex);
        result = propLawir->giveRecordKeywordField(propLawName);

        mpPropagationLaw = classFactory.createPropagationLaw( propLawName.c_str() );
        if ( mpPropagationLaw != NULL ) {
            mpPropagationLaw->initializeFrom(propLawir);
        } else {
            OOFEM_ERROR( "Failed to create propagation law (%s)", propLawName.c_str() );
        }
    }

    // Set start of the enrichment dof pool for the given EI
    int xDofPoolAllocSize = this->giveDofPoolSize();
    this->startOfDofIdPool = this->giveDomain()->giveNextFreeDofID(xDofPoolAllocSize);
    this->endOfDofIdPool = this->startOfDofIdPool + xDofPoolAllocSize - 1;


    XfemManager *xMan = this->giveDomain()->giveXfemManager();
    //    mpEnrichmentDomain->CallNodeEnrMarkerUpdate(* this, * xMan);
    this->updateNodeEnrMarker(* xMan);


    writeVtkDebug();

    return 1;
}

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

        for ( int i = 1; i <= dofManNumbers.giveSize(); i++ ) {
            //std::list< int > :: iterator p;
            std :: vector< int > :: iterator p;
            p = std :: find( this->dofManList.begin(), this->dofManList.end(), dofManNumbers.at(i) );
            if ( p == this->dofManList.end() ) {          // if new node
                this->dofManList.push_back( dofManNumbers.at(i) );
            }
        }

        std :: sort( dofManList.begin(), this->dofManList.end() );
    }
}

void Delamination :: evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl) const
{
    if ( iLocalCoord.giveSize() != 3 ) {
        OOFEM_ERROR("iLocalCoord.giveSize() != 3")
    }

    double levelSet = iLocalCoord.at(3) - giveDelamXiCoord();

    oEnrFunc.resize(1, 0.0);
    mpEnrichmentFunc->evaluateEnrFuncAt(oEnrFunc [ 0 ], iGlobalCoord, levelSet);
}


void Delamination :: evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl, const FloatArray &iN, const IntArray &iElNodes) const
{
    evaluateEnrFuncAt(oEnrFunc, iGlobalCoord, iLocalCoord, iNodeInd, iEl);
}


IRResultType Delamination :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    result = EnrichmentItem :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    // Compute the delamination xi-coord
    IR_GIVE_FIELD(ir, this->interfaceNum, _IFT_Delamination_interfacenum); // interface number from the bottom
    IR_GIVE_FIELD(ir, this->crossSectionNum, _IFT_Delamination_csnum);

    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveDomain()->giveCrossSection(this->crossSectionNum) );
    if ( layeredCS == NULL ) {
        OOFEM_WARNING("Delamination EI requires a valid layered cross section number input: see record '%s'.", _IFT_Delamination_csnum);
        return IRRT_BAD_FORMAT;
    } else if ( this->interfaceNum.giveSize() < 1 || this->interfaceNum.giveSize() > 2 ) {
        OOFEM_WARNING("Size of record 'interfacenum' must be 1 or 2");
        return IRRT_BAD_FORMAT;
    }

    // check that interface numbers are valid
    interfaceNum.printYourself("interface num");
    for ( int i = 1; i <= this->interfaceNum.giveSize(); i++ ) {
        if ( this->interfaceNum.at(i) < 1 || this->interfaceNum.at(i) >= layeredCS->giveNumberOfLayers() ) {
            OOFEM_WARNING( "Cross section does not contain the interface number (%d) specified in the record '%s' since number of layers is %d.", this->interfaceNum.at(i), _IFT_Delamination_interfacenum, layeredCS->giveNumberOfLayers() );
            return IRRT_BAD_FORMAT;
        }
    }

    // compute xi-coord of the delamination
    this->delamXiCoord = -1.0;
    double totalThickness = layeredCS->give(CS_Thickness, FloatArray(), NULL, false); // no position available
    for ( int i = 1; i <= this->interfaceNum.at(1); i++ ) {
        this->delamXiCoord += layeredCS->giveLayerThickness(i) / totalThickness * 2.0;
        this->xiBottom += layeredCS->giveLayerThickness(i) / totalThickness * 2.0;
    }

    if ( this->interfaceNum.giveSize() == 2 ) {
        if ( this->interfaceNum.at(1) >= this->interfaceNum.at(2) ) {
            OOFEM_WARNING("second intercfacenum must be greater than the first one");
            return IRRT_BAD_FORMAT;
        }
        for ( int i = 1; i <= this->interfaceNum.at(2); i++ ) {
            this->xiTop += layeredCS->giveLayerThickness(i) / totalThickness * 2.0;
        }
    } else {
        this->xiTop = 1.0; // default is the top surface
    }


    IR_GIVE_OPTIONAL_FIELD(ir, this->matNum, _IFT_Delamination_CohesiveZoneMaterial);
    if ( this->matNum > 0 ) {
        this->mat = this->giveDomain()->giveMaterial(this->matNum);
    }

    return IRRT_OK;
}



void
Delamination :: appendInputRecords(DynamicDataReader &oDR)
{
    ///@todo almost everything is copied from EnrichmentItem :: giveInputRecord, should be written in a better way
    DynamicInputRecord *eiRec = new DynamicInputRecord();
    FEMComponent :: giveInputRecord(* eiRec);

    eiRec->setField(mEnrFrontIndex, _IFT_EnrichmentItem_front);
    eiRec->setField(mPropLawIndex, _IFT_EnrichmentItem_propagationlaw);

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

    DynamicInputRecord *geoRec = new DynamicInputRecord();
    geoRec->setRecordKeywordField(this->giveInputRecordName(), 1);

    IntArray idList;
    idList.resize( dofManList.size() );
    for ( size_t i = 0; i < dofManList.size(); i++ ) {
        idList.at(i + 1) = dofManList [ i ];
    }

    geoRec->setField(idList, _IFT_ListBasedEI_list);

    oDR.insertInputRecord(DataReader :: IR_geoRec, geoRec);

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

void Delamination :: evalLevelSetNormal(double &oLevelSet, const FloatArray &iGlobalCoord, const FloatArray &iN, const IntArray &iNodeInd) const
{
    // TODO: For consistency, this should be evaluated based on giveDelamXiCoord() /ES
    //    interpLevelSet(oLevelSet, iN, iNodeInd);
}

void Delamination :: evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iPos, const double &iLevelSet) const
{
    oEnrFunc.resize(1, 0.0);
    mpEnrichmentFunc->evaluateEnrFuncAt(oEnrFunc [ 0 ], iPos, iLevelSet);
}
} // end namespace oofem
