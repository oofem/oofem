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

#include "xfemmanager.h"
#include "inputrecord.h"
#include "intarray.h"
#include "connectivitytable.h"
#include "floatarray.h"
#include "alist.h"
#include "domain.h"
#include "enrichmentdomain.h"
#include "element.h"
#include "dofmanager.h"
#include "cltypes.h"
#include "xfemelementinterface.h"
#include "classfactory.h"
#include "masterdof.h"
#include "datareader.h"
#include "datastream.h"
#include "contextioerr.h"
#include "dynamicinputrecord.h"
#include "internalstatevaluetype.h"

namespace oofem {

XfemManager :: XfemManager(Domain *domain)
{
    this->domain = domain;
    this->enrichmentItemList = new AList< EnrichmentItem >(0);
    numberOfEnrichmentItems = -1;
    mNumGpPerTri = 12;
    doVTKExport = false;
    vtkExportFields.resize(0);
}

XfemManager :: ~XfemManager()
{
    delete enrichmentItemList;
}


InternalStateValueType 
XfemManager :: giveXFEMStateValueType(XFEMStateType type)
{
    switch ( type ) {
    case XFEMST_Enrichment:
    case XFEMST_LevelSetPhi:
    case XFEMST_LevelSetGamma:
    case XFEMST_NumIntersecPoints:
    case XFEMST_NodeEnrMarker:
        return ISVT_SCALAR;
    default:
        return ISVT_UNDEFINED;

    }
}


bool XfemManager :: isElementEnriched(const Element *elem)
{
    // Loop over all EI which asks if el is enriched.
    for ( int i = 1; i <= this->giveNumberOfEnrichmentItems(); i++ ) {
        if ( this->giveEnrichmentItem(i)->isElementEnriched(elem) ) {
            return true;
        }
    }

    return false;
}

EnrichmentItem *XfemManager :: giveEnrichmentItem(int n)
{
    // Returns the n-th enrichment item.
    if ( enrichmentItemList->includes(n) ) {
        return enrichmentItemList->at(n);
    } else {
        OOFEM_ERROR2("giveEnrichmentItem: undefined enrichmentItem (%d)", n);
    }

    return NULL;
}


void
XfemManager :: createEnrichedDofs()
{
    // Creates new dofs due to enrichment and appends them to the dof managers
    IntArray dofIdArray;

    for ( int j = 1; j <= this->giveNumberOfEnrichmentItems(); j++ ) {
        EnrichmentItem *ei = this->giveEnrichmentItem(j);
        ei->createEnrichedDofs();
    }

}


void 
XfemManager :: addEnrichedDofsTo( DofManager *dMan, IntArray &dofIdArray )
{
    int nDofs = dMan->giveNumberOfDofs();
    for ( int j = 1; j <= dofIdArray.giveSize(); j++ ) {                      
            dMan->appendDof( new MasterDof( nDofs + j, dMan, ( DofIDItem ) ( dofIdArray.at(j) ) ) );   
    }

}


IRResultType XfemManager :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, numberOfEnrichmentItems, _IFT_XfemManager_numberOfEnrichmentItems);

    IR_GIVE_OPTIONAL_FIELD(ir, mNumGpPerTri, _IFT_XfemManager_numberOfGpPerTri);

    IR_GIVE_OPTIONAL_FIELD(ir, doVTKExport, _IFT_XfemManager_VTKExport);
    if ( doVTKExport ) {
        IntArray exportFields;
        IR_GIVE_FIELD(ir, this->vtkExportFields, _IFT_XfemManager_VTKExportFields);
    }

    return IRRT_OK;
}


void XfemManager :: giveInputRecord(DynamicInputRecord &input)
{
    input.setRecordKeywordField(_IFT_XfemManager_Name, 1);
    input.setField(numberOfEnrichmentItems, _IFT_XfemManager_numberOfEnrichmentItems);
    input.setField(mNumGpPerTri, _IFT_XfemManager_numberOfGpPerTri);
    input.setField(doVTKExport, _IFT_XfemManager_VTKExport);
    input.setField(vtkExportFields, _IFT_XfemManager_VTKExportFields);

}

int XfemManager :: instanciateYourself(DataReader *dr)
{
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

        EnrichmentItem *ei = classFactory.createEnrichmentItem( name.c_str(), i, this, this->giveDomain() );
        if ( ei == NULL ) {
            OOFEM_ERROR2( "XfemManager::instanciateYourself: unknown enrichment item (%s)", name.c_str() );
        }

        ei->initializeFrom(mir);
        ei->instanciateYourself(dr);
        this->enrichmentItemList->put(i, ei);
    }

    return 1;
}

void XfemManager :: setDomain(Domain *ipDomain)
{
    domain = ipDomain;

    int numEI = enrichmentItemList->giveSize();

    for(int i = 1; i <= numEI; i++) {
        enrichmentItemList->at(i)->setDomain(ipDomain);
    }

}

contextIOResultType XfemManager :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( mode & CM_Definition ) {
        if ( !stream->write(& this->numberOfEnrichmentItems, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    for ( int i = 1; i <= this->numberOfEnrichmentItems; i++ ) {
        EnrichmentItem *obj = this->giveEnrichmentItem(i);
        if ( ( mode & CM_Definition ) ) {
            if ( !stream->write( obj->giveInputRecordName() ) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }

        if ( ( iores = obj->saveContext(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


contextIOResultType XfemManager :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( mode & CM_Definition ) {
        if ( !stream->read(& this->numberOfEnrichmentItems, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    if ( mode & CM_Definition ) {
        this->enrichmentItemList->growTo(this->numberOfEnrichmentItems);
    }

    for ( int i = 1; i <= this->numberOfEnrichmentItems; i++ ) {
        EnrichmentItem *obj;
        if ( mode & CM_Definition ) {
            std :: string name;
            if ( !stream->read(name) ) {
                THROW_CIOERR(CIO_IOERR);
            }

            obj = classFactory.createEnrichmentItem(name.c_str(), i, this, this->domain);
            enrichmentItemList->put(i, obj);
        } else {
            obj = this->giveEnrichmentItem(i);
        }

        if ( ( iores = obj->restoreContext(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    return CIO_OK;
}

void XfemManager :: updateYourself()
{
    // Update level sets
    for ( int i = 1; i <= enrichmentItemList->giveSize(); i++ ) {
        enrichmentItemList->at(i)->updateGeometry();
    }
}

void XfemManager :: propagateFronts()
{
    for ( int i = 1; i <= enrichmentItemList->giveSize(); i++ ) {
        enrichmentItemList->at(i)->propagateFronts();
    }
}

bool XfemManager :: hasPropagatingFronts()
{
    for ( int i = 1; i <= enrichmentItemList->giveSize(); i++ ) {
    	if( enrichmentItemList->at(i)->hasPropagatingFronts() ) {
    		return true;
    	}
    }

	return false;
}

} // end namespace oofem
