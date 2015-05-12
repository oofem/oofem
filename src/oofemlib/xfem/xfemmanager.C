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
#include "domain.h"
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
#include "XFEMDebugTools.h"
#include "xfemtolerances.h"

namespace oofem {
REGISTER_XfemManager(XfemManager)

XfemManager :: XfemManager(Domain *domain)
{
    this->domain = domain;
    numberOfEnrichmentItems = -1;
    mNumGpPerTri = 12;

    // Default is no refinement of triangles.
    mNumTriRef = 0;

    // Default is no scaling of enrichment dofs.
    mEnrDofScaleFac = 1.0;

    doVTKExport = false;
    mDebugVTK = false;
    vtkExportFields.clear();

    mNodeEnrichmentItemIndices.resize(0);
    mElementEnrichmentItemIndices.clear();
    mMaterialModifyingEnrItemIndices.clear();
}

XfemManager :: ~XfemManager()
{}


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
#if 0
    // Loop over all EI which asks if el is enriched.
    for ( int i = 1; i <= this->giveNumberOfEnrichmentItems(); i++ ) {
        if ( this->giveEnrichmentItem(i)->isElementEnriched(elem) ) {
            return true;
        }
    }
#else
    // An element is enriched if one of its nodes is enriched.
    for ( int n: elem->giveDofManArray() ) {
        if ( mNodeEnrichmentItemIndices [ n - 1 ].size() > 0 ) {
            return true;
        }
    }

#endif

    return false;
}

void
XfemManager :: createEnrichedDofs()
{
    // Creates new dofs due to enrichment and appends them to the dof managers
    IntArray dofIdArray;

    for ( auto &ei: enrichmentItemList ) {
        ei->createEnrichedDofs();
    }
}

IRResultType XfemManager :: initializeFrom(InputRecord *ir)
{
    IRResultType result; // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, numberOfEnrichmentItems, _IFT_XfemManager_numberOfEnrichmentItems);

    IR_GIVE_OPTIONAL_FIELD(ir, mNumGpPerTri, _IFT_XfemManager_numberOfGpPerTri);
    IR_GIVE_OPTIONAL_FIELD(ir, mNumTriRef, _IFT_XfemManager_numberOfTriRefs);

    IR_GIVE_OPTIONAL_FIELD(ir, mEnrDofScaleFac, _IFT_XfemManager_enrDofScaleFac);
    //printf("mEnrDofScaleFac: %e\n", mEnrDofScaleFac );

    IR_GIVE_OPTIONAL_FIELD(ir, doVTKExport, _IFT_XfemManager_VTKExport);
    if ( doVTKExport ) {
        IntArray exportFields;
        IR_GIVE_FIELD(ir, this->vtkExportFields, _IFT_XfemManager_VTKExportFields);
    }

    int vtkDebug = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, vtkDebug, _IFT_XfemManager_debugVTK);
    if ( vtkDebug == 1 ) {
        mDebugVTK = true;
    }

    // TODO: Read as input.
    XfemTolerances :: setCharacteristicElementLength(0.001);

    return IRRT_OK;
}


void XfemManager :: giveInputRecord(DynamicInputRecord &input)
{
    input.setRecordKeywordField(giveInputRecordName(), 1);
    input.setField(numberOfEnrichmentItems, _IFT_XfemManager_numberOfEnrichmentItems);
    input.setField(mNumGpPerTri, _IFT_XfemManager_numberOfGpPerTri);
    input.setField(mNumTriRef, _IFT_XfemManager_numberOfTriRefs);
    input.setField(mEnrDofScaleFac, _IFT_XfemManager_enrDofScaleFac);
    input.setField(doVTKExport, _IFT_XfemManager_VTKExport);
    input.setField(vtkExportFields, _IFT_XfemManager_VTKExportFields);

    if ( mDebugVTK ) {
        input.setField(1, _IFT_XfemManager_debugVTK);
    }
}

int XfemManager :: instanciateYourself(DataReader *dr)
{
    IRResultType result; // Required by IR_GIVE_FIELD macro
    std :: string name;

    enrichmentItemList.resize(numberOfEnrichmentItems);
    for ( int i = 1; i <= numberOfEnrichmentItems; i++ ) {
        InputRecord *mir = dr->giveInputRecord(DataReader :: IR_enrichItemRec, i);
        result = mir->giveRecordKeywordField(name);

        if ( result != IRRT_OK ) {
            mir->report_error(this->giveClassName(), __func__, "", result, __FILE__, __LINE__);
        }

        std :: unique_ptr< EnrichmentItem >ei( classFactory.createEnrichmentItem( name.c_str(), i, this, this->giveDomain() ) );
        if ( ei.get() == NULL ) {
            OOFEM_ERROR( "unknown enrichment item (%s)", name.c_str() );
        }

        ei->initializeFrom(mir);
        ei->instanciateYourself(dr);
        this->enrichmentItemList [ i - 1 ] = std :: move(ei);
    }

    updateNodeEnrichmentItemMap();

    return 1;
}

void XfemManager :: setDomain(Domain *ipDomain)
{
    domain = ipDomain;

    for ( auto &ei : enrichmentItemList ) {
        ei->setDomain(ipDomain);
    }
}

contextIOResultType XfemManager :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( mode & CM_Definition ) {
        if ( !stream.write(this->numberOfEnrichmentItems) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    for ( int i = 1; i <= this->numberOfEnrichmentItems; i++ ) {
        EnrichmentItem *object = this->giveEnrichmentItem(i);
        if ( ( mode & CM_Definition ) ) {
            if ( !stream.write( object->giveInputRecordName() ) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }

        if ( ( iores = object->saveContext(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


contextIOResultType XfemManager :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( mode & CM_Definition ) {
        if ( !stream.read(this->numberOfEnrichmentItems) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        this->enrichmentItemList.resize(this->numberOfEnrichmentItems);
    }

    for ( int i = 1; i <= this->numberOfEnrichmentItems; i++ ) {
        EnrichmentItem *obj;
        if ( mode & CM_Definition ) {
            std :: string name;
            if ( !stream.read(name) ) {
                THROW_CIOERR(CIO_IOERR);
            }

            std :: unique_ptr< EnrichmentItem >ei( classFactory.createEnrichmentItem(name.c_str(), i, this, this->domain) );
            obj = ei.get();
            enrichmentItemList.insert( enrichmentItemList.begin() + i - 1, std :: move(ei) );
        } else {
            obj = this->giveEnrichmentItem(i);
        }

        if ( ( iores = obj->restoreContext(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    return CIO_OK;
}

void XfemManager :: updateYourself(TimeStep *tStep)
{
    // Update level sets
    for ( auto &ei: enrichmentItemList ) {
        ei->updateGeometry();
    }

    updateNodeEnrichmentItemMap();
}

void XfemManager :: propagateFronts(bool &oAnyFronHasPropagated)
{
    oAnyFronHasPropagated = false;

    for ( auto &ei: enrichmentItemList ) {

        bool eiHasPropagated = false;
        ei->propagateFronts(eiHasPropagated);

        if(eiHasPropagated) {
            oAnyFronHasPropagated = true;
        }
#if 0
        if ( giveVtkDebug() ) {
            GeometryBasedEI *geoEI = dynamic_cast< GeometryBasedEI * >( ei );
            if ( geoEI != NULL ) {
                std :: vector< FloatArray >points;
                geoEI->giveSubPolygon(points, -0.1, 1.1);

                std :: vector< double >x, y;
                for ( size_t j = 0; j < points.size(); j++ ) {
                    x.push_back( points [ j ].at(1) );
                    y.push_back( points [ j ].at(2) );
                }


                char fileName [ 200 ];
                sprintf( fileName, "crack%d.dat", ei->giveNumber() );
                XFEMDebugTools :: WriteArrayToGnuplot(fileName, x, y);
            }
        }
#endif
    }

    updateNodeEnrichmentItemMap();
}

bool XfemManager :: hasPropagatingFronts()
{
    for ( auto &ei: enrichmentItemList ) {
        if ( ei->hasPropagatingFronts() ) {
            return true;
        }
    }

    return false;
}

void XfemManager :: updateNodeEnrichmentItemMap()
{
    Domain *domain = giveDomain();
    int nDMan = domain->giveNumberOfDofManagers();
    mNodeEnrichmentItemIndices.clear();
    mNodeEnrichmentItemIndices.resize(nDMan);

    int nElem = domain->giveNumberOfElements();
    mElementEnrichmentItemIndices.clear();

    for ( int i = 1; i <= nElem; i++ ) {
        int elIndex = domain->giveElement(i)->giveGlobalNumber();
        int elPlaceInArray = domain->giveElementPlaceInArray(elIndex);
        if ( i != elPlaceInArray ) {
            printf("i != elPlaceInArray.\n");
            exit(0);
        }
        mElementEnrichmentItemIndices [ elPlaceInArray ].clear();
    }

    int nEI = giveNumberOfEnrichmentItems();

    for ( int eiIndex = 1; eiIndex <= nEI; eiIndex++ ) {
        EnrichmentItem *ei = giveEnrichmentItem(eiIndex);

        const std :: unordered_map< int, NodeEnrichmentType > &enrNodeInd = ei->giveEnrNodeMap();

        //for(size_t i = 0; i < enrNodeInd.size(); i++) {
        for ( auto &nodeEiPair: enrNodeInd ) {
            mNodeEnrichmentItemIndices [ nodeEiPair.first - 1 ].push_back(eiIndex);

            ConnectivityTable *ct = domain->giveConnectivityTable();
            //const IntArray *nodeElements = ct->giveDofManConnectivityArray(nodeEiPair.first);
            IntArray nodeElements;
            IntArray nodeList = {
                nodeEiPair.first
            };
            ct->giveNodeNeighbourList(nodeElements, nodeList);

            for ( int i = 1; i <= nodeElements.giveSize(); i++ ) {
                int elInd = nodeElements.at(i);

                bool found = false;
                for ( size_t j = 0; j < mElementEnrichmentItemIndices [ elInd ].size(); j++ ) {
                    if ( mElementEnrichmentItemIndices [ elInd ] [ j ] == eiIndex ) {
                        found = true;
                        break;
                    }
                }

                if ( !found ) {
                    mElementEnrichmentItemIndices [ elInd ].push_back(eiIndex);
                }
            }
        }
    }




    mMaterialModifyingEnrItemIndices.clear();
    for ( int eiIndex = 1; eiIndex <= nEI; eiIndex++ ) {
        EnrichmentItem *ei = giveEnrichmentItem(eiIndex);

        if ( ei->canModifyMaterial() ) {
            mMaterialModifyingEnrItemIndices.push_back(eiIndex);
        }
    }
}

void XfemManager :: giveElementEnrichmentItemIndices(std :: vector< int > &oElemEnrInd, int iElementIndex) const
{
    auto res = mElementEnrichmentItemIndices.find(iElementIndex);
    if ( res != mElementEnrichmentItemIndices.end() ) {
        oElemEnrInd = res->second;
    }
}
} // end namespace oofem
