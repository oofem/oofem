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
#include "nucleationcriterion.h"
#include "sm/Elements/Shells/shell7basexfem.h"
#include "sm/EngineeringModels/structengngmodel.h"

namespace oofem {
REGISTER_XfemManager(XfemManager)

XfemManager :: XfemManager(Domain *domain)
{
    this->domain = domain;
    numberOfEnrichmentItems = -1;
    numberOfNucleationCriteria = 0;
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
	mXFEMPotentialDofIDs.clear();

    for ( auto &ei: enrichmentItemList ) {
        IntArray dofIdArray;
        ei->createEnrichedDofs();
        ei->givePotentialEIDofIdArray(dofIdArray);
//        printf("dofIdArray: "); dofIdArray.printYourself();
        mXFEMPotentialDofIDs.followedBy(dofIdArray);
    }
}

IntArray XfemManager :: giveEnrichedDofIDs(const DofManager &iDMan) const
{
    IntArray dofIdArray;

    for(int id : mXFEMPotentialDofIDs) {
    	if(iDMan.hasDofID( DofIDItem(id) )) {
    		dofIdArray.followedBy(id);
    	}
    }

    return dofIdArray;
}

void XfemManager :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_FIELD(ir, numberOfEnrichmentItems, _IFT_XfemManager_numberOfEnrichmentItems);

    IR_GIVE_OPTIONAL_FIELD(ir, numberOfNucleationCriteria, _IFT_XfemManager_numberOfNucleationCriteria);
//    printf("numberOfNucleationCriteria: %d\n", numberOfNucleationCriteria);

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
}


void XfemManager :: giveInputRecord(DynamicInputRecord &input)
{
    input.setRecordKeywordField(giveInputRecordName(), 1);

    numberOfEnrichmentItems = giveNumberOfEnrichmentItems();
    input.setField(numberOfEnrichmentItems, _IFT_XfemManager_numberOfEnrichmentItems);

    numberOfNucleationCriteria = giveNumberOfNucleationCriteria();
    input.setField(numberOfNucleationCriteria, _IFT_XfemManager_numberOfNucleationCriteria);


    input.setField(mNumGpPerTri, _IFT_XfemManager_numberOfGpPerTri);
    input.setField(mNumTriRef, _IFT_XfemManager_numberOfTriRefs);
    input.setField(mEnrDofScaleFac, _IFT_XfemManager_enrDofScaleFac);
    input.setField(doVTKExport, _IFT_XfemManager_VTKExport);
    input.setField(vtkExportFields, _IFT_XfemManager_VTKExportFields);

    if ( mDebugVTK ) {
        input.setField(1, _IFT_XfemManager_debugVTK);
    }
}

int XfemManager :: instanciateYourself(DataReader &dr)
{
    std :: string name;

    enrichmentItemList.resize(numberOfEnrichmentItems);
    for ( int i = 1; i <= numberOfEnrichmentItems; i++ ) {
        auto &mir = dr.giveInputRecord(DataReader :: IR_enrichItemRec, i);
        mir.giveRecordKeywordField(name);

        std :: unique_ptr< EnrichmentItem >ei( classFactory.createEnrichmentItem( name.c_str(), i, this, this->giveDomain() ) );
        if ( ei.get() == NULL ) {
            OOFEM_ERROR( "unknown enrichment item (%s)", name.c_str() );
        }

        ei->initializeFrom(mir);
        ei->instanciateYourself(dr);
        this->enrichmentItemList [ i - 1 ] = std :: move(ei);
    }

    mNucleationCriteria.resize(numberOfNucleationCriteria);
    for ( int i = 1; i <= numberOfNucleationCriteria; i++ ) {
        auto &mir = dr.giveInputRecord(DataReader :: IR_crackNucleationRec, i);
        mir.giveRecordKeywordField(name);

        std :: unique_ptr< NucleationCriterion >nc( classFactory.createNucleationCriterion( name.c_str(), this->giveDomain() ) );
        if ( nc.get() == NULL ) {
            OOFEM_ERROR( "Unknown nucleation criterion: (%s)", name.c_str() );
        }

        nc->initializeFrom(mir);
        nc->instanciateYourself(dr);
        this->mNucleationCriteria [ i - 1 ] = std :: move(nc);
    }


   
    return 1;
}

void XfemManager :: postInitialize()
{
    for ( int i = 1; i <= numberOfEnrichmentItems; i++ ) {
        this->giveEnrichmentItem(i)->postInitialize();
    }

     for ( int i = 1; i <= numberOfNucleationCriteria; i++ ) {
         giveNucleationCriterion(i)->postInitialize();
     }

    updateNodeEnrichmentItemMap();

}



void XfemManager :: setDomain(Domain *ipDomain)
{
    domain = ipDomain;

    for ( auto &ei : enrichmentItemList ) {
        ei->setDomain(ipDomain);
    }
}

void XfemManager :: saveContext(DataStream &stream, ContextMode mode)
{
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

        object->saveContext(stream, mode);
    }
}


void XfemManager :: restoreContext(DataStream &stream, ContextMode mode)
{
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

        obj->restoreContext(stream, mode);
    }
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

void XfemManager :: initiateFronts(bool &oAnyFronHasPropagated, TimeStep *tStep)
{
#ifdef __SM_MODULE
  oAnyFronHasPropagated = false;
        
    // Loop over EI:s and collect cross sections which have delaminaion EI:s
    IntArray CSnumbers;
    std :: vector < FloatArray > initiationFactors; initiationFactors.resize(this->domain->giveNumberOfCrossSectionModels());
    for ( auto &ei: enrichmentItemList ) {
        if ( Delamination *dei =  dynamic_cast< Delamination * >( ei.get() ) ) {
            int CSinterfaceNumber = dei->giveDelamInterfaceNum();
            for (int CSnumber : dei->giveDelamCrossSectionNum()) {
                CSnumbers.insertSortedOnce(CSnumber);
                if (initiationFactors[CSnumber-1].giveSize() < CSinterfaceNumber) {
                    initiationFactors[CSnumber-1].resizeWithValues(CSinterfaceNumber);
                }
                initiationFactors[CSnumber-1].at(CSinterfaceNumber) = dei->giveInitiationFactor();
            }
        }
    }
    
    bool failureChecked = false;
    std :: vector < IntArray > CSinterfaceNumbers; CSinterfaceNumbers.resize(CSnumbers.giveSize());
    std :: vector < IntArray > CSDofManNumbers; CSDofManNumbers.resize(CSnumbers.giveSize());
    
    for ( auto &ei: enrichmentItemList ) {

        bool eiHasPropagated = false;

        if ( Delamination *dei =  dynamic_cast< Delamination * >( ei.get() ) ) {         
            
            if ( !failureChecked ) {
                dei->findInitiationFronts(failureChecked, CSnumbers, CSinterfaceNumbers, CSDofManNumbers, initiationFactors, tStep);
            }
            
            for (int CSnum : dei->giveDelamCrossSectionNum()) {
                                
                int iCS = CSnumbers.findSorted(CSnum);
                int iInt = CSinterfaceNumbers[iCS-1].findSorted(dei->giveDelamInterfaceNum());
                if ( iInt ) {
                    // Check if nodes are viable for enrichment
                    ///TODO: this should actually not inlcude the nodes at the boundary of the delamination, since this will propagate the delamination outside.
                    IntArray delamNodes, propNodes;
                    Set *elSet = this->giveDomain()->giveSet(this->giveDomain()->giveCrossSection(CSnum)->giveSetNumber());
                    for (int elID : elSet->giveElementList() ) {
                        delamNodes.followedBy(this->giveDomain()->giveElement(elID)->giveDofManArray());
                    } 
                    delamNodes.sort();
                    delamNodes.findCommonValuesSorted(CSDofManNumbers[iCS-1], propNodes);
                    dei->initiateFronts(eiHasPropagated,propNodes);
                }
            }
        } else {
            OOFEM_ERROR(" XfemManager :: initiateFronts not implemented for other than Delamination.")
        }
        if(eiHasPropagated) {
            oAnyFronHasPropagated = true;
        }
    }
    updateNodeEnrichmentItemMap();

#else
    OOFEM_ERROR(" XfemManager :: initiateFronts not implemented for other than Delamination.")
#endif


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

bool XfemManager :: hasInitiationCriteria()
{
    for ( auto &ei: enrichmentItemList ) {
        if ( ei->hasInitiationCriteria() ) {
            return true;
        }
    }

    return false;
}

void XfemManager :: clearEnrichmentItems()
{
	enrichmentItemList.clear();
    updateNodeEnrichmentItemMap();
}

void XfemManager :: appendEnrichmentItems(std :: vector< std :: unique_ptr< EnrichmentItem > > &iEIlist)
{
	for( auto &ei : iEIlist ) {
		enrichmentItemList.push_back(std::move(ei));
	}

	numberOfEnrichmentItems = enrichmentItemList.size();
    updateNodeEnrichmentItemMap();
}

void XfemManager :: nucleateEnrichmentItems(bool &oNewItemsWereNucleated)
{
//	printf("Entering XfemManager :: nucleateEnrichmentItems\n");

	for(auto &nucCrit : mNucleationCriteria) {
		std::vector<std::unique_ptr<EnrichmentItem>> eiList = std::move(nucCrit->nucleateEnrichmentItems());

		if(eiList.size() > 0) {
//			printf("eiList.size(): %lu\n", eiList.size() );

//			if(giveNumberOfEnrichmentItems() == 0) {
//				printf("giveNumberOfEnrichmentItems() == 0\n");

				for(auto &ei : eiList) {
					enrichmentItemList.push_back(std::move(ei));
				}

				//				enrichmentItemList.push_back(std::move(ei));
				numberOfEnrichmentItems = enrichmentItemList.size();
				oNewItemsWereNucleated = true;
			    updateNodeEnrichmentItemMap();

				return;
//			}
		}
	}

	oNewItemsWereNucleated = false;
	return;
}

bool XfemManager :: hasNucleationCriteria()
{
	return ( mNucleationCriteria.size() > 0 );
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
