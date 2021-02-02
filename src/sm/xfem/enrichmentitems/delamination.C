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
#include "sm/CrossSections/layeredcrosssection.h"
#include "dynamicinputrecord.h"
#include "dynamicdatareader.h"
#include "xfem/enrichmentfunction.h"
#include "xfem/propagationlaw.h"
#include "xfem/xfemmanager.h"
#include "xfem/enrichmentfronts/enrichmentfrontdonothing.h"
#include "sm/Elements/Shells/shell7basexfem.h"
#include "spatiallocalizer.h"

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
    this->crossSectionNum.clear();
    this->matNum = 0;
    this->xiBottom = -1.0;
    this->xiTop = -1.0;
    this->initiationFactor = 1e6; // a high number
    this->initiationRadius = 0.0;
    this->recoverStresses = true; // default is to recover transverse stresses
}

int Delamination :: instanciateYourself(DataReader &dr)
{
    std :: string name;

    // Instantiate enrichment function
    {
        auto &mir = dr.giveInputRecord(DataReader :: IR_enrichFuncRec, 1);
        mir.giveRecordKeywordField(name);

        mpEnrichmentFunc = classFactory.createEnrichmentFunction( name.c_str(), 1, this->giveDomain() );
        if ( mpEnrichmentFunc ) {
            mpEnrichmentFunc->initializeFrom(mir);
        } else {
            OOFEM_ERROR( "failed to create enrichment function (%s)", name.c_str() );
        }
    }


    // Instantiate enrichment domain
    {
        auto &mir = dr.giveInputRecord(DataReader :: IR_geoRec, 1);
        mir.giveRecordKeywordField(name);

        IntArray idList;
        IR_GIVE_FIELD(mir, idList, _IFT_ListBasedEI_list);
        for ( int i = 1; i <= idList.giveSize(); i++ ) {
            this->dofManList.push_back( idList.at(i) );
        }

        std :: sort( dofManList.begin(), this->dofManList.end() );
        //IR_GIVE_FIELD(ir, this->xi, _IFT_DofManList_DelaminationLevel);
    }

    // Instantiate EnrichmentFront
    if ( mEnrFrontIndex == 0 ) {
        mpEnrichmentFrontStart = std::make_unique<EnrFrontDoNothing>(this->giveNumber());
        mpEnrichmentFrontEnd = std::make_unique<EnrFrontDoNothing>(this->giveNumber());
    } else {
        std :: string enrFrontNameStart, enrFrontNameEnd;

        auto &enrFrontStartIr = dr.giveInputRecord(DataReader :: IR_enrichFrontRec, mEnrFrontIndex);
        enrFrontStartIr.giveRecordKeywordField(enrFrontNameStart);

        mpEnrichmentFrontStart = classFactory.createEnrichmentFront( enrFrontNameStart.c_str() );
        if ( mpEnrichmentFrontStart ) {
            mpEnrichmentFrontStart->initializeFrom(enrFrontStartIr);
            //printf("EnrichmentFrontStart : %s \n", mpEnrichmentFrontStart->giveClassName()); 
        } else {
            OOFEM_ERROR( "Failed to create enrichment front (%s)", enrFrontNameStart.c_str() );
        }

        auto &enrFrontEndIr = dr.giveInputRecord(DataReader :: IR_enrichFrontRec, mEnrFrontIndex);
        enrFrontEndIr.giveRecordKeywordField(enrFrontNameEnd);

        mpEnrichmentFrontEnd = classFactory.createEnrichmentFront( enrFrontNameEnd.c_str() );
        if ( mpEnrichmentFrontEnd ) {
            mpEnrichmentFrontEnd->initializeFrom(enrFrontEndIr);
            //printf("EnrichmentFrontEnd   : %s \n", mpEnrichmentFrontEnd->giveClassName()); 
        } else {
            OOFEM_ERROR( "Failed to create enrichment front (%s)", enrFrontNameEnd.c_str() );
        }
    }


    // Instantiate PropagationLaw
    if ( mPropLawIndex == 0 ) {
        mpPropagationLaw = std::make_unique<PLDoNothing>();
    } else {
        std :: string propLawName;

        auto &propLawir = dr.giveInputRecord(DataReader :: IR_propagationLawRec, mPropLawIndex);
        propLawir.giveRecordKeywordField(propLawName);

        mpPropagationLaw = classFactory.createPropagationLaw( propLawName.c_str() );
        if ( mpPropagationLaw ) {
            mpPropagationLaw->initializeFrom(propLawir);
        } else {
            OOFEM_ERROR( "Failed to create propagation law (%s)", propLawName.c_str() );
        }
    }

    // Set start of the enrichment dof pool for the given EI
    int xDofPoolAllocSize = this->giveDofPoolSize();
    this->startOfDofIdPool = this->giveDomain()->giveNextFreeDofID(xDofPoolAllocSize);
    this->endOfDofIdPool = this->startOfDofIdPool + xDofPoolAllocSize - 1;
   
    //writeVtkDebug();

    return 1;
}

void Delamination :: postInitialize() {
    XfemManager *xMan = this->giveDomain()->giveXfemManager();
    //    mpEnrichmentDomain->CallNodeEnrMarkerUpdate(* this, * xMan);
    this->updateNodeEnrMarker(* xMan);

}


#if 0
void
Delamination :: updateGeometry(FailureCriteriaStatus *fc, TimeStep *tStep)
{
    OOFEM_ERROR( "UpdateGeometry in Delamination not used" );
    return;
    
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
#endif

bool 
Delamination :: hasInitiationCriteria()
{
    if ( this->initiationFactor < 1e6 ) { // less than default
        return true;
    }

    return false;
}


void 
Delamination :: propagateFronts(bool &oFrontsHavePropagated)
{
    oFrontsHavePropagated = false;
    
    Domain *d = this->giveDomain();
    TipPropagation tipProp;
    if ( mpPropagationLaw->propagateInterface(* giveDomain(), * mpEnrichmentFrontStart, tipProp) ) {
        // Propagate front
        
        // Check if nodes are viable for enrichment
        ///TODO: this should actually not inlcude the nodes at the boundary of the delamination, since this will propagate the delamination outside.
        IntArray delamNodes, propNodes;
        for (int CSnum : this->giveDelamCrossSectionNum()) {
            Set *elSet = d->giveSet(d->giveCrossSection(CSnum)->giveSetNumber());
            for (int elID : elSet->giveElementList() ) {
                delamNodes.followedBy(d->giveElement(elID)->giveDofManArray());
            } 
        }
        delamNodes.sort();
        
        delamNodes.findCommonValuesSorted(tipProp.mPropagationDofManNumbers, propNodes);
        //propNodes.printYourself("propNodes");
        
        bool printed = false;
        for ( int inode : propNodes ) {
            //std::list< int > :: iterator p;
            std :: vector< int > :: iterator p;
            p = std :: find( this->dofManList.begin(), this->dofManList.end(), inode );
            if ( p == this->dofManList.end() ) {          // if new node
                if ( !printed ) {
                    printf("\n Enrichment %i - The following nodes will be expanded to:",this->giveNumber());
                    printed = true;
                }
                printf(" %i", inode );
                this->dofManList.push_back( inode );
            }
        }
        printf(" \n");

        std :: sort( dofManList.begin(), this->dofManList.end() );

        oFrontsHavePropagated = true;
    }
    
    this->updateGeometry();
    
}

void
Delamination :: findInitiationFronts(bool &failureChecked, const IntArray &CSnumbers, std :: vector< IntArray > &CSinterfaceNumbers, std :: vector< IntArray > &CSDofManNumbers, std :: vector< FloatArray > &initiationFactors, TimeStep *tStep)
{           
    // Loop through all cross sections associated with delaminations 
    // Returns 
    // CSinterfaceNumbers: the failed interface number associated with the cross sections. 
    // CSDofManNumbers: the dofmanagers that should be enriched (associated with the cross sections)
    IntArray failedElementInterfaces;  
    IntArray elementNumbers; 
    SpatialLocalizer *localizer = this->giveDomain()->giveSpatialLocalizer();
    
    // NB: Assumes that elements can only be included in one cross section. 
    
    for ( int iCS = 1 ; iCS <= CSnumbers.giveSize() ; iCS++ ) {
        
        int eltSetNumber = this->giveDomain()->giveCrossSection(CSnumbers.at(iCS))->giveSetNumber();
        //printf("Cross section No. %i, set No. %i \n",CSnumbers.at(iCS),eltSetNumber);
        IntArray elementNumbers = this->giveDomain()->giveSet(eltSetNumber)->giveElementList();
        
        for ( auto eltNumber : elementNumbers ) {
            Element *elt = this->giveDomain()->giveGlobalElement(eltNumber);
            
            if ( Shell7BaseXFEM *shellElt = dynamic_cast < Shell7BaseXFEM * > (elt) ) {
                
                //bool recoverStresses = true;
                shellElt->giveFailedInterfaceNumber(failedElementInterfaces, initiationFactors[iCS-1], tStep, this->recoverStresses);
                //failedElementInterfaces.printYourself("failedElementInterfaces");
                for (int eltInt : failedElementInterfaces ) {
                    CSinterfaceNumbers[iCS-1].insertSortedOnce(eltInt);
                }
                if ( !failedElementInterfaces.isEmpty() ) {
                    for (int iDF : shellElt->giveDofManArray() ) {
                        //printf("element node %d \n",iDF);
                        if ( this->initiationRadius > 0.0 ) {
                            const auto &gCoords = this->giveDomain()->giveNode(iDF)->giveCoordinates();
                            std :: list< int > nodeList;
                            localizer->giveAllNodesWithinBox(nodeList, gCoords, initiationRadius);
                                
                            for ( int jNode : nodeList ) {
                                //printf("nodeList node %d \n",jNode);
                                CSDofManNumbers[iCS-1].insertSortedOnce(jNode);
                            }
                        } else {
                            CSDofManNumbers[iCS-1].insertSortedOnce(iDF);
                        }
                    }
                }
            }
        }
    } 
    
    failureChecked = true;
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


void Delamination :: initializeFrom(InputRecord &ir)
{
    EnrichmentItem :: initializeFrom(ir);

    // Compute the delamination xi-coord
    IR_GIVE_FIELD(ir, this->interfaceNum, _IFT_Delamination_interfacenum); // interface number from the bottom
    IR_GIVE_FIELD(ir, this->crossSectionNum, _IFT_Delamination_csnum);
    if ( ir.hasField(_IFT_Delamination_averageStresses) ) {
        this->recoverStresses = false;
        //printf("averageStresses");
    }
    
#if 1
    // update: csnum is now an IntArray of cross sections viable for delamination
    bool checkCS = false; 
    double totalThickness(0.0);
    FloatArray layerThicknesses;
    int numberOfLayers(0);
    for (int iCS : this->crossSectionNum) {
        LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveDomain()->giveCrossSection(iCS) );
        if ( !layeredCS ) {
            throw ValueInputException(ir, _IFT_Delamination_csnum, "Delamination EI requires a valid layered cross section number input");
        } else if ( this->interfaceNum.giveSize() < 1 || this->interfaceNum.giveSize() > 2 ) {
            throw ValueInputException(ir, _IFT_Delamination_interfacenum, "Size must be 1 or 2");
        }

        // check that interface numbers are valid
        //interfaceNum.printYourself("interface num");
        for ( int i = 1; i <= this->interfaceNum.giveSize(); i++ ) {
            if ( this->interfaceNum.at(i) < 1 || this->interfaceNum.at(i) >= layeredCS->giveNumberOfLayers() ) {
                throw ValueInputException(ir, _IFT_Delamination_interfacenum, "Cross section does not contain the interface number");
            }
        }
        
        
        if (checkCS) {
            if ( layeredCS->give(CS_Thickness, FloatArray(), nullptr, false) != totalThickness ) {
                throw ValueInputException(ir, _IFT_Delamination_csnum, "Delamination cross section have different totalThickness");
            }
            if ( layeredCS->giveNumberOfLayers() != numberOfLayers ) {
                throw ValueInputException(ir, _IFT_Delamination_csnum, "Delamination cross section have different number of layers");
            }
            
        } else 
        numberOfLayers = layeredCS->giveNumberOfLayers();
        totalThickness = layeredCS->give(CS_Thickness, FloatArray(), nullptr, false); // no position available
        for ( int i = 1 ; i <= numberOfLayers ; i++) {
            double layerThickness = layeredCS->giveLayerThickness(i);
            if (checkCS) {
                if ( layerThickness != layerThicknesses.at(i) ) {
                    throw ValueInputException(ir, _IFT_Delamination_csnum, "Delamination cross section have different layer thicknesses");
                }
                layerThicknesses.at(i) = layerThickness;
            } else {
                layerThicknesses.append(layeredCS->giveLayerThickness(i));
            }
        }

        // compute xi-coord of the delamination
        this->delamXiCoord = -1.0;
        this->xiBottom = -1.0;
        for ( int i = 1; i <= this->interfaceNum.at(1); i++ ) {
            this->delamXiCoord += layeredCS->giveLayerThickness(i) / totalThickness * 2.0;
            this->xiBottom += layeredCS->giveLayerThickness(i) / totalThickness * 2.0;
        }

        if ( this->interfaceNum.giveSize() == 2 ) {
            if ( this->interfaceNum.at(1) >= this->interfaceNum.at(2) ) {
                throw ValueInputException(ir, _IFT_Delamination_interfacenum, "second intercfacenum must be greater than the first one");
            }
            this->xiTop = -1.0;
            for ( int i = 1; i <= this->interfaceNum.at(2); i++ ) {
                this->xiTop += layeredCS->giveLayerThickness(i) / totalThickness * 2.0;
            }
        } else {
            this->xiTop = 1.0; // default is the top surface
        }
        checkCS = true;
    }
#else
    // old csnum (int version). NB: this was a bug since element nodes that where not part of the cross section could be enriched. 
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveDomain()->giveCrossSection(this->crossSectionNum) );
    if ( !layeredCS ) {
        throw ValueInputException(ir, _IFT_Delamination_csnum, "Delamination EI requires a valid layered cross section number input");
    } else if ( this->interfaceNum.giveSize() < 1 || this->interfaceNum.giveSize() > 2 ) {
        throw ValueInputException(ir, _IFT_Delamination_interfacenum, "size must be 1 or 2");
    }

    // check that interface numbers are valid
    //interfaceNum.printYourself("interface num");
    for ( int i = 1; i <= this->interfaceNum.giveSize(); i++ ) {
        if ( this->interfaceNum.at(i) < 1 || this->interfaceNum.at(i) >= layeredCS->giveNumberOfLayers() ) {
            throw ValueInputException(ir, _IFT_Delamination_interfacenum, "Cross section does not contain the interface number");
        }
    }

    // compute xi-coord of the delamination
    this->delamXiCoord = -1.0;
    double totalThickness = layeredCS->give(CS_Thickness, FloatArray(), nullptr, false); // no position available
    for ( int i = 1; i <= this->interfaceNum.at(1); i++ ) {
        this->delamXiCoord += layeredCS->giveLayerThickness(i) / totalThickness * 2.0;
        this->xiBottom += layeredCS->giveLayerThickness(i) / totalThickness * 2.0;
    }

    if ( this->interfaceNum.giveSize() == 2 ) {
        if ( this->interfaceNum.at(1) >= this->interfaceNum.at(2) ) {
            throw ValueInputException(ir, "second intercfacenum must be greater than the first one");
        }
        for ( int i = 1; i <= this->interfaceNum.at(2); i++ ) {
            this->xiTop += layeredCS->giveLayerThickness(i) / totalThickness * 2.0;
        }
    } else {
        this->xiTop = 1.0; // default is the top surface
    }
#endif

    IR_GIVE_OPTIONAL_FIELD(ir, this->matNum, _IFT_Delamination_CohesiveZoneMaterial);
    if ( this->matNum > 0 ) {
        this->mat = this->giveDomain()->giveMaterial(this->matNum);
    }
    
    IR_GIVE_OPTIONAL_FIELD(ir, this->initiationFactor, _IFT_Delamination_initiationFactor);
    if ( this->initiationFactor <= 0 ) {
        throw ValueInputException(ir, _IFT_Delamination_initiationFactor, "must be positive");
    }
    
    IR_GIVE_OPTIONAL_FIELD(ir, this->initiationRadius, _IFT_Delamination_initiationRadius);
    if ( this->initiationRadius < 0 ) {
        throw ValueInputException(ir, _IFT_Delamination_initiationRadius, "must be positive");
    }
}



void
Delamination :: appendInputRecords(DynamicDataReader &oDR)
{
    ///@todo almost everything is copied from EnrichmentItem :: giveInputRecord, should be written in a better way
    auto eiRec = std::make_unique<DynamicInputRecord>();
    FEMComponent :: giveInputRecord(* eiRec);

    eiRec->setField(mEnrFrontIndex, _IFT_EnrichmentItem_front);
    eiRec->setField(mPropLawIndex, _IFT_EnrichmentItem_propagationlaw);

    // Delamination specific records
    eiRec->setField(this->interfaceNum, _IFT_Delamination_interfacenum);
    eiRec->setField(this->crossSectionNum, _IFT_Delamination_csnum);
    eiRec->setField(this->matNum, _IFT_Delamination_CohesiveZoneMaterial);
    eiRec->setField(this->initiationFactor, _IFT_Delamination_initiationFactor);
    eiRec->setField(this->initiationRadius, _IFT_Delamination_initiationRadius);    
    if ( !eiRec->hasField(_IFT_Delamination_averageStresses) ) {
        recoverStresses = false;
    }


    oDR.insertInputRecord(DataReader :: IR_enrichItemRec, std::move(eiRec));

    // Enrichment function
    auto efRec = std::make_unique<DynamicInputRecord>();
    mpEnrichmentFunc->giveInputRecord(* efRec);
    oDR.insertInputRecord(DataReader :: IR_enrichFuncRec, std::move(efRec));


    // Enrichment domain

    auto geoRec = std::make_unique<DynamicInputRecord>();
    geoRec->setRecordKeywordField(this->giveInputRecordName(), 1);

    IntArray idList;
    idList.resize( dofManList.size() );
    for ( size_t i = 0; i < dofManList.size(); i++ ) {
        idList.at(i + 1) = dofManList [ i ];
    }

    geoRec->setField(idList, _IFT_ListBasedEI_list);

    oDR.insertInputRecord(DataReader :: IR_geoRec, std::move(geoRec));

    // Enrichment front
    if ( mEnrFrontIndex != 0 ) {
        auto efrRecStart = std::make_unique<DynamicInputRecord>();
        mpEnrichmentFrontStart->giveInputRecord(* efrRecStart);
        oDR.insertInputRecord(DataReader :: IR_enrichFrontRec, std::move(efrRecStart));

        auto efrRecEnd = std::make_unique<DynamicInputRecord>();
        mpEnrichmentFrontEnd->giveInputRecord(* efrRecEnd);
        oDR.insertInputRecord(DataReader :: IR_enrichFrontRec, std::move(efrRecEnd));
    }

    if ( mPropLawIndex != 0 ) {
        // Propagation law
        auto plRec = std::make_unique<DynamicInputRecord>();
        this->mpPropagationLaw->giveInputRecord(* plRec);
        oDR.insertInputRecord(DataReader :: IR_propagationLawRec, std::move(plRec));
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
