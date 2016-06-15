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

#include "xfem/propagationlaws/plczdamageradius.h"
#include "xfem/enrichmentitems/delamination.h"
#include "Elements/Shells/shell7basexfem.h"
#include "Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"
#include <Materials/InterfaceMaterials/structuralinterfacematerial.h>
#include <Materials/structuralms.h>

#include "xfem/propagationlaw.h"
#include "xfem/tipinfo.h"
#include "xfem/enrichmentfronts/enrichmentfront.h"
#include "xfem/xfemmanager.h"
#include "classfactory.h"
#include "domain.h"
#include "mathfem.h"
#include "dynamicinputrecord.h"
#include "spatiallocalizer.h"
#include "connectivitytable.h"

namespace oofem {
REGISTER_PropagationLaw(PLCZdamageRadius)

/////////////////////////////////////////////
IRResultType PLCZdamageRadius :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    
    IR_GIVE_FIELD(ir, mIncrementRadius,         _IFT_PLCZdamageRadius_IncRadius);
    IR_GIVE_FIELD(ir, mDamageThreshold, _IFT_PLCZdamageRadius_DamageThreshold);
    IR_GIVE_OPTIONAL_FIELD(ir, this->mPropCS, _IFT_PLCZdamageRadius_PropagationCS);

    return IRRT_OK;
}

void PLCZdamageRadius :: giveInputRecord(DynamicInputRecord &input)
{
    int number = 1;
    input.setRecordKeywordField(this->giveInputRecordName(), number);

    input.setField(mIncrementRadius,            _IFT_PLCZdamageRadius_IncRadius);
    input.setField(mDamageThreshold,      _IFT_PLCZdamageRadius_DamageThreshold);
    input.setField(mPropCS,                 _IFT_PLCZdamageRadius_PropagationCS);

}

bool PLCZdamageRadius :: propagateInterface(Domain &iDomain, EnrichmentFront &iEnrFront, TipPropagation &oTipProp)
{
    if ( !iEnrFront.propagationIsAllowed() ) {
        printf("EnrichmentFront.propagationIsAllowed is false \n");
        return false;
    }
    
    const TipInfo &iTipInfo = iEnrFront.giveTipInfo();   // includes the dofman numbers which represent the boundary of the EI. 
    //tipInfo.mTipDofManNumbers.printYourself();
        
    // No listbased tip (or EI) present, so nothing to propagate. 
    if ( iTipInfo.mTipDofManNumbers.giveSize() == 0 ) {
        printf("No dofmans in tip; nothing to propagate. \n");
        return false;  
    }
    
    // Find enriched elements (only implmented for delamination and shell7xfembase)
    
    IntArray propagationDF; 
    int EIindex(iEnrFront.mEIindex);
    IntArray tipDFnumbers(iTipInfo.mTipDofManNumbers); tipDFnumbers.sort();
    EnrichmentItem *ei = iDomain.giveXfemManager()->giveEnrichmentItem(EIindex);
    
    if ( Delamination* dei = dynamic_cast < Delamination* > ( ei ) ) {
        
        // Find elements to propagate from (remove any elements not in delamination cross sections)
        IntArray EIelements, tempEIelements, CSelements;
        iDomain.giveConnectivityTable()->giveNodeNeighbourList(tempEIelements,tipDFnumbers);
        for (int CSnum : dei->giveDelamCrossSectionNum()) {
            CrossSection *CS = iDomain.giveCrossSection(CSnum);
            if (this->mPropCS.contains(CS->giveNumber()) || this->mPropCS.containsOnlyZeroes() ) {
                CSelements.followedBy(iDomain.giveSet(CS->giveSetNumber())->giveElementList()); 
            }
        }
        CSelements.sort();
        tempEIelements.findCommonValuesSorted(CSelements,EIelements);
        
        for ( int iElt : EIelements ) {
            
            bool CZdamageThresholdMet = false;
            double maxDamage(0.0);
            double CZdamage(0.0);
            Element *elt = iDomain.giveElement(iElt);
            
            if ( Shell7BaseXFEM *shellElt = dynamic_cast < Shell7BaseXFEM * > (elt) ) {
                
                int interfaceMatNumber(shellElt->giveLayeredCS()->giveInterfaceMaterialNum(dei->giveDelamInterfaceNum()));
                
                if (interfaceMatNumber) {
                    
                    StructuralInterfaceMaterial *intMat = dynamic_cast < StructuralInterfaceMaterial * > (shellElt->giveLayeredCS()->giveInterfaceMaterial(dei->giveDelamInterfaceNum()) );
                    if (intMat == 0) {
                        OOFEM_ERROR("NULL pointer to material, interface %i",dei->giveDelamInterfaceNum());
                    }
                    
                    for (GaussPoint *gp: *shellElt->czIntegrationRulesArray[ dei->giveDelamInterfaceNum() - 1 ]) {
                        
                        StructuralInterfaceMaterialStatus *intMatStatus = static_cast < StructuralInterfaceMaterialStatus * >( intMat->giveStatus(gp) );
                        if (intMatStatus == 0) {
                            OOFEM_ERROR("NULL pointer to material status");
                        }
                        CZdamage = intMatStatus->giveTempDamage();
                        if (CZdamage > maxDamage) {maxDamage = CZdamage;}
                        if (CZdamage >= this->mDamageThreshold) {
                            CZdamageThresholdMet = true;
                            break;
                        }
                    }
                    //printf(" Max damage in element %i: %f \n",shellElt->giveNumber(),maxDamage);
                } else {
                    // No interface material. Treat interface as fully damaged.
                    CZdamageThresholdMet = true;
                    CZdamage = 1.1;
                    //printf(" interface %i has no material in element %i. Interface treated as fully damaged \n",dei->giveDelamInterfaceNum(),shellElt->giveNumber());
                }
                
                if (CZdamageThresholdMet) {
                    for (int iDF : shellElt->giveDofManArray() ) {
                        //if ( tipDFnumbers.findSorted(iDF) ) {
                            propagationDF.insertSortedOnce(iDF);
                        //}
                    }
//                     if (CZdamage < 1.1) {
//                         printf(" Damage threshold (%f) met in element %i \n",CZdamage,shellElt->giveNumber());
//                     }
                }
                
            } else {
                OOFEM_ERROR("Propagation law CZ-damage not implemented for element type %s",elt->giveClassName() );
            }
        }
        
    } else {
        OOFEM_ERROR("Propagation law CZ-damage not implemented for enrichment type %s",ei->giveClassName() );
    }
    
    
    // Localise nodes within certain radius from tip nodes
    oTipProp.mPropagationDofManNumbers.clear();
    SpatialLocalizer *localizer = iDomain.giveSpatialLocalizer();
    //propagationDF.printYourself("propagationDofMans");
    
    for ( int i = 1 ; i <= propagationDF.giveSize() ; i++ ) {
        
        Node *iNode = iDomain.giveNode(propagationDF.at(i));
        const FloatArray gCoords = iNode->giveNodeCoordinates();
        
        std :: list< int > nodeList;
        localizer->giveAllNodesWithinBox(nodeList,gCoords,mIncrementRadius);
        for ( int jNode : nodeList ) {
            //printf("nodeList node %d \n",jNode);
            oTipProp.mPropagationDofManNumbers.insertSortedOnce(jNode);
        }
        
    }
    //oTipProp.mPropagationDofManNumbers.printYourself(" The following noded will be propagated to:");
    
    return true;
}
} // end namespace oofem
