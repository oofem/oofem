/* $Header: /home/cvs/bp/oofem/oofemlib/src/nonlocalmaterialext.C,v 1.14.4.1 2004/04/05 15:19:43 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/


// file: nonlocalmaterialext.C

#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "element.h"
#include "timestep.h"
#include "integrationrule.h"
#include "dynalist.h"
#include "nonlocalmaterialext.h"
#include "crosssection.h"
#include "spatiallocalizer.h"
#include "domain.h"
#include "nonlocalbarrier.h"
#ifndef __MAKEDEPEND
#include <math.h>
#endif


#define NonlocalMaterialZeroWeight 1.e-3

// initialize class variable
// StateCounterType NonlocalMaterialExtensionInterface :: lastUpdatedStateCounter = 0;

void 
NonlocalMaterialExtensionInterface :: updateDomainBeforeNonlocAverage (TimeStep* atTime)
{
 int nelem, i;
 Domain *d = this->giveDomain();
 nelem = d->giveNumberOfElements();

 //if (lastUpdatedStateCounter == atTime->giveSolutionStateCounter()) return; // already updated
 if (d->giveNonlocalUpdateStateCounter() == atTime->giveSolutionStateCounter()) return; // already updated
 for (i=1; i<=nelem; i++) {
  d->giveElement(i)->updateBeforeNonlocalAverage(atTime);
 }
 // mark last update counter to prevent multiple updates
 //this->lastUpdatedStateCounter = atTime->giveSolutionStateCounter();
 d->setNonlocalUpdateStateCounter (atTime->giveSolutionStateCounter());
}

void 
NonlocalMaterialExtensionInterface :: buildNonlocalPointTable (GaussPoint* gp)
{
 //int nelem, i, j;
 int j;
 double radius, weight, elemVolume, integrationVolume = 0.;

 NonlocalMaterialStatusExtensionInterface* statusExt = 
  (NonlocalMaterialStatusExtensionInterface*) gp->giveMaterialStatus()->
   giveInterface(NonlocalMaterialStatusExtensionInterfaceType);
 dynaList<localIntegrationRecord>* iList ;

 Element* ielem;
 GaussPoint* jGp;
 IntegrationRule* iRule;

 //nelem = this->giveDomain()->giveNumberOfElements();
 
 if (!statusExt) {
   OOFEM_ERROR ("NonlocalMaterialExtensionInterface::buildNonlocalPointTable : local material status encountered");
 }

 
 // test for bounded support - if no bounded support, the nonlocal point table is 
 // big vasting of space.
 if (!this->hasBoundedSupport()) return;


 if (!statusExt->giveIntegrationDomainList()->isEmpty()) return; // already done
 iList = statusExt->giveIntegrationDomainList();

 FloatArray gpCoords, jGpCoords;
 SpatialLocalizer::elementContainerType elemSet;
 if (gp->giveElement()->computeGlobalCoordinates (gpCoords, *(gp->giveCoordinates())) == 0) {
  OOFEM_ERROR ("NonlocalMaterialExtensionInterface::buildNonlocalPointTable: computeGlobalCoordinates of target failed");
 }
 // ask for radius of influence
 this->giveSupportRadius (radius);
 // ask domain spatial localizer for list of elements with IP within this zone
 this->giveDomain()->giveSpatialLocalizer()->giveAllElementsWithIpWithinBox (elemSet, gpCoords, radius);
 // initialize iList
 

 SpatialLocalizer::elementContainerType::iterator pos;

 for (pos=elemSet.begin(); pos !=  elemSet.end(); ++pos) {
  ielem = this->giveDomain()->giveElement(*pos);
  if (regionMap.at(ielem->giveRegionNumber()) == 0) {
   iRule = ielem->giveDefaultIntegrationRulePtr ();
   for (j=0 ; j < iRule->getNumberOfIntegrationPoints() ; j++) {
    jGp = iRule->getIntegrationPoint(j) ;
    if (ielem->computeGlobalCoordinates (jGpCoords, *(jGp->giveCoordinates()))) {
     weight = this->computeWeightFunction (gpCoords, jGpCoords);
     
     /*
       if ((weight > NonlocalMaterialZeroWeight) && (!this->isBarrierActivated(gpCoords, jGpCoords))) {
     */
     this->applyBarrierConstraints (gpCoords, jGpCoords, weight);
     if (weight > NonlocalMaterialZeroWeight) {
       localIntegrationRecord ir;
       ir.nearGp = jGp;                   // store gp 
       elemVolume = weight * jGp->giveElement()->computeVolumeAround (jGp);
       ir.weight = elemVolume;            // store gp weight
       iList->pushBack(ir); // store own copy in list
       integrationVolume += elemVolume;
     }
    } else {
      OOFEM_ERROR ("NonlocalMaterialExtensionInterface::buildNonlocalPointTable: computeGlobalCoordinates of target failed");
    }
   }
  }
 } // loop over elements
 statusExt->setIntegrationScale (integrationVolume); // remember scaling factor

 /*
 // Old implementation without spatial localizer

 FloatArray jGpCoords;
 for (i=1; i<=nelem; i++) {
  ielem = this->giveDomain()->giveElement(i);
  if (regionMap.at(ielem->giveRegionNumber()) == 0) {
   iRule = ielem->giveDefaultIntegrationRulePtr ();
   for (j=0 ; j < iRule->getNumberOfIntegrationPoints() ; j++) {
    jGp = iRule->getIntegrationPoint(j) ;
    if (ielem->computeGlobalCoordinates (jGpCoords, *(jGp->giveCoordinates()))) {
      weight = this->computeWeightFunction (gpCoords, jGpCoords);
      if (weight > NonlocalMaterialZeroWeight) {
       localIntegrationRecord ir;
       ir.nearGp = jGp;                   // store gp 
       elemVolume = weight * jGp->giveElement()->computeVolumeAround (jGp);
       ir.weight = elemVolume;            // store gp weight
       iList->append(ir); // store own copy in list
       integrationVolume += elemVolume;
      }
    } else _error ("buildNonlocalPointTable: computeGlobalCoordinates failed"); 
   }
  }
 } // loop over elements
 statusExt->setIntegrationScale (integrationVolume); // remember scaling factor
 */

}

void 
NonlocalMaterialExtensionInterface :: rebuildNonlocalPointTable (GaussPoint* gp, IntArray* contributingElems)
{
  //int nelem, i, j;
  int j;
  double radius, weight, elemVolume, integrationVolume = 0.;
  
  NonlocalMaterialStatusExtensionInterface* statusExt = 
    (NonlocalMaterialStatusExtensionInterface*) gp->giveMaterialStatus()->
    giveInterface(NonlocalMaterialStatusExtensionInterfaceType);
  dynaList<localIntegrationRecord>* iList ;
  
  Element* ielem;
  GaussPoint* jGp;
  IntegrationRule* iRule;
  
  if (!statusExt) {
    OOFEM_ERROR ("NonlocalMaterialExtensionInterface::buildNonlocalPointTable : local material status encountered");
  }
  
  // test for bounded support - if no bounded support, the nonlocal point table is 
  // big vasting of space.
  if (!this->hasBoundedSupport()) return;

  iList = statusExt->giveIntegrationDomainList();
  iList ->clear();
  
  if (contributingElems == NULL) {
    // no element table provided, use standart method
    this->buildNonlocalPointTable (gp);
  } else {

    FloatArray gpCoords, jGpCoords;
    int _e, _size = contributingElems->giveSize();
    if (gp->giveElement()->computeGlobalCoordinates (gpCoords, *(gp->giveCoordinates())) == 0) {
      OOFEM_ERROR ("NonlocalMaterialExtensionInterface::buildNonlocalPointTable: computeGlobalCoordinates of target failed");
    }
    // ask for radius of influence
    this->giveSupportRadius (radius);

    // initialize iList
    for (_e=1; _e<=_size; _e++) {
      ielem = this->giveDomain()->giveElement(contributingElems->at(_e));
      if (regionMap.at(ielem->giveRegionNumber()) == 0) {
        iRule = ielem->giveDefaultIntegrationRulePtr ();
        for (j=0 ; j < iRule->getNumberOfIntegrationPoints() ; j++) {
          jGp = iRule->getIntegrationPoint(j) ;
          if (ielem->computeGlobalCoordinates (jGpCoords, *(jGp->giveCoordinates()))) {
            weight = this->computeWeightFunction (gpCoords, jGpCoords);
            
            this->applyBarrierConstraints (gpCoords, jGpCoords, weight);
            if (weight > NonlocalMaterialZeroWeight) {
              localIntegrationRecord ir;
              ir.nearGp = jGp;                   // store gp 
              elemVolume = weight * jGp->giveElement()->computeVolumeAround (jGp);
              ir.weight = elemVolume;            // store gp weight
              iList->pushBack(ir); // store own copy in list
              integrationVolume += elemVolume;
            }
          } else {
            OOFEM_ERROR ("NonlocalMaterialExtensionInterface::buildNonlocalPointTable: computeGlobalCoordinates of target failed");
          }
        }
      }
    } // loop over elements
    statusExt->setIntegrationScale (integrationVolume); // remember scaling factor
  }
}



IRResultType
NonlocalMaterialExtensionInterface :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 if (ir->hasField (IFT_NonlocalMaterialExtensionInterface_regionmap, "regionmap")) {
  IR_GIVE_FIELD (ir, regionMap, IFT_NonlocalMaterialExtensionInterface_regionmap, "regionmap"); // Macro
  if (regionMap.giveSize() != this->giveDomain()->giveNumberOfRegions()) {
   OOFEM_ERROR("NonlocalMaterialExtensionInterface::instanciateFrom: regionMap size mismatch");
  }
 } else regionMap.zero();
  
 return IRRT_OK;

}


int
NonlocalMaterialExtensionInterface::giveInputRecordString(std::string &str, bool keyword)
{
 char buff[1024];
 int i;

 if(regionMap.isEmpty() == 0){
  sprintf(buff, " regionmap %d", regionMap.giveSize());
  str += buff;
  for(i = 1; i <= regionMap.giveSize(); i++){
   sprintf(buff, " %d", regionMap.at(i));
   str += buff;
  }
 }

 return 1;
}


//int NonlocalMaterialExtension :: giveElementRegion (Element* element) {return element->giveCrossSection()->giveNumber();} 
//int NonlocalMaterialExtension :: giveNumberOfRegions () {return domain->giveNumberOfCrossSectionModels();}

/*
bool
NonlocalMaterialExtensionInterface::isBarrierActivated (const FloatArray& c1, const FloatArray& c2) const
{
 int ib, nbarrier = domain->giveNumberOfNonlocalBarriers();

 for (ib=1; ib<=nbarrier; ib++) {
  if (domain->giveNonlocalBarrier(ib)->isActivated(c1,c2)) 
   return true;
 }

 return false;
}
*/

void 
NonlocalMaterialExtensionInterface::applyBarrierConstraints (const FloatArray& gpCoords, const FloatArray& jGpCoords, double &weight) 
{
 int ib, nbarrier = domain->giveNumberOfNonlocalBarriers();
 bool shieldFlag = false;

 for (ib=1; ib<=nbarrier; ib++) {
   domain->giveNonlocalBarrier(ib)->applyConstraint(gpCoords, jGpCoords, weight, shieldFlag, this);
   if (shieldFlag) {
     weight = 0.0;
     return;
   }
 }
 
 return;
}




////////////////////////////////////////////////////////////////////////////////////////////////////

//NonlocalMaterialStatus :: NonlocalMaterialStatus(int n, Domain* d, GaussPoint* g) : 
//MaterialStatus(n,d,g),  integrationDomainList() 
NonlocalMaterialStatusExtensionInterface :: NonlocalMaterialStatusExtensionInterface() : Interface (), integrationDomainList() 
{
 integrationScale = 0.;
}

NonlocalMaterialStatusExtensionInterface :: ~NonlocalMaterialStatusExtensionInterface() 
{
 ;
}

