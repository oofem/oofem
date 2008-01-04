/* $Header: /home/cvs/bp/oofem/sm/src/directerrorindicatorrc.C,v 1.7.4.1 2004/04/05 15:19:46 bp Exp $ */
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

#include "errorestimator.h"
#include "directerrorindicatorrc.h"
#include "domain.h"
#include "conTable.h"
#include "mathfem.h"
#include "timestep.h"


void
DirectErrorIndicatorRC :: giveNodeChar (int inode, TimeStep* tStep, double& indicatorVal, double& currDensity) {

 int isize, i;
 const IntArray *con;
 Domain* d = this->giveDomain();
 ConnectivityTable* ct = d->giveConnectivityTable();
 DirectErrorIndicatorRCInterface* interface;

 con = ct->giveDofManConnectivityArray (inode);
 isize = con->giveSize();
 
 for (i=1; i<=isize; i++) {
  // ask for indicator variable value and determine current mesh density
  interface = (DirectErrorIndicatorRCInterface*) 
   d->giveElement(con->at(i))->giveInterface (DirectErrorIndicatorRCInterfaceType);
  if (!interface) {
    OOFEM_WARNING2 ("DirectErrorIndicatorRC::giveRequiredDofManDensity: element %d does not support DirectErrorIndicatorRCInterface", con->at(i));
  }
  if (i==1) {
   indicatorVal = ee->giveElementError (indicatorET, d->giveElement(con->at(i)), tStep);
   currDensity = interface->DirectErrorIndicatorRCI_giveCharacteristicSize ();
  } else {
   indicatorVal = max (indicatorVal, ee->giveElementError (indicatorET, d->giveElement(con->at(i)), tStep));
   currDensity =  max (currDensity, interface->DirectErrorIndicatorRCI_giveCharacteristicSize ());
  }
 } 
}


double
DirectErrorIndicatorRC :: giveDofManDensity (int num) {

 int isize, i;
 double currDensity=0.0;
 const IntArray *con;
 Domain* d = this->giveDomain();
 ConnectivityTable* ct = d->giveConnectivityTable();
 DirectErrorIndicatorRCInterface* interface;

 con = ct->giveDofManConnectivityArray (num);
 isize = con->giveSize();
 
 for (i=1; i<=isize; i++) {
  // ask for indicator variable value and determine current mesh density
  interface = (DirectErrorIndicatorRCInterface*) 
   d->giveElement(con->at(i))->giveInterface (DirectErrorIndicatorRCInterfaceType);
  if (!interface) {
    OOFEM_WARNING2 ("DirectErrorIndicatorRC::giveRequiredDofManDensity: elem %d does not support DirectErrorIndicatorRCInterface", con->at(i));
  }
  if (i==1) {
   currDensity = interface->DirectErrorIndicatorRCI_giveCharacteristicSize ();
  } else {
   currDensity =  max (currDensity, interface->DirectErrorIndicatorRCI_giveCharacteristicSize ());
  }
 } 
 return currDensity;
}



int
DirectErrorIndicatorRC :: estimateMeshDensities(TimeStep* tStep) {
 
 Domain* d = this->giveDomain();
 int inode, nnodes = d->giveNumberOfDofManagers();
 double indicatorVal, currDensity, proposedDensity;

 if (stateCounter == tStep->giveSolutionStateCounter()) return 1;

 this->currStrategy  = NoRemeshing_RS;
 this->nodalDensities.resize (nnodes);
   
 for (inode = 1; inode <= nnodes; inode ++) {
  
  this->giveNodeChar (inode, tStep, indicatorVal, currDensity);

  if (indicatorVal < minIndicatorLimit) {
   this->nodalDensities.at(inode) = zeroIndicatorDensity; //currDensity;
  } else if (indicatorVal >= maxIndicatorLimit) {
   this->nodalDensities.at(inode) = proposedDensity = maxIndicatorDensity;
   
   if (proposedDensity < (currDensity*this->remeshingDensityRatioToggle)) {
    this->currStrategy  = RemeshingFromCurrentState_RS;
   }

  } else {
   // evaluate the required size
   proposedDensity = minIndicatorDensity + 
    (indicatorVal-minIndicatorLimit) * (maxIndicatorDensity-minIndicatorDensity)/(maxIndicatorLimit-minIndicatorLimit);
   //proposedDensity = min (currDensity, proposedDensity);
   this->nodalDensities.at(inode) = proposedDensity;
   
   if (proposedDensity < (currDensity*this->remeshingDensityRatioToggle)) {
    this->currStrategy  = RemeshingFromCurrentState_RS;
   }
  }
 }
 // remember time stamp
 stateCounter = tStep->giveSolutionStateCounter();
 return 1;
}


double
DirectErrorIndicatorRC::giveRequiredDofManDensity (int num, TimeStep* tStep, int relative) {
 
 this->estimateMeshDensities (tStep);
 if (relative) {
  double currDens;
  currDens = this->giveDofManDensity (num);
  return this->nodalDensities.at(num) / currDens;
 } else return this->nodalDensities.at(num);
}
 

IRResultType
DirectErrorIndicatorRC :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 IR_GIVE_FIELD (ir, minIndicatorLimit, IFT_DirectErrorIndicatorRC_minlim, "minlim"); // Macro
 IR_GIVE_FIELD (ir, maxIndicatorLimit, IFT_DirectErrorIndicatorRC_maxlim, "maxlim"); // Macro
 IR_GIVE_FIELD (ir, minIndicatorDensity, IFT_DirectErrorIndicatorRC_mindens, "mindens"); // Macro
 IR_GIVE_FIELD (ir, maxIndicatorDensity, IFT_DirectErrorIndicatorRC_maxdens, "maxdens"); // Macro
 IR_GIVE_FIELD (ir, zeroIndicatorDensity, IFT_DirectErrorIndicatorRC_defdens, "defdens"); // Macro

 remeshingDensityRatioToggle = 0.80;
 IR_GIVE_OPTIONAL_FIELD (ir, remeshingDensityRatioToggle, IFT_DirectErrorIndicatorRC_remeshingdensityratio, "remeshingdensityratio"); // Macro

 return IRRT_OK;

}

RemeshingStrategy 
DirectErrorIndicatorRC::giveRemeshingStrategy (TimeStep* tStep)
{
 this->estimateMeshDensities (tStep);
 return this->currStrategy;
}
