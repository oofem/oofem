/* $Header: /home/cvs/bp/oofem/oofemlib/src/material.C,v 1.11.4.1 2004/04/05 15:19:43 bp Exp $ */
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


//   file MATERIAL.CC

#include "material.h"
#include "crosssection.h"
#include "domain.h"
#include "verbose.h"
#include "isolinearelasticmaterial.h"
//#include "ortholinearelasticmaterial.h"
//#include "perfectlyplasticmaterial.h"
//#include "steel1.h"
//#include "concrete2.h"
//#include "concrete3.h"
//#include "cebfip78.h"
//#include "doublepowerlaw.h"
//#include "b3mat.h"
#include "isolinearheatmat.h"

#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "cltypes.h"
#include "mathfem.h"
#include "usrdefsub.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <math.h>
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
#endif


void
Material :: giveCharacteristicMatrix (FloatMatrix& answer,
                   MatResponseForm form, MatResponseMode rMode,
                   GaussPoint* gp,
                   TimeStep* atTime)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
  _error ("Material::giveCharacteristicMatrix is fully abstract, no implementation");
  return ;
}


double
Material :: giveCharacteristicValue (MatResponseMode rMode,
                                     GaussPoint* gp,
                                     TimeStep* atTime)
//
// Returns characteristic value of the receiver
//
{
  _error ("Material :: giveCharacteristicValue is purely abstract");
  return 0.0;
}


double  
Material :: give (int aProperty)
   // Returns the value of the property aProperty (e.g. the Young's modulus
   // 'E') of the receiver.
 // atTime allows time dependent behaviour to be taken into account
{
   double  value = 0.0;
 
   if (propertyDictionary -> includes(aProperty))
   value = propertyDictionary -> at(aProperty) ;
   else {                                         
     _error ("give: property not defined");
   }
   return value ;
}


IRResultType
Material :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 double value ;
 
#  ifdef VERBOSE
 // VERBOSE_PRINT1 ("Instanciating material ",this->giveNumber()) 
#  endif
 
 IR_GIVE_FIELD (ir, value, IFT_Material_density, "d"); // Macro
 propertyDictionary -> add('d',value);
 
 this->castingTime = -1.e6;
 IR_GIVE_OPTIONAL_FIELD (ir, castingTime, IFT_Material_castingtime, "castingtime"); // Macro
 
 return IRRT_OK;
}



int
Material::giveInputRecordString(std::string &str, bool keyword)
{
 char buff[1024];

 FEMComponent::giveInputRecordString(str, keyword);
 sprintf(buff, " d %e", this -> give('d'));
 str += buff;

 return 1;
}


int
Material :: hasMaterialModeCapability (MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
 return 0;
}


void  
Material :: printYourself ()
   // Prints the receiver on screen.
{
  printf ("Material with properties : \n") ;
  propertyDictionary -> printYourself() ;
}



Material* 
Material :: ofType (char* aClass)
// Returns a new material, which has the same number than the receiver,
// but belongs to aClass (Hook, Microplane, ...).
{
  Material* newMaterial ;
  
  if (! strncasecmp(aClass,"isole",5))
  newMaterial = new IsotropicLinearElasticMaterial (this->giveNumber(),domain);
/*  else if (! strncasecmp(aClass,"orthole",7))
  newMaterial = new OrthotropicLinearElasticMaterial(this->giveNumber(),domain);
  else if (! strncasecmp(aClass,"steel1",6))
  newMaterial = new Steel1 (number,domain) ;
  else if (! strncasecmp(aClass,"concrete2",9))
  newMaterial = new Concrete2 (number,domain) ;
  else if (! strncasecmp(aClass,"concrete3",9))
  newMaterial = new Concrete3 (number,domain) ;
  else if (! strncasecmp(aClass,"cebfip78",8))
  newMaterial = new CebFip78Material (number,domain) ;
  else if (! strncasecmp(aClass,"doublepowerlaw",14))
  newMaterial = new DoublePowerLawMaterial (number,domain) ;
  else if (! strncasecmp(aClass,"b3mat",5))
  newMaterial = new B3Material (number,domain) ;  */
  else if (! strncasecmp(aClass,"isolhm",6))
    newMaterial = new IsotropicLinearHeatMaterial (number,domain) ;
  
  else {
  // last resort - call aditional user defined subroutine
  newMaterial = ::CreateUsrDefMaterialOfType (aClass,number,domain);
  if (newMaterial == NULL) {
   
    _error2 ("ofType:  unknown material (%s)\n", aClass) ;
    exit(0) ;}
 }
 return newMaterial ;
}




// 
// store & restore context - materialinfo in gp not saved now!
//

contextIOResultType
Material :: saveContext (DataStream* stream, ContextMode mode, void *obj)
//
// saves full material status (saves state variables, that completely describe
// current state) stored in gp->matstatusDict with key =  (int)this->giveClassID()
// storing of correspondin context if it is defined for current material in 
// gp status dictionary should be performed here by overloading this function.
// (such code should invoke also cooresponding function for yield conditions,
//  submaterials and so on)
//

//
{
 contextIOResultType iores;

 if ((iores = FEMComponent::saveContext (stream, mode, obj)) != CIO_OK) THROW_CIOERR(iores);
 // corresponding gp is passed in obj
 GaussPoint* gp = (GaussPoint*) obj;
 if (gp == NULL) THROW_CIOERR(CIO_BADOBJ);
 
 // write raw data - we save status there for this
 MaterialStatus* status =  this->giveStatus(gp);

 if (status) if ((iores = status->saveContext (stream, mode, obj)) != CIO_OK) THROW_CIOERR(iores);

 return CIO_OK;
}

contextIOResultType
Material :: restoreContext (DataStream* stream, ContextMode mode, void *obj)
//
// restores full material status (saves state variables, that completely describe
// current state) stored in gp->matstatusDict with key =  (int)this->giveClassID()
// resstoring of correspondin context if it is defined for current material in 
// gp status dictionary should be performed here by overloading this function.
// (such code should invoke also cooresponding function for yield conditions,
//  submaterials and so on)
//

//
{
 contextIOResultType iores;
 // invoke base class service
 if ((iores = FEMComponent::restoreContext (stream, mode, obj)) != CIO_OK) THROW_CIOERR(iores);

 // corresponding gp is passed in obj
 GaussPoint* gp = (GaussPoint*) obj;
 if (gp == NULL) THROW_CIOERR(CIO_BADOBJ);
 
 // read raw data - context
 MaterialStatus* status =  this->giveStatus(gp);
 if (status) if ((iores = status->restoreContext (stream, mode, obj)) != CIO_OK) THROW_CIOERR(iores);

 return CIO_OK;
}



MaterialStatus*
Material :: giveStatus (GaussPoint* gp) const
/* 
 returns material status in gp corresponding to specific material class
*/
{
 MaterialStatus* status;
 
 status = gp->giveMaterialStatus();
 if (status==NULL) {
  // create a new one
  status = this -> CreateStatus (gp);

  if(this -> giveClassID() != gp -> giveElement() -> giveMaterial() -> giveClassID()) {
    _warning2 ("giveStatus: Created material status at element %d is of different material type", 
              gp -> giveElement() -> giveNumber());
  }

  // if newly created status is null
  // dont include it. specific instance
  // does not have status.
  if (status != NULL) gp ->setMaterialStatus (status);
 }
 return status;
}


void 
Material :: initTempStatus (GaussPoint* gp)
//
// Initialize MatStatus (respective it's temporary variables at the begining 
// of integrating incremental constitutive relations) to correct values
//
{
 MaterialStatus *status = this -> giveStatus (gp);
 if (status)
  status -> initTempStatus ();
}
 

void
Material :: initGpForNewStep (GaussPoint* gp)
//
// initialize gp record at the begining of new load Increment
// initialize gp status using this->initTempStatus(gp);
//
// this means: kepping Stress, Strain PlasticStrain Vectors.
// zeroing Stress, Strain PlasticStrain Increment Vectors.
//
{
 // initialize status
 this->initTempStatus (gp);
}



void
Material :: updateYourself(GaussPoint* gp, TimeStep* atTime) 
//
//
// We call MaterialStatus->updateYourself()
//
{
 MaterialStatus* status = this->giveStatus(gp);
 if (status) status->updateYourself(atTime);
}



