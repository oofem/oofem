/* $Header: /home/cvs/bp/oofem/oofemlib/src/integrationrule.C,v 1.7.4.1 2004/04/05 15:19:43 bp Exp $ */
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


// file integrationRule.C

#include "integrationrule.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"

// initialize class member 

IntegrationRule::IntegrationRule (int n, Domain* domain, int startIndx, int endIndx, bool dynamic)
: FEMComponent (n,domain)
{
  numberOfIntegrationPoints = 0;
  gaussPointArray = NULL;
  firstLocalStrainIndx = startIndx;
  lastLocalStrainIndx  = endIndx;
  isDynamic = dynamic;
}


IntegrationRule::~IntegrationRule ()
{
  this->clear();
}


void
IntegrationRule::clear()
{
  int i;
  
  if (gaussPointArray) {
    for (i=0 ; i<numberOfIntegrationPoints ; i++)
      delete gaussPointArray[i] ;
    delete gaussPointArray ;}
  gaussPointArray = NULL;
  numberOfIntegrationPoints = 0;
}

GaussPoint* 
IntegrationRule :: getIntegrationPoint (int i)
{
#  ifdef DEBUG
 if ((i < 0) || (i >= numberOfIntegrationPoints)) {
   _error2 ("IntegrationRule::getIntegrationPoint - request out of bounds (%d)",i);
 }
#  endif
 return gaussPointArray[i];
}

void
IntegrationRule :: printOutputAt (FILE* file, TimeStep* stepN)
// Performs end-of-step operations.
{
 int i;

    for (i=0 ; i < numberOfIntegrationPoints ; i++)
    gaussPointArray[i] -> printOutputAt(file, stepN);
}

void
IntegrationRule :: updateYourself (TimeStep* tStep)
{
 // Updates the receiver at end of step.
 int i;

 for (i=0 ; i < numberOfIntegrationPoints ; i++)
  gaussPointArray[i] -> updateYourself(tStep) ;
}

 
void
IntegrationRule :: initForNewStep ()
{
 // initializes receiver to new time step or can be used
 // if current time step must be restarted
 // 
 // call material->initGpForNewStep() for all GPs.
 //

 int i;
 
 for (i=0 ; i < numberOfIntegrationPoints ; i++)
  gaussPointArray[i]->giveMaterial()->initGpForNewStep(gaussPointArray[i]);
}


contextIOResultType
IntegrationRule :: saveContext (FILE* stream, void *obj)
{
 //
 // saves full  context (saves state variables, that completely describe
 // current state)
 //

  contextIOResultType iores;
  int         i ;
  GaussPoint* gp ;

  if (stream == NULL) _error ("saveContex : can't write into NULL stream");

  if (!isDynamic) {
    for (i=0 ; i < numberOfIntegrationPoints ; i++) {
      gp = gaussPointArray[i] ;
      if ((iores = gp->giveCrossSection()->saveContext(stream,gp)) != CIO_OK) THROW_CIOERR(iores);
    }
  } else {
    // write size 
    if (fwrite(&numberOfIntegrationPoints,sizeof(int),1,stream)!=1) THROW_CIOERR (CIO_IOERR);
    for (i=0 ; i < numberOfIntegrationPoints ; i++) {
      gp = gaussPointArray[i] ;
      // write gp weight, coordinates, element number, and material mode
      double dval = gp->giveWeight();
      if (fwrite(&dval,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
      if ((iores = gp->giveCoordinates()->storeYourself(stream))!=CIO_OK) THROW_CIOERR(iores);
      int ival = gp->giveElement()->giveNumber();
      if (fwrite(&ival,sizeof(int),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
      MaterialMode m = gp->giveMaterialMode();
      if (fwrite(&m,sizeof(MaterialMode),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
      // write gp data
      if ((iores = gp->giveCrossSection()->saveContext(stream,gp)) != CIO_OK) THROW_CIOERR(iores);
    }
  }    
  return CIO_OK;
}

contextIOResultType
IntegrationRule :: restoreContext (FILE* stream, void *obj)
{
 //
 // restores full element context (saves state variables, that completely describe
 // current state)
 //
 
  contextIOResultType iores;
 int         i ;
  GaussPoint* gp ;
 
  if (stream == NULL) _error ("restoreContex : can't write into NULL stream");

  if (!isDynamic) {
    for (i=0 ; i < numberOfIntegrationPoints ; i++) {
      gp = gaussPointArray[i] ;
      if ((iores = gp -> giveCrossSection()->restoreContext(stream,gp)) != CIO_OK) THROW_CIOERR(iores);
    }
  } else {
    // read size
    bool __create = false;
    int size;

    if (fread (&size,sizeof(int),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
    if (numberOfIntegrationPoints != size) {
      this->clear();
      __create = true;
      gaussPointArray = new GaussPoint* [size];
      numberOfIntegrationPoints = size;

    }
    
    for (i=0 ; i < numberOfIntegrationPoints ; i++) {
      // read weight
      double w;
      if (fread (&w,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
      // read coords
      FloatArray c;
      if ((iores = c.restoreYourself(stream))!= CIO_OK) THROW_CIOERR(iores);
      // read element number and material mode
      int n;
      if (fread(&n,sizeof(int),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
      MaterialMode m;
      if (fread(&m,sizeof(MaterialMode),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
      if (__create) {
        (gaussPointArray)[i] = new GaussPoint(domain->giveElement(n),i+1,c.GiveCopy(),w,m);
      } else {
        gp = gaussPointArray[i] ;
        gp->setWeight (w);
        gp->setCoordinates(c);
        gp->setElement (domain->giveElement(n));
        gp->setMaterialMode (m);
      }
      // read gp data
      gp = gaussPointArray[i] ;
      if ((iores = gp -> giveCrossSection()->restoreContext(stream,gp)) != CIO_OK) THROW_CIOERR(iores);
    }
  }
 
  return CIO_OK;
}


int 
IntegrationRule :: setUpIntegrationPoints (integrationDomain mode, int nPoints, Element* elem,
               MaterialMode matMode)
{
  switch (mode) {
  case _Line:
    return  (numberOfIntegrationPoints = this->SetUpPointsOnLine (nPoints , elem, matMode, &gaussPointArray)) ;

  case _Triangle:
    return  (numberOfIntegrationPoints = this->SetUpPointsOnTriagle (nPoints , elem, matMode, &gaussPointArray)) ;

  case _Square:
    return  (numberOfIntegrationPoints = this->SetUpPointsOnSquare  (nPoints , elem, matMode, &gaussPointArray)) ;

  case _Cube:
    return  (numberOfIntegrationPoints = this->SetUpPointsOnCube    (nPoints , elem, matMode, &gaussPointArray)) ;

  case _Tetrahedra:
    return  (numberOfIntegrationPoints = this->SetUpPointsOnTetrahedra (nPoints , elem, matMode, &gaussPointArray)) ;

  default:
    _error ("IntegrationRule::setUpIntegrationPoints - unknown mode");
  }
  return 0;
}

int 
IntegrationRule::setUpEmbeddedIntegrationPoints (integrationDomain mode, int nPoints, Element *elem, MaterialMode matMode,
						 const FloatArray** coords)
{
  switch (mode) {
  case _Embedded2dLine:
    return  (numberOfIntegrationPoints = this->SetUpPointsOn2DEmbeddedLine (nPoints , elem, matMode, &gaussPointArray, coords)) ;

  default:
    _error ("IntegrationRule::setUpEmbeddedIntegrationPoints - unknown mode");
  }
  return 0;
}


