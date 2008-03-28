/* $Header: /home/cvs/bp/oofem/oofemlib/src/Attic/lobattoir.h,v 1.1.2.1 2004/04/05 15:19:43 bp Exp $ */
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


//
// class lobattoIntegrationRule
//

#ifndef lobattoir_h
#define lobattoir_h

#include "integrationrule.h"

/**
 Class representing Lobatto-quadrature integration rule. 
 The number of integration points and their coordinates and integration weights depends on
 integration rule type (rule for intagration in 1d, 2d, 3d) and required  acurracy.
 */
class LobattoIntegrationRule : public IntegrationRule
{
/*
DESCRIPTION:
   Implements integration rule class.
  Stores integration points used for integration
  of necesary terms (for example computation of  stiffness matrix 
  or computation of element nodal force vector ) 
  and it  corresponds to some local strains 
  on finite element level. Finite element can have many 
  integration rules corresponding to  different strains.

TASKS:
   instanciating yourself
  returning number of integration points used
  returning requested integration point - method getIntegrationPoint
  returning inteval of components (i.e.of local strain vector), where apply
  printing yourself
  updating yourself
  initializing for new time step
  saving & restoring context
*/
private:

public:
 /**
  Constructor.
  @param n number associated with receiver
  @param domain reference to domain.
  @param startIndx first component, for which rule applies
  @param endIndx last component, for which rule applies
  @param dynamic flag indicating that receiver can change
  */
 LobattoIntegrationRule (int , Element*, int, int, bool dynamic);
 /// Destructor
  ~LobattoIntegrationRule();

 ///Returns classType id of receiver.
  classType giveClassID () const     {return LobattoIntegrationRuleClass;}
 ///Returns class name of the receiver.
  const char*  giveClassName () const {return "LobattoIntegrationRule" ;}
 IRResultType initializeFrom (InputRecord* ir) {return IRRT_OK;}

 /**
  Returns requred number of integration points to exactly integrate
  polynomial of order approxOrder on given domain.
  When approxOrder is too large and is not supported by implementation
  method returns -1.
  */
 int getRequiredNumberOfIntegrationPoints (integrationDomain dType, int approxOrder) ;

protected:
 /**
  Sets up receiver's  integration points on unit line integration domain.
  @returns number of integration points. 
  */
  int  SetUpPointsOnLine    (int , Element*, MaterialMode, GaussPoint***) ;
 /**
  Sets up receiver's  integration points on triangular (area coords) integration domain.
  @returns number of integration points.
  */
  int  SetUpPointsOnTriagle (int , Element*, MaterialMode, GaussPoint***) ;
 /**
  Sets up receiver's  integration points on unit square integration domain.
  @returns number of integration points.
  */
  int  SetUpPointsOnSquare  (int , Element*, MaterialMode, GaussPoint***) ;
 /**
  Sets up receiver's  integration points on unit cube integration domain.
  @returns number of integration points.
  */
  int  SetUpPointsOnCube    (int , Element*, MaterialMode, GaussPoint***) ;
 /**
  Sets up receiver's  integration points on tetrahedra (volume coords) integration domain.
  @returns number of integration points.
  */
  int  SetUpPointsOnTetrahedra (int , Element*, MaterialMode, GaussPoint***) ;
  
};

#endif // lobattoir_h
