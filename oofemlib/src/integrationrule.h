/* $Header: /home/cvs/bp/oofem/oofemlib/src/integrationrule.h,v 1.9 2003/04/06 14:08:24 bp Exp $ */
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
// class IntegrationRule
//

#ifndef integrationrule_h
#define integrationrule_h

#include "cltypes.h"
#include "element.h"
#include "femcmpnn.h"

class GaussPoint;

/**
 Abstract base class representing integration rule. The integration rule is 
 a collection of integration points used to  numerically integrate some formula.
 The number of integration points and their coordinates and integration weights depends on
 integration rule type (rule for intagration in 1d, 2d, 3d) and required  acurracy.
 General services for inicialization are declared. Services for integration point retrieval are provided.
 
 In general, finite elements can have multiple integration rules, for diferrent tasks or
 when some components are integrated using reduced or selective integration.
 Therefore, first and last index variables are introduced to characterize components
 for which given integration rule applies.

 The integration rule is a rather passive object. It does not perform numerical
 integration - it just provide way how to set up correct integration points 
 and weights. 

 Because integration points contain related history parameters (using matarial status),
 the unique copy of integration rule must exist on each element. The integration rule is 
 exclusively possessed by particular finite element.
*/
class IntegrationRule 
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
-   instanciating yourself
-  returning number of integration points used
-  returning requested integration point - method getIntegrationPoint
-  returning inteval of components (i.e.of local strain vector), where apply
-   returning array of gauss points, acording to specific 
     integration rule (Gauss -rule, Newton-Cortes rule ...).
     integration points and corresponding weights are stored in 
     Gauss point class.

  printing yourself
  updating yourself
  initializing for new time step
  saving & restoring context


REMARK
   The Integrator is a rather passive object : it does not perform numerical
   integration - it just provide way how to set up correct integration points 
   and weights. Integration is performed by elements.
*/
private:

  /// number
  int number;
  /// Pointer to element  
  Element* elem;

 /// Array containing integration points
 GaussPoint** gaussPointArray;
 /// Number of integration point of receiver
 int numberOfIntegrationPoints;
 /**
  firstLocalStrainIndx and lastLocalStrainIndx indexes describe range of components (strains for example)
  for which receiver intagration points apply.
  */
 int firstLocalStrainIndx, lastLocalStrainIndx;

 /** flag indicating that rule is dynamic, ie, its gauss points (their number, coordinates, weights) can change during
     computation. Then some more data should be stored/restored from context file to reflect such dynamic feature */
 bool isDynamic;
public:

 /**
  Constructor.
  @param n number associated with receiver
  @param domain reference to domain.
  @param startIndx first component, for which rule applies
  @param endIndx last component, for which rule applies
  @param dynamic flag indicating that receiver can change
  */
  IntegrationRule (int n, Element* e, int startIndx, int endIndx, bool dynamic);
  IntegrationRule (int n, Element* e);
 /// Destructor.
  virtual ~IntegrationRule();
  
 /**
  Returns number of integration points of receiver.
  */
  int getNumberOfIntegrationPoints () {return numberOfIntegrationPoints;}
 /**
  Access particular integration point of receiver.
  @param n integration point number (should be in range 0,.., getNumberOfIntegrationPoints()-1).
  */
  GaussPoint* getIntegrationPoint  (int n);
 /**
  Returns starting component index, for which receiver applies.
  */
  int getStartIndexOfLocalStrainWhereApply () {return firstLocalStrainIndx;}
 /**
  Returns last component index, for which receiver applies.
  */
  int getEndIndexOfLocalStrainWhereApply () {return lastLocalStrainIndx;}
 /**
  Initializes the receiver. Receiver integration points are created acording to given parameters.
  @param mode describes integration domain
  @param nPoints required number of integration points of receiver
  @param matMode material mode of receiver's intagration points
  @return nPoints
  */
  int setUpIntegrationPoints (integrationDomain mode, int nPoints, MaterialMode matMode);
  /**
  Initializes the receiver. Receiver integration points are created acording to given parameters.
  @param mode describes integration domain
  @param nPoints required number of integration points of receiver
  @param matMode material mode of receiver's intagration points
  @return nPoints
  */
  int setUpEmbeddedIntegrationPoints (integrationDomain mode, int nPoints, MaterialMode matMode,
				      const FloatArray** coords);

 /**
  Prints receiver's output to given stream.
  Invokes printOutputAt service on all receiver's integration points.
  */
  virtual void printOutputAt (FILE* file, TimeStep* stepN);
 /**
  Updates receiver state.
   Calls updateYourself service of all receiver's integration points.
  */
  void updateYourself (TimeStep* tStep);
 /**
  Initializes receiver.
  Calls initForNewStep service of all receiver's integration points.
  */
  void initForNewStep ();

  /** Returns reference to element containing receiver */
  Element* giveElement () {return elem;}

  /**
     Abstract service.
     Returns requred number of integration points to exactly integrate
     polynomial of order approxOrder on given domain.
     When approxOrder is too large and is not supported by implementation
     method returns -1. Must be overloaded by derived classes.
  */
  virtual int getRequiredNumberOfIntegrationPoints (integrationDomain dType, int approxOrder) {return 0;}
  
 /**
  Saves receiver's context to stream.
  Calls saveContext service for all receiver's integration points.
  Note: does not call the FEMComponent::saveContext service, in order not
  to write class id info for each integration rule.
  @exception throws an ContextIOERR exception if error encountered.
  */
  virtual contextIOResultType saveContext (DataStream* stream, ContextMode mode, void *obj);
 /**
  Restores receiver's context to stream.
  Calls restoreContext service for all receiver's integration points.
  Note: does not call the FEMComponent::restoreContext service, in order not
  to write class id info for each integration rule.
  @param obj should be a pointer to invoking element, ie., to which the receiver will belong to.
  @exception throws an ContextIOERR exception if error encountered.
  */
  virtual contextIOResultType restoreContext (DataStream* stream, ContextMode mode, void *obj);
  /**
     Clears the receiver, ie dealocates all integration points
  */
  void clear ();

  ///Returns classType id of receiver.
  virtual classType giveClassID () const  {return IntegrationRuleClass;}
  /// Returns class name of the receiver.
  virtual const char*  giveClassName () const {return "IntegrationRule" ;}
  virtual IRResultType initializeFrom (InputRecord* ir) {return IRRT_OK;}

protected:
 /** 
  Sets up receiver's  integration points on unit line integration domain.
  Default implementaion does not sets up any integration points and returns 0.
  Must be overloaded by deived classes.
  @returns number of integration points. 
  */
  virtual int  SetUpPointsOnLine       (int , Element*, MaterialMode, GaussPoint***) {return 0;}
 /** 
  Sets up receiver's  integration points on triangular (area coords) integration domain.
  Default implementaion does not sets up any integration points and returns 0.
  Must be overloaded by deived classes.
  @returns number of integration points.
  */
  virtual int  SetUpPointsOnTriagle    (int , Element*, MaterialMode, GaussPoint***) {return 0;}
 /** 
  Sets up receiver's  integration points on unit square integration domain.
  Default implementaion does not sets up any integration points and returns 0.
  Must be overloaded by deived classes.
  @returns number of integration points.
  */
  virtual int  SetUpPointsOnSquare     (int , Element*, MaterialMode, GaussPoint***) {return 0;}
 /** 
  Sets up receiver's  integration points on unit cube integration domain.
  Default implementaion does not sets up any integration points and returns 0.
  Must be overloaded by deived classes.
  @returns number of integration points.
  */
  virtual int  SetUpPointsOnCube       (int , Element*, MaterialMode, GaussPoint***) {return 0;}
 /** 
  Sets up receiver's  integration points on tetrahedra (volume coords) integration domain.
  Default implementaion does not sets up any integration points and returns 0.
  Must be overloaded by deived classes.
  @returns number of integration points.
  */
  virtual int  SetUpPointsOnTetrahedra (int , Element*, MaterialMode, GaussPoint***) {return 0;}
  /**
     Sets up integration points on 2D embedded line inside 2D volume (the list of local coordinates 
     should be provided).
  */
  virtual int SetUpPointsOn2DEmbeddedLine (int nPoints, Element* elem, MaterialMode mode, GaussPoint***, 
                                           const FloatArray** coords) {return 0;}

};

#endif
