/* $Header: /home/cvs/bp/oofem/sm/src/zzerrorestimator.h,v 1.5 2003/04/06 14:08:32 bp Exp $ */
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

//   *************************************************************************
//   *** CLASS ERROR ESTIMATOR (INDICATOR) ACORDING TO ZIENKIEWICZ AND ZHU ***
//   *************************************************************************

#ifndef zzerrorestimator_h 

#include "errorestimator.h"
#include "interface.h"


#define ZZErrorEstimator_ElementResultCashed

class Element;
class GaussPoint;

/**
 The implementation of Zienkiewicz Zhu Error Estimator.
 The basic task is to evaluate the stress error on associated domain. 
 The algorithm is written in general way, so it is possible to to evaluate 
 different errors (for example temperature error). Then corresponding
 member attribute identifying the type of quantity used should be declared and initialized
 (for example in instanciateYourself() service). Then the modification is required 
 only when requesting element contributions.

 This task requires the special element algorithms, which are supported at element level 
 using interface concept.
 This estimator also provides the compatible Remeshing Criteria, which 
 based on error measure will evaluate the required mesh density of a new domain.
*/
class ZZErrorEstimator : public ErrorEstimator {

public:
 // type of norm used
 enum NormType {L2Norm, EnergyNorm};
  // nodal recovery type
  enum NodalRecoveryType {ZZRecovery, SPRRecovery};

protected:
 /// global error norm
 double globalENorm;  
  /// global norm of quantity which error is evaluated
  double globalSNorm;
#ifdef ZZErrorEstimator_ElementResultCashed
  /// cache storing element norms 
  FloatArray eNorms;
#endif
  /// type of norm used
  NormType normType;
  /// nodal recovery type
  NodalRecoveryType nodalRecoveryType;

 /// actual state counter.
 StateCounterType        stateCounter;

public:
 /// Constructor
 ZZErrorEstimator (int n, Domain *d) : ErrorEstimator(n,d) {eeType = EET_ZZEE; stateCounter = 0; 
                               normType = L2Norm; nodalRecoveryType = ZZRecovery;}
 /// Destructor
 ~ZZErrorEstimator() {}
  /** Returns the element error of requested type. The estimateError service should be called before.
  @param type error type
  @param elem element for which error requested
  @param tStep time step
  */
 virtual double giveElementError (EE_ErrorType type, Element* elem, TimeStep* tStep);
 /** Returns the characteristic value of given type. 
  The estimateError service should be called before. Intended to be used by remeshingCriterias to querry 
  various values provided by specific error estimator.
  @param type value type
  @param tStep time step
  */
 virtual double giveValue (EE_ValueType type, TimeStep* tStep);

 /**
  Estimates the error on associated domain at given timeSte.
  @param tStep time step
  */
 virtual int estimateError (EE_ErrorMode mode, TimeStep* tStep);
 /** Returns reference to associated remeshing criteria.
  */
 virtual RemeshingCriteria* giveRemeshingCrit ();
 /** Initializes receiver acording to object description stored in input record.
  This function is called immediately after creating object using
  constructor. InitString can be imagined as data record in component database
  belonging to receiver. Receiver may use value-name extracting functions 
  to extract particular field from record. 
  @see readInteger, readDouble and similar functions */
 virtual IRResultType initializeFrom (InputRecord* ir);
 /// Returns class name of the receiver.
 const char* giveClassName () const { return "ZZErrorEstimator" ;}
 /** Returns classType id of receiver.
  @see FEMComponent::giveClassID 
  */
 classType                giveClassID () const { return ZZErrorEstimatorClass; }

protected:
};


/**
 The element interface corresponding to ZZErrorEstimator.
 It declares necessary services provided by element to be compatible with ZZErrorEstimator.
*/
class ZZErrorEstimatorInterface : public Interface
{
public:
 /// Constructor
  ZZErrorEstimatorInterface () {}

  /// Returns reference to corresponding element 
  virtual Element* ZZErrorEstimatorI_giveElement() = 0;
  /// Computes the interpolation matrix used for recovered values
  virtual void ZZErrorEstimatorI_computeEstimatedStressInterpolationMtrx (FloatMatrix& answer, GaussPoint* gp, 
                                      InternalStateType) = 0;
  /**
   Computes the element contributions to global norms.
   @param eNorm element contribution to error norm
   @param sNorm element contribution to norm of quantity which error is evaluated
   @param norm determines the type of norm to evaluate
   @param tStep time Step
   */
  virtual void ZZErrorEstimatorI_computeElementContributions (double& eNorm, double& sNorm, ZZErrorEstimator::NormType, 
                                InternalStateType type, TimeStep*);
};


/**
 The class representing Zienkiewicz-Zhu remeshing criteria.
 (Assumes that error is equally distributed between elements, then the requirement for max. permissible error
 can be translated into placing a limit on the error on each element.)
 The basic task is to evaluate the required mesh density (at nodes) on given domain,
 based on informations provided by the compatible error ertimator. 
 

 The remeshing criteria is maintained by the corresponding error estimator. This is mainly due to fact, that is
 necessary for given EE to create compatible RC. In our concept, the EE is responsible. 
 
*/
class ZZRemeshingCriteria : public RemeshingCriteria {

public:
 /// mode of receiver, allows to use it in more general situations
 enum ZZRemeshingCriteriaModeType {stressBased};

protected:
 /// Array of nodal mesh densities
 FloatArray nodalDensities;
  /// Remeshing strategy proposed
  RemeshingStrategy remeshingStrategy;
 /// actual values (densities) state counter.
 StateCounterType        stateCounter;
  /// mode of receiver
  ZZRemeshingCriteriaModeType mode;
  /// required error to obtain
  double requiredError;
  /// minimum element size alloved
 double minElemSize;


public:
 /// Constructor
 ZZRemeshingCriteria (int n, ErrorEstimator *e) ;
  /// Destructor
  ~ZZRemeshingCriteria() {}

 /** Returns the required mesh size n given dof manager.
  The mesh density is defined as a required element size 
  (in 1D the element length, in 2D the square from element area).
  @param num dofman  number
  @param tStep time step
  @param relative if zero, then actual density is returned, otherwise the relative density to current is returned.
  */
 virtual double giveRequiredDofManDensity (int num, TimeStep* tStep, int relative=0);
 /**
   Returns existing mesh size for given dof manager.
   @param num dofMan number
 */
 virtual double giveDofManDensity (int num);
  /**
  Determines, if the remeshing is needed, and if nedded, the type of strategy used
  */
  virtual RemeshingStrategy giveRemeshingStrategy (TimeStep* tStep);
 /**
  Estimates the nodal densities.
  @param tStep time step
  */
 virtual int estimateMeshDensities (TimeStep* tStep);
 /** Initializes receiver acording to object description stored in input record.
  This function is called immediately after creating object using
  constructor. InitString can be imagined as data record in component database
  belonging to receiver. Receiver use value-name extracting functions 
  to extract particular field from record. 
  @see readInteger, readDouble and similar functions */
 virtual IRResultType initializeFrom (InputRecord* ir);

  /// Returns "ZZErrorEstimator" - class name of the receiver.
 const char* giveClassName () const { return "ZZErrorEstimator" ;}
 /** Returns ZZRemeshingCriteriaClass - classType id of receiver.
  @see FEMComponent::giveClassID 
  */
  classType                giveClassID () const { return ZZRemeshingCriteriaClass; }

protected:
  
};

/**
 The corresponding element interface to ZZRemeshingCriteria class.
 Declares the necessary services, which have to be provided by particular elements.
*/

class ZZRemeshingCriteriaInterface : public Interface {
 
public:
 /// Constructor
 ZZRemeshingCriteriaInterface() : Interface () {}
 /**
   Determines the characteristic size of element. This quantity is defined as follows:
   For 1D it is the element length, for 2D it is the square root of element area.
  */
 virtual double ZZRemeshingCriteriaI_giveCharacteristicSize () = 0;
 /**
  Returns the polynomial order of receiver trial functions.
  */
 virtual int ZZRemeshingCriteriaI_givePolynOrder () = 0;

};




#define zzerrorestimator_h
#endif

