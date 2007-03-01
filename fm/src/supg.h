/* $Header: /home/cvs/bp/oofem/tm/src/nltransienttransportproblem.h,v 1.1 2003/04/14 16:01:39 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2002   Borek Patzak                                       



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

#ifndef supg_h

#include "engngm.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"
#include "materialinterface.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

/**
 This class represents transient incopressible flow problem. Solution is based on 
 algorithm with SUPG/PSPG stabilization.
*/
class SUPG : public EngngModel
{ 
protected:
  /// Numerical method used to solve the problem
  SparseLinearSystemNM *nMethod;

  LinSystSolverType solverType;
  SparseMtrxType sparseMtrxType;
  
  SparseMtrx* lhs;
  PrimaryField* VelocityPressureField;
  //PrimaryField VelocityField;
  FloatArray accelerationVector; //, previousAccelerationVector;
  FloatArray incrementalSolutionVector; 
  
  double deltaT;
  int deltaTLTF;
  /// covergence tolerance
  double atolv, rtolv; 
  /// max number of iterations
  int maxiter;
  /** flag if set to true (default), then when max number of iteration reached, computation stops
      otherwise computation continues with next step
  */
  bool stopmaxiter;
  // integration constants
  double alpha;

  int initFlag;
  int consistentMassFlag;

  bool equationScalingFlag;
  // length, velocity, and density scales
  double lscale, uscale, dscale;
  // Reynolds number
  double Re;
  // indicates if equation renumbering requiered
  bool renumberFlag;


  // material interface representation for multicomponent flows
  MaterialInterface* materialInterface;
  // map of active dofmans for problems with free surface and only one fluid considered
  // IntArray __DofManActivityMask;
  // free surface flag -> we solve free surface problem by single reference fluid
  // int fsflag; 

public:
  /*  SUPG (int i, EngngModel* _master = NULL) : EngngModel (i,_master), VelocityPressureField(this,1,FBID_VelocityPressureField, EID_MomentumBalance_ConservationEquation, 1),accelerationVector()
   */
  SUPG (int i, EngngModel* _master = NULL) : EngngModel (i,_master), accelerationVector() {
    initFlag=1; lhs=NULL; ndomains=1;nMethod=NULL;
    VelocityPressureField = NULL;
    consistentMassFlag = 0;
    equationScalingFlag = false; lscale=uscale=dscale=1.0;renumberFlag=false;
    materialInterface = NULL;
  }
  ~SUPG () {
    if (VelocityPressureField) delete VelocityPressureField;
    if (materialInterface) delete materialInterface;
  }

  void solveYourselfAt (TimeStep *);
 /**
  Updates nodal values
  (calls also this->updateDofUnknownsDictionary for updating dofs unknowns dictionaries
  if model supports changes of static system). The element internal state update is also forced using
  updateInternalState service.
  */
  virtual void               updateYourself (TimeStep*) ;

  double   giveUnknownComponent ( EquationID, ValueModeType, TimeStep*, Domain*, Dof*);
  double   giveUnknownComponent ( UnknownType, ValueModeType, TimeStep*, Domain*, Dof*);
  /**
     Returns characteristic vector of element. The Element::GiveCharacteristicVector function
     should not be called directly, because EngngModel may require some special modification
     of characteristic vectors supported on element level. 
     @param answer characteristic matrix
     @param num element number
     @param type type of CharMatrix requsted
     @param tStep time step when response is computed
     @param domain source domain
  */
  virtual void giveElementCharacteristicVector (FloatArray& answer, int num, CharType type, ValueModeType mode, TimeStep* tStep, Domain* domain) ;
  /**
     Returns characteristic matrix of element. The Element::GiveCharacteristicMatrix function
     should not be called directly, because EngngModel may require some special modification
     of characteristic matrices supported on element level. But default implementation does 
     the direct call to element level. 
     @param answer characteristic matrix
     @param num element number
     @param type type of CharMatrix requsted
     @param tStep time step when response is computed
     @param domain source domain
  */
  virtual void giveElementCharacteristicMatrix (FloatMatrix& answer, int num, CharType type, TimeStep* tStep, Domain* domain) ;


  contextIOResultType saveContext (DataStream* stream, ContextMode mode, void *obj = NULL) ;
  contextIOResultType restoreContext (DataStream* stream, ContextMode mode, void *obj = NULL);
 
  void   updateDomainLinks();
  
  TimeStep* giveNextStep ();
  TimeStep*  giveSolutionStepWhenIcApply();
  NumericalMethod* giveNumericalMethod (TimeStep*);
  
  /// Initialization from given input record
  IRResultType initializeFrom (InputRecord* ir);

  // consistency check
  virtual int checkConsistency (); // returns nonzero if o.k.
  // identification
  const char* giveClassName () const { return "SUPG";}
  classType giveClassID ()      const { return SUPGClass;}
  fMode giveFormulation () { return TL; }

  /** DOF printing routine. Called by DofManagers to print Dof specific part.
      Dof class provides component printing routines, but emodel is responsible
      for what will be printed at DOF level.
      @param stream output stream
      @param iDof dof to be processed
      @param atTime solution step
  */
  virtual void printDofOutputAt (FILE* stream, Dof* iDof, TimeStep* atTime);

   virtual int                giveNumberOfEquations (EquationID);
   virtual int                giveNumberOfPrescribedEquations (EquationID);
   virtual int                giveNumberOfDomainEquations (int, EquationID);
   virtual int                giveNumberOfPrescribedDomainEquations (int, EquationID);

  virtual int      giveNewEquationNumber (int domain, DofIDItem) ;
  virtual int      giveNewPrescribedEquationNumber (int domain, DofIDItem) ;

  /// Returns the Equation scaling flag, which is used to indicate that governing equation(s) are scaled, or non-dimensionalized
  virtual bool giveEquationScalingFlag () {return equationScalingFlag;}
  /// Returns the scale factor for given variable type
  virtual double giveVariableScale (VarScaleType varId); 
 /** 
  Prints output of receiver to ouput domain stream, for given time step.
  Corresponding function for element gauss points is invoked
  (gaussPoint::printOutputAt).
  */
  //virtual void                  printOutputAt (FILE *, TimeStep*) ;

  virtual int       requiresUnknowsDictionaryUpdate () {return renumberFlag;}
  virtual bool      requiresEquationRenumbering(TimeStep*) {return renumberFlag;}
  virtual void      updateDofUnknownsDictionary (DofManager*, TimeStep*);
  /*
    Here we store only total and inceremental value; so hash is computed from mode value only
  */
  virtual int       giveUnknownDictHashIndx (EquationID type, ValueModeType mode, TimeStep* stepN);
  
 /**
  Forces equation renumbering on given domain. All equation numbers in all dofManagers are invalidated,
  and new equation numbers are generated starting from domainNeqs entry corresponding to given domain. 
  It will update numberOfEquations variable accordingly.
  Should be used at startup to force equation numbering and therefore sets numberOfEquations.
  Must be used if model supports changes of static system to assign  new valid equation numbers
  to dofManagers.
  */
  virtual int       forceEquationNumbering (int i);
 /**
  Forces equation renumbering on all domains associated to engng model. 
  All equation numbers in all domains for all dofManagers are invalidated,
  and new equation numbers are generated starting from 1 on each domain. 
  It will update numberOfEquations variable accordingly.
  Should be used at startup to force equation numbering and therefore sets numberOfEquations.
  Must be used if model supports changes of static system to assign  new valid equation numbers
  to dofManagers.
  */
  virtual int       forceEquationNumbering () {return EngngModel::forceEquationNumbering();}



protected:
  /**
    Updates nodal values
    (calls also this->updateDofUnknownsDictionary for updating dofs unknowns dictionaries
    if model supports changes of static system). The element internal state update is also forced using
    updateInternalState service.
    */
  void updateInternalState (TimeStep *);
  void applyIC (TimeStep*);
  void assembleAlgorithmicPartOfRhs (FloatArray& rhs, EquationID ut, TimeStep* tStep, int nite);
  void evaluateElementStabilizationCoeffs (TimeStep*) ;
  void updateElementsForNewInterfacePosition (TimeStep*) ;

  void updateDofUnknownsDictionary_predictor (TimeStep* tStep);
  void updateDofUnknownsDictionary_corrector (TimeStep* tStep);

  //void initDofManActivityMap ();
  //void updateDofManActivityMap (TimeStep* tStep);
  void updateDofManVals (TimeStep* tStep);
  //void imposeAmbientPressureInOuterNodes(SparseMtrx* lhs, FloatArray* rhs, TimeStep* stepN);
  //void    __debug  (TimeStep* atTime);

} ;

#define supg_h
#endif
