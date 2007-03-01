/* $Header: /home/cvs/bp/oofem/tm/src/nltransienttransportproblem.h,v 1.1 2003/04/14 16:01:39 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2005   Borek Patzak                                       



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
// Class NonLinearTransientTransportProblem
//

#ifndef cbs_h

#include "engngm.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"
//<RESTRICTED_SECTION>
#include "materialinterface.h"
//</RESTRICTED_SECTION>
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

/**
 This class represents CBS algorithm for solving incompressible Navier-Stokes equations
*/
class CBS : public EngngModel
{ 
protected:
  /// Numerical method used to solve the problem
  SparseLinearSystemNM *nMethod;

  LinSystSolverType solverType;
  SparseMtrxType sparseMtrxType;
  
  SparseMtrx* lhs;
  /// Pressure field
  PrimaryField PressureField;
  /// Velocity field
  PrimaryField VelocityField;
  FloatArray deltaAuxVelocity;
  FloatArray prescribedTractionPressure;
  FloatArray nodalPrescribedTractionPressureConnectivity;
  
  /** lumped mass matrix */
  FloatArray mm;
  /** Sparse consistent mass */
  SparseMtrx* mss;

  /// time step and its minimal value
  double deltaT, minDeltaT;
  // integration constants
  double theta[2];

  int initFlag;
  /// consistent mass flag
  int consistentMassFlag;

  int numberOfMomentumEqs, numberOfConservationEqs;
  int numberOfPrescribedMomentumEqs, numberOfPrescribedConservationEqs;

  bool equationScalingFlag;
  /// length, velocity, and density scales
  double lscale, uscale, dscale;
  /// Reynolds number
  double Re;

//<RESTRICTED_SECTION>
  // material interface representation for multicomponent flows
  MaterialInterface* materialInterface;
//</RESTRICTED_SECTION>
public:
  CBS (int i, EngngModel* _master = NULL) : EngngModel (i,_master), PressureField(this,1,FBID_PressureField, EID_ConservationEquation, 1),
    VelocityField (this,1,FBID_VelocityField, EID_MomentumBalance, 1) {initFlag=1; lhs=NULL; ndomains=1;nMethod=NULL;
    numberOfMomentumEqs=numberOfConservationEqs=numberOfPrescribedMomentumEqs=numberOfPrescribedConservationEqs=0;consistentMassFlag = 0;
    equationScalingFlag = false; lscale=uscale=dscale=1.0;
//<RESTRICTED_SECTION>
    materialInterface = NULL;
//</RESTRICTED_SECTION>
  }
  ~CBS () {
//<RESTRICTED_SECTION>
    if (materialInterface) delete materialInterface;
//</RESTRICTED_SECTION>
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
  contextIOResultType saveContext (FILE* stream, void *obj = NULL) ;
  contextIOResultType restoreContext (FILE* stream, void *obj = NULL);
 
  void   updateDomainLinks();
  
  TimeStep* giveNextStep ();
  TimeStep*  giveSolutionStepWhenIcApply();
  NumericalMethod* giveNumericalMethod (TimeStep*);
  
  /// Initialization from given input record
  IRResultType initializeFrom (InputRecord* ir);

  // consistency check
  virtual int checkConsistency (); // returns nonzero if o.k.
  // identification
  const char* giveClassName () const { return "CBS";}
  classType giveClassID ()      const { return CBSClass;}
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
} ;

#define cbs_h
#endif
