/* $Header: /home/cvs/bp/oofem/tm/src/nonstationarytransportproblem.h,v 1.2 2003/05/19 13:04:10 bp Exp $ */
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

//
// Class NonStartionaryTransportProblem
//

#ifndef nonstationarytransportproblem_h
#define nonstationarytransportproblem_h

#include "structengngmodel.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

/**
 This class represents linear nonstationary transport problem.
*/
class NonStationaryTransportProblem : public EngngModel
{ 
protected:
  SparseMtrx* lhs;
  FloatArray rhs, bcRhs;
  //FloatArray solutionVector, previousSolutionVector;
 PrimaryField FluxField;
 
 LinSystSolverType solverType;
 SparseMtrxType sparseMtrxType;
 /// Numerical method used to solve the problem
 SparseLinearSystemNM *nMethod;

 double deltaT;
 double alpha;

 // if set then stabilization using lumped capacity will be used
 int lumpedCapacityStab;

 /// if set, the receiver flux field will be exported using FieldManager
 int exportFieldFlag;

 int initFlag;
 /// Associated time function for time step increment
 int dtTimeFunction ;

public:
  NonStationaryTransportProblem (int i, EngngModel* _master = NULL) : EngngModel (i,_master), 
    rhs(), bcRhs(), FluxField (this, 1, FBID_FluxField, EID_ConservationEquation, 1)
    {lhs = NULL; ndomains = 1; nMethod = NULL; lumpedCapacityStab = 0; exportFieldFlag = 0; initFlag=1; dtTimeFunction=0;}
  ~NonStationaryTransportProblem () 
    {if (lhs) delete  lhs; if (nMethod) delete nMethod;}

  void solveYourselfAt (TimeStep *);
 /**
  Updates nodal values
  (calls also this->updateDofUnknownsDictionary for updating dofs unknowns dictionaries
  if model supports changes of static system). The element internal state update is also forced using
  updateInternalState service.
  */
  virtual void               updateYourself (TimeStep*) ;
  double   giveUnknownComponent ( EquationID, ValueModeType, TimeStep*, Domain*, Dof*);
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
  const char* giveClassName () const { return "NonStationaryTransportProblem";}
  classType giveClassID ()      const { return NonStationaryTransportProblemClass;}
  fMode giveFormulation () { return TL; }

 virtual void giveElementCharacteristicMatrix (FloatMatrix& answer, int num, 
                        CharType type, TimeStep* tStep, Domain* domain) ;

  /** DOF printing routine. Called by DofManagers to print Dof specific part.
  Dof class provides component printing routines, but emodel is responsible
  for what will be printed at DOF level.
  @param stream output stream
  @param iDof dof to be processed
  @param atTime solution step
  */
 virtual void printDofOutputAt (FILE* stream, Dof* iDof, TimeStep* atTime);

 /**
    Returns time function for time step increment.
    Used time function should provide step lengths as function of step number. 
    Initial step with number 0 is considered as [ -dt(0), 0 ], first step is [ 0, dt(1) ], ...
 */
 LoadTimeFunction*  giveDtTimeFunction () ;
 
 /**
    Returns the timestep length for given step number n, initial step is number 0
 */
 double giveDeltaT (int n);
 

#ifdef __PETSC_MODULE
   /**
      Creates Petsc contexts. Must be implemented by derived classes since the governing equation type is reqired 
      for context creation.
    */
  virtual void initPetscContexts ();
#endif



protected:
  virtual void assembleAlgorithmicPartOfRhs (FloatArray& rhs, EquationID ut, TimeStep* tStep);
  virtual void applyIC (TimeStep*);
  /**
    Assembles part of rhs due to Dirichlet boundary conditions.
    @param answer global vector where the contribution will be added
    @param tStep solution step
    @param mode CharTypeMode of result
    @param lhsType type of element matrix to be multiplied by vector of prescribed. 
    The giveElementCharacteristicMatrix service is used to get/compute element matrix.
    @param d domain
    */
  virtual void assembleDirichletBcRhsVector (FloatArray& answer, TimeStep* tStep, EquationID ut, ValueModeType mode, 
                                             CharType lhsType, Domain* d);
} ;

#endif // nonstationarytransportproblem_h
