/* $Header: /home/cvs/bp/oofem/sm/src/linearstatic.h,v 1.7.4.1 2004/04/05 15:19:47 bp Exp $ */
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
// Class LinearStatic
//

#ifndef linearstatic_h
#define linearstatic_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "structengngmodel.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"

class LinearStatic : public StructuralEngngModel
{ 
/*
   This class implements LinearStatic Engineering problem.
 Multiple loading works only if linear elastic material (such as isoLE)  is used.
 (Other non-linear materials kepp load history, so such multiple loading
 will cause that next step will be assumed as new load increment, 
 not the total new load). Because they always copute real stresses acording
 to reached strain state, they are not able to respond to linear analysis.
  
DESCRIPTION:
   Solution of this problem is series of loading cases, maintained as sequence of
   time-steps. This solution is in form of linear equation system Ax=b
TASK:
   Creating Numerical method for solving Ax=b
   Interfacing Numerical method to Elements
   Managing time  steps
*/

 protected:
  SparseMtrx* stiffnessMatrix;
  FloatArray loadVector;
  FloatArray displacementVector;
  
  LinSystSolverType solverType;
  SparseMtrxType sparseMtrxType;
  /// Numerical method used to solve the problem
  SparseLinearSystemNM *nMethod;

  int initFlag;

#ifdef __PETSC_MODULE
  Vec _loadVec, _dispVec;
#endif

 public:
  LinearStatic (int i, EngngModel* _master = NULL) ;
  ~LinearStatic () ;
// solving
  void solveYourself ();
  void solveYourselfAt (TimeStep *);
  //int requiresNewLhs () {return 0;}
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
  NumericalMethod* giveNumericalMethod (TimeStep*);
  //void printReactionForces (TimeStep *);
  void               terminate (TimeStep*);

 IRResultType initializeFrom (InputRecord* ir);

 // consistency check
 virtual int checkConsistency (); // returns nonzero if o.k.

  /** DOF printing routine. Called by DofManagers to print Dof specific part.
  Dof class provides component printing routines, but emodel is responsible
  for what will be printed at DOF level.
  @param stream output stream
  @param iDof dof to be processed
  @param atTime solution step
  */
 virtual void printDofOutputAt (FILE* stream, Dof* iDof, TimeStep* atTime);

// identification
  const char* giveClassName () const { return "LinearStatic";}
  classType giveClassID ()      const { return LinearStaticClass;}
  fMode giveFormulation () { return TL; }
#ifdef __PARALLEL_MODE
  /**
     Determines the space necessary for send/receive buffer.
     It uses related communication map pattern to determine the maximum size needed.
     @param commMap communication map used to send/receive messages
     @param buff communication buffer
     @return upper bound of space needed
  */
  int estimateMaxPackSize (IntArray& commMap, CommunicationBuffer& buff, int packUnpackType) ;
#endif
 
} ;

#endif // linearstatic_h
