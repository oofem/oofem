/* $Header: /home/cvs/bp/oofem/sm/src/incrementallinearstatic.h,v 1.5 2003/05/19 13:03:59 bp Exp $ */
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
// Class Incremental LinearStatic
//

#ifndef incrementallinearstatic_h
#define incrementallinearstatic_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "structengngmodel.h"
#include "sparselinsystemnm.h"
#include "sparselinsystemnm.h"
#include "cltypes.h"

/**
  This class implements Incremental LinearStatic Engineering problem.
 problem is solved as series of linear solutions. This class is intended to 
 be used for solving linear creep problems or incremental perfect plasticity.
 Supports the changes of static scheme (applying, removing and changing  boundary conditions) 
 during the analysis.
*/
class IncrementalLinearStatic : public StructuralEngngModel
{ 
/*
  This class implements Incremental LinearStatic Engineering problem.
 problem is solved as series of linear solutions. This class is intended to 
 be used for solving linear creep problems or incremental perfect plasticity.
 Of course approprite material models must be used.
 
  
DESCRIPTION:
   Solution of this problem is series of loading steps, maintained as sequence of
   time-steps. This solution is in form of linear equation system Ax=b
TASK:
   Creating Numerical method for solving Ax=b
   Interfacing Numerical method to Elements
   Managing time  steps
*/

 protected:
  SparseMtrx* stiffnessMatrix;
  FloatArray incrementOfLoadVector;
  FloatArray incrementOfDisplacementVector;
  FloatArray totalDisplacementVector;
  FloatArray discreteTimes;
  
  double endOfTimeOfInterest;
 /// Numerical method used to solve the problem
 SparseLinearSystemNM *nMethod;

 LinSystSolverType solverType;
 SparseMtrxType sparseMtrxType;
  

 public:
  IncrementalLinearStatic (int i, EngngModel* _master = NULL);
  ~IncrementalLinearStatic () ;
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
  virtual void               updateYourself (TimeStep *) ;
  //double   giveUnknownComponent (CharType, int);
  contextIOResultType saveContext (DataStream* stream, ContextMode mode, void *obj = NULL) ;
  contextIOResultType restoreContext (DataStream* stream, ContextMode mode, void *obj = NULL);
  TimeStep* giveNextStep ();

  double             giveDiscreteTime (int);
  IRResultType initializeFrom (InputRecord* ir);
  virtual int        giveNumberOfFirstStep () {return 1;}
  virtual double     giveEndOfTimeOfInterest () {return endOfTimeOfInterest;}

  NumericalMethod* giveNumericalMethod (TimeStep*);
  /** DOF printing routine. Called by DofManagers to print Dof specific part.
  Dof class provides component printing routines, but emodel is responsible
  for what will be printed at DOF level.
  @param stream output stream
  @param iDof dof to be processed
  @param atTime solution step
  */
 virtual void printDofOutputAt (FILE* stream, Dof* iDof, TimeStep* atTime);
  

// identification
  const char* giveClassName () const { return "IncrementalLinearStatic";}
  classType giveClassID () const {return IncrementalLinearStaticClass;}
  fMode giveFormulation () { return TL; } 
  // virtual  LoadResponseMode giveLoadResponseMode () {return IncrementOfLoad;}
 virtual int       requiresUnknownsDictionaryUpdate () {return 1;}
 virtual bool      requiresEquationRenumbering(TimeStep*) {return true;}
 virtual void      updateDofUnknownsDictionary (DofManager*, TimeStep*);
 /*
   Here we store only total and inceremental value; so hash is computed from mode value only
 */
 virtual int       giveUnknownDictHashIndx (EquationID type, ValueModeType mode, TimeStep* stepN)
   {return (int) mode;}


  
} ;

#endif // incrementallinearstatic_h
