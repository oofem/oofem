/* $Header: /home/cvs/bp/oofem/sm/src/eigenvaluedynamic.h,v 1.7 2003/04/06 14:08:30 bp Exp $ */
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
// Class EigenValueDynamic
//

#ifndef eigenvaluedynamic_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "engngm.h"
#include "sparsegeneigenvalsystemnm.h"
#include "sparsemtrx.h"
#include "flotmtrx.h"
#include "flotarry.h"

class EigenValueDynamic : public EngngModel
{ 
/*
   This class implements way for examining eigen values and eigen vectors in 
   dynamic problems.
DESCRIPTION:
   Solution of this problem is base on equation in the form of: Ky=wMy
   Currently eigen value problem is solved using subspace iteration.
TASK:
   Assembling the governing equation in the form  Ky=wMy
   Creating Numerical method for Ky=wMy
   Interfacing Numerical method to Elements
*/

 private:
  SparseMtrx* stiffnessMatrix;
  SparseMtrx* massMatrix;
  FloatMatrix eigVec;
  FloatArray  eigVal;
  int      numberOfRequiredEigenValues ;
  int      activeVector ; 
 int restoreFlag;
  double   rtolv;           // precision
 /// Numerical method used to solve the problem
 SparseGeneralEigenValueSystemNM *nMethod;
 GenEigvalSolverType solverType;
 

 public:
  EigenValueDynamic (int i, EngngModel* _master = NULL) : EngngModel (i,_master) 
    {stiffnessMatrix = NULL; massMatrix = NULL; 
   numberOfSteps = 1;  ndomains = 1; nMethod = NULL;}
  ~EigenValueDynamic () {delete  stiffnessMatrix; delete massMatrix; 
             if (nMethod) delete nMethod;}
// solving
  //void solveYourself ();
  void solveYourselfAt (TimeStep *);
  void terminate (TimeStep *);
  int requiresNewLsh () {return 0;}
  virtual void               updateYourself (TimeStep *) ;
 /*
  The active eigen value is identified by inrinsic time of TimeStep.
  */
  double   giveUnknownComponent ( EquationID, ValueModeType, TimeStep*, Domain*, Dof*);
  double   giveUnknownComponent ( UnknownType, ValueModeType, TimeStep*, Domain*, Dof*);
  IRResultType initializeFrom (InputRecord* ir);
  contextIOResultType saveContext (FILE* stream, void *obj = NULL) ;
  contextIOResultType restoreContext (FILE* stream, void *obj = NULL);
  TimeStep* giveNextStep ();
  NumericalMethod* giveNumericalMethod (TimeStep*);
  void   setActiveVector (int i) {activeVector = i;}
 int resolveCorrespondingEigenStepNumber (void* obj);

  /** DOF printing routine. Called by DofManagers to print Dof specific part.
  Dof class provides component printing routines, but emodel is responsible
  for what will be printed at DOF level.
  @param stream output stream
  @param iDof dof to be processed
  @param atTime solution step
  */
 virtual void printDofOutputAt (FILE* stream, Dof* iDof, TimeStep* atTime);

// identification
  const char* giveClassName () const { return "EigenValueDynamic";}
  classType giveClassID ()      const { return EigenValueDynamicClass;}
} ;

#define linearstatic_h
#endif
