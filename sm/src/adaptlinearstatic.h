/* $Header: /home/cvs/bp/oofem/sm/src/adaptlinearstatic.h,v 1.5 2003/04/06 14:08:30 bp Exp $ */
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
// Class AdaptiveLinearStatic
//

#ifndef adaptivelinearstatic_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "linearstatic.h"

class AdaptiveLinearStatic : public LinearStatic
{ 
/*
   This class implements Adaptive LinearStatic Engineering problem.
  Multiple loading cases are not supported.
  Due to linearity of a problem, the complete reanalysis from the beginning 
  is done after adaptive remeshing.
  Solution Steps represent a seriaes of adaptive analyses.
*/

 protected:

  /// Type determining used mesh package
  enum MeshPackageType {MPT_T3D, MPT_TARGE2, MPT_FREEM};

  ErrorEstimator* ee;
  MeshPackageType meshPackage;

 public:
  AdaptiveLinearStatic (int i, EngngModel* _master = NULL) : LinearStatic (i,_master) {ee=NULL;}
  ~AdaptiveLinearStatic () {}
// solving
  void solveYourselfAt (TimeStep *);
  /**
  Initializes the newly generated discretization state acording to previous solution.
  This process should typically include restoring old solution, instanciating newly
  generated domain(s) and by mapping procedure. 
  */
  virtual int                initializeAdaptive (int stepNumber);
 /**
  Restores the  state of model from output stream. Restores not only the receiver state,
  but also same function is invoked for all DofManagers and Elements in associated
  domain. Note that by restoring element  context also contexts of all associated
  integration points (and material statuses) are restored.
  */  
  virtual contextIOResultType                restoreContext (FILE* stream, void* obj = NULL) ;

  void updateDomainLinks();
 
  IRResultType initializeFrom (InputRecord* ir);
  /** Service for accessing ErrorEstimator corresponding to particular domain */
  ErrorEstimator* giveDomainErrorEstimator (int n) {return ee;}

// identification
   const char* giveClassName () const { return "AdaptiveLinearStatic";}
  classType giveClassID ()      const { return AdaptiveLinearStaticClass;}
} ;

#define adaptivelinearstatic_h
#endif
