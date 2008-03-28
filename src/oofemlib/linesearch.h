/* $Header: /home/cvs/bp/oofem/oofemlib/src/linesearch.h,v 1.3 2003/04/06 14:08:25 bp Exp $ */
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


//   **************************
//   *** CLASS LineSearchNM ***
//   **************************

 
#ifndef linesearch_h
#define linesearch_h


#include "nummet.h"
#include "cltypes.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

class EngngModel ; class SparseMtrx; class FloatArray; 


/**
This base class is an abstraction/implementation for numerical method solving 
line search optimisation problem.
*/
class LineSearchNM : public NumericalMethod
{
public:
 enum LS_status {ls_ok, ls_failed};

 int max_iter;
 double ls_tolerance;
 double amplifFactor;
 double maxEta, minEta;
 FloatArray eta;
 FloatArray prod;

protected:
public:
 /// Constructor
 LineSearchNM (int i, Domain* d,EngngModel* m);
 /// Destructor
 ~LineSearchNM ();

 // identification 
 /// Returns class name of the receiver.
 const char*  giveClassName () const { return "LineSearchNM" ;}
 /** Returns classType id of receiver.
  @see FEMComponent::giveClassID 
  */
 classType giveClassID () const { return NumericalMethodClass ;}

 /**
  Solves the line search optimalisation problem in the form of $g(r)=0; r_{new}=r_{old}+\eta\delta r; 0 < \eta < 1$, 
  The aim is to find $\eta$ so that the $g(r)$ has decresed sufficiently.
  The total solution vector is updated at exit as well as InternalRhs vector.
  @param r  old total solution (total displacement)
  @param dr increment of solution (incremental displacaments)
  @param F  old InternalRhs (real internal forces)
  @param R  reference incremental Rhs (incremental load)
  @param R0 initial Rhs (initial load)
  @param etaValue reached eta value
  @param status linesearch status
  */
 virtual NM_Status solve (FloatArray* r, FloatArray* dr, FloatArray* F, FloatArray* R, FloatArray* R0,
              IntArray& eqnmask, double lambda, double&etaValue, LS_status& status, TimeStep*);

 // management  components
 IRResultType initializeFrom (InputRecord* ir) ;
 contextIOResultType    saveContext (DataStream* stream, ContextMode mode, void *obj = NULL) {return CIO_OK;}
 contextIOResultType    restoreContext(DataStream* stream, ContextMode mode, void *obj = NULL) {return CIO_OK;}



protected:
 void search (int istep, FloatArray& prod, FloatArray& eta, double amp, double maxeta, double mineta, int& status);
 };

#endif // linesearch_h









