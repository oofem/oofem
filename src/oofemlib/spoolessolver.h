/* $Header: /home/cvs/bp/oofem/oofemlib/src/spoolessolver.h,v 1.1 2003/04/06 14:08:26 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2003   Borek Patzak                                       



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
#ifndef spoolessolver_h
#define spoolessolver_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "flotarry.h"
#include "cltypes.h"
#include "spoolesinterface.h"

class Domain; class EngngModel; class FloatMatrix;

/**
 Implements the solution of linear system of equation in the form Ax=b using solvers 
 from SPOOLES library. Can work with only spooles sparse matrix implementation. 
*/
class SpoolesSolver : public SparseLinearSystemNM
{
private:
#ifdef __SPOOLES_MODULE
  /// last mapped Lhs matrix
  SparseMtrx*    Lhs;
  /// last mapped matrix version
 SparseMtrx::SparseMtrxVersionType lhsVersion;
  int            msglvl;
  FILE           *msgFile;
  int            msgFileCloseFlag;

  FrontMtx        *frontmtx ;
  IV              *oldToNewIV, *newToOldIV;
  ETree           *frontETree ;
 IVL             *adjIVL, *symbfacIVL ;
 SubMtxManager   *mtxmanager  ; 
  Graph           *graph ;
#endif
public:
 
 /**
  Constructor - creates new instance of LDLTFactorization, with number i, belonging to domain d and Engngmodel m.
  */
 SpoolesSolver (int i, Domain* d,EngngModel* m);

  ///Destructor
 ~SpoolesSolver () ;// destructor

  /**
    Solves the given linear system by LDL^T factorization. 
    Implementation rely on factorization support provided by mapped sparse matrix.
    It calls Lhs->factorized()->backSubstitutionWith(*solutionArray). Sets solved flag to 1 if o.k.
    @param A coefficient matrix 
    @param b right hand side
    @param x solution array
    @return NM_Status value
    @param tNow time step
    */
  NM_Status solve (SparseMtrx* A, FloatArray* b, FloatArray* x);

  /// Initializes receiver from given record. Empty implementation.
 IRResultType initializeFrom (InputRecord* ir);

  // identification 
 const char*  giveClassName () const { return "SpoolesSolver" ;}
  classType giveClassID () const { return  SpoolesSolverClass ;}
  LinSystSolverType giveLinSystSolverType() const {return ST_Spooles;}
 };


#endif // spoolessolver_h









