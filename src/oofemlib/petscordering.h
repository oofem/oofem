/* $Header: /home/cvs/bp/oofem/oofemlib/src/Attic/petscordering.h,v 1.1.2.1 2004/04/05 15:19:43 bp Exp $ */
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

#ifndef petscordering_h
#define petscordering_h
#ifdef __PARALLEL_MODE

#include "appordering.h"
#include "intarray.h"
#include "dofmanager.h"
#ifndef __MAKEDEPEND
#include <map>
#endif

class PetscOrdering_Base : public ApplicationOrdering {
protected:
public:
  PetscOrdering_Base() : ApplicationOrdering() {}  
  
  /// Returns true if given DofManager is local (ie maintained by the receiver processor)
  bool isLocal (DofManager*);
  /// Returns true if given DofManager is shared between partitions
  bool isShared (DofManager*);
};

/**
Ordering from oofem natural ordering (includes all local and shared eqs)
to global ordering
*/
class PetscNatural2GlobalOrdering : public PetscOrdering_Base {
protected:
  IntArray locGlobMap; // old->new; uses 0-baseg global eq ordering; 1-based local ordering
  std::map<int,int> globLocMap; // new->old

  // number of local and global eqs.
  int l_neqs, g_neqs;

public:
  PetscNatural2GlobalOrdering();
  ~PetscNatural2GlobalOrdering() {}

  void init (EngngModel*, EquationID ut, int di, EquationType et = et_standard);

  virtual int giveNewEq (int leq);
  virtual int giveOldEq (int eq);

  virtual void map2New (IntArray& answer, const IntArray& src, int baseOffset = 0) ;
  virtual void map2Old (IntArray& answer, const IntArray& src, int baseOffset = 0) ;

  int giveNumberOfLocalEqs() {return l_neqs;}
  int giveNumberOfGlobalEqs() {return g_neqs;}

  IntArray* giveN2Gmap() {return &locGlobMap;}
};

/**
Ordering from oofem natural ordering (includes all local and shared eqs)
to local ordering, where only locally maintained eqs are considered.
*/

class PetscNatural2LocalOrdering : public PetscOrdering_Base {
protected:
  IntArray n2l; // natural->local
public:
  PetscNatural2LocalOrdering();
  ~PetscNatural2LocalOrdering() {}

  void init (EngngModel*, EquationID ut, int di, EquationType et = et_standard);

  virtual int giveNewEq (int leq);
  virtual int giveOldEq (int eq);

  virtual void map2New (IntArray& answer, const IntArray& src, int baseOffset = 0) ;
  virtual void map2Old (IntArray& answer, const IntArray& src, int baseOffset = 0) ;

  IntArray* giveN2Lmap() {return &n2l;}

};

#endif
#endif
