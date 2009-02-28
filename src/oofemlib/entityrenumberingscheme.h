/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2008   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

//
// FILE: entityrenumberingscheme.h
//

#ifndef entityrenumberingscheme_h
#define entityrenumberingscheme_h

/**
 * Type allowing to specify the required renumbering scheme;
 * One can have a renumbering scheme for dof managers
 * and another one for elements;
 */
enum EntityRenumberingScheme {
    ERS_DofManager,
    ERS_Element
};


class EntityRenumberingFunctor {
 public:
  // two possible functions to call member function. virtual cause derived
  // classes will use a pointer to an object and a pointer to a member function
  // to make the function call
  virtual int operator()(int, EntityRenumberingScheme)=0;  // call using operator
  virtual int call(int, EntityRenumberingScheme)=0;        // call using function
};

// derived template class
template <class TClass> class SpecificEntityRenumberingFunctor : public EntityRenumberingFunctor
{
 private:
  int (TClass::*fpt)(int, EntityRenumberingScheme);   // pointer to member function
  TClass* pt2Object;                                   // pointer to object
  
 public:
  
  // constructor - takes pointer to an object and pointer to a member and stores
  // them in two private variables
  SpecificEntityRenumberingFunctor(TClass* _pt2Object, int (TClass::*_fpt)(int, EntityRenumberingScheme))
    { pt2Object = _pt2Object;  fpt=_fpt; };
  
  // override operator "()"
  virtual int operator()(int n, EntityRenumberingScheme ers)
  { return (*pt2Object.*fpt)(n, ers);};              // execute member function
  
  // override function "Call"
  virtual int call(int n, EntityRenumberingScheme ers)
  { return (*pt2Object.*fpt)(n, ers);};             // execute member function
};


#endif // entityrenumberingscheme_h

