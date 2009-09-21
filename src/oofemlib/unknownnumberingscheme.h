/* $Id$ */
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

#ifndef unknownnumberingscheme_h
#define unknownnumberingscheme_h

#include "dof.h"

/**
   Abstract base class allowing to controll the way, how equations are assigned to individual DOFs.
   The instances are typically used in EngngModel to asseble characteristic contributions and they
   allow to control the numbering of unknowns.
 */
class UnknownNumberingScheme {
 protected:
 public:
  UnknownNumberingScheme (void) {};
  virtual ~UnknownNumberingScheme () {}
  
  /**
     Initializes the receiver, if necessary
  */
  virtual void init () = 0;
  /**
     Returns true, if receiver is the default engngModel equation numbering scheme;
     This is useful for some components (typically elements), that cache their code numbers
     for default numbering to avoid repeated evaluation.
  */
  virtual bool isDefault() const {return false;}
  /**
     Returns the equation number for corresponding DOF. The numbering should return nonzero value if
     the equation is assigned to the given DOF, zero otherwise.
  */
  virtual int giveDofEquationNumber (Dof* dof) const = 0;
  
  /**
    Returns reqired number of domain equation. Number is always less or equal to the sum of all DOFs gathered from all nodes.
  */
  virtual int giveRequiredNumberOfDomainEquation () const {return 0;}
};

/**
   The representation of EngngModel default unknown numbering. The equation numbers are assigned 
   by the engng model itself to individual DOFs. Therefore, this call is a simple shell around
   DofEquationNumbering interface, forwarding all the reqests to individual DOFs.
 */
class EModelDefaultEquationNumbering : public UnknownNumberingScheme {
 protected:
 public:
  
  EModelDefaultEquationNumbering (void) : UnknownNumberingScheme () {}
  
  virtual void init () {}
  virtual bool isDefault() const {return true;}
  virtual int giveDofEquationNumber (Dof* dof) const {
    return dof->__giveEquationNumber();
  }
};


/**
   The representation of EngngModel default prescribed unknown numbering. 
   The equation numbers are assigned by the engng model itself to individual DOFs. 
   Therefore, this call is a simple shell around
   DofEquationNumbering interface, forwarding all the reqests to individual DOFs.
 */
class EModelDefaultPrescribedEquationNumbering : public UnknownNumberingScheme {
 public:
  
  EModelDefaultPrescribedEquationNumbering (void) : UnknownNumberingScheme () {}
  
  virtual void init () {}
  virtual int giveDofEquationNumber (Dof* dof) const {
    return dof->__givePrescribedEquationNumber();
  }
};
 

#endif // unknownnumberingscheme_h
