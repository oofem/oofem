/* $Header: /home/cvs/bp/oofem/sm/src/doublepowerlaw.h,v 1.4 2003/04/06 14:08:30 bp Exp $ */
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

//   *********************************************
//   *** CLASS RHEOLOGIC DOUBLE POWER LAW Material
//   *********************************************
#ifndef doublepowerlaw_h 
#define doublepowerlaw_h 

#include "cltypes.h"
#include "maxwellChM.h"

class DoublePowerLawMaterial : public MaxwellChainMaterial
{
/*
   This class implements a rheologic Maxwelll chain model in a finite
 element problem. 
 
 DESCRIPTION
 TASK
*/
protected:
 double E28;    // Young modulus at age of 28 days [MPa]
 double fi1;    // basic creep coefficient
 double m,n;
 double alpha;

public:
 DoublePowerLawMaterial (int n,Domain* d) : MaxwellChainMaterial (n, d) {}
 ~DoublePowerLawMaterial () {}
 


 const char* giveClassName () const { return "DoublePowerLawMaterial" ;}
 classType giveClassID ()         const      {return DoublePowerLawMaterialClass;}
  IRResultType initializeFrom (InputRecord* ir);
protected:

 virtual double  computeCreepFunction (GaussPoint* gp, double atTime, double ofAge);
} ;


#endif // doublepowerlaw_h 
