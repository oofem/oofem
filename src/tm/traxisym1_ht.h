/* $Header: /home/cvs/bp/oofem/tm/src/traxisym1_ht.h,v 1.2 2003/04/23 14:22:15 bp Exp $ */
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

//   *************************************************************************************************
//   *** CLASS TrAxisym1_ht: Axisymmtric Triangle(2d), linear approximation, Heat Transfer element ***
//   *************************************************************************************************
 
#ifndef traxisym1_ht_h
#define traxisym1_ht_h

#include "tr1_ht.h"

class TrAxisym1_ht : public Tr1_ht
{
 
protected:

public:
 
 // constructor
 TrAxisym1_ht (int,Domain*, ElementMode em = HeatTransferEM) ;
 ~TrAxisym1_ht () ;                        // destructor

 double                computeVolumeAround (GaussPoint*);
 const char* giveClassName () const { return "TrAxisym1_htElement" ;}
 classType                giveClassID () const { return TrAxisym1_htClass; }

protected:
 double computeEdgeVolumeAround(GaussPoint* gp, int iEdge) ;
 double computeRadiusAt (GaussPoint* gp);
} ;

#endif // traxisym1_ht_h
