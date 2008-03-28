/* $Header: /home/cvs/bp/oofem/oofemlib/src/range.h,v 1.3 2003/04/06 14:08:25 bp Exp $ */
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

#ifndef range_h
#define range_h

/**
 Class Range is an abstraction for interval of integer numbers. It is decribed using its start and end values of interval 
 it represents. The inteval is defined to represent all valus between start and end values, including start and end values.
 Function for testing if number is in interval is provided. Used by OutputManager
 to eficiently maintain intervals.
 */
class Range 
{
protected:
 /// interval start value
 int startIndx;
 /// interval end value
 int endIndx;

public:
 /// Constructor. Creates Range containing only given single number
 Range (int indx) {startIndx = endIndx = indx;}
 /// Consructor. Creates range <li, hi>
 Range (int li, int hi) {startIndx = li; endIndx = hi;}
 /// Empty range constructor
 Range () {startIndx = 0; endIndx = -1;}


 /// tests if number is in range
 int test (int i) {return ((i>=startIndx) && (i<=endIndx))?1:0;}
};

#endif // range_h
