/* $Header: /home/cvs/bp/oofem/oofemlib/src/rowcol.h,v 1.2 2003/04/06 14:08:25 bp Exp $ */
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

//   ************************
//   *** CLASS ROW-COLUMN ***
//   ************************


#ifndef rowcol_h
#define rowcol_h

#include "debug.h"

class FloatArray ; class IntArray ;

/**
 This class implements a segment of a unsymmetric matrix stored in
 segmented form (skyline).
 A row-column segment i contains the following items :
 - the i-th row of the lower half of the matrix ('row')
 - the i-th column of the upper part of the matrix ('column')
 - the i-th diagonal coefficient of the matrix ('diag').
 Since the profile of the matrix is supposed to be symmetric (but not
 its coeefficients), the row and the column have the same length.
 The row stores coefficients from left to right, the column from up to
 down ; both exclude the diagonal coefficient. Both start at the first
 non-zero coefficient, whose position is 'start'.
*/
class RowColumn
/*
   This class implements a segment of a unsymmetric matrix stored in
   segmented form (skyline).
 DESCRIPTION :
   A row-column segment i contains the following items :
     - the i-th row of the lower half of the matrix ('row')
     - the i-th column of the upper part of the matrix ('column')
     - the i-th diagonal coefficient of the matrix ('diag').
   Since the profile of the matrix is supposed to be symmetric (but not
   its coeefficients), the row and the column have the same length.
   The row stores coefficients from left to right, the column from up to
   down ; both exclude the diagonal coefficient. Both start at the first
   non-zero coefficient, whose position is 'start'.
 TASKS :
   - storing and returning coefficents. For segment k :
     .any coefficient A(i,k) of the upper part of the matrix (i<k) is acces-
      sed through method 'atU(i)'. It belongs to the column of the segment ;
     .any coefficient A(k,i) of the lower part of the matrix (i<k) is acces-
      sed through method 'atL(i)'. It belongs to the row of the segment ;
     .the coefficient A(k,k) is obtained through method 'atDiag()' ;
   - enlarging itself in order to accomodate more coefficients (method
     'growTo') ;
   - resetting to zero all of its coefficients (method 'reinitialized').
*/
{
protected:

 int        number ;
 int        start ;
 double*    row ;
 double*    column ;
 double     diag ;

public:
 RowColumn (int,int) ;
 ~RowColumn () ;

#     ifdef DEBUG
  double&  atU (int) ;
  double&  atL (int) ;
#     else
  double&  atU (int i)                    { return column[i-start] ;}
  double&  atL (int i)                    { return row[i-start] ;}
#     endif
 double&     atDiag ()                      { return diag ;}
 void        checkBounds (int) ;
 void        checkSizeTowards (const IntArray&) ;
 double      dot (const FloatArray&,char,int,int) ;
 int         giveStart ()                   { return start ;}
 void        growTo (int) ;
 void        zero() ;
 void        printYourself () ;
 int         giveSize()                     { return 1+2*(number-start); }

 RowColumn*  GiveCopy () ;

protected:
 RowColumn (int, int, double*, double*, double);

} ;

#endif // rowcol_h

