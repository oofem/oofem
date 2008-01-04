/* $Header: /home/cvs/bp/oofem/oofemlib/src/matconst.h,v 1.4 2003/04/06 14:08:25 bp Exp $ */
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


/*
This file matconst.h defines several material constant (respective their
representative number).

These constant are used to acces material parameters  at material level from 
property dictionary 

*/

/* poiss coeficients */
#define NYxz 300
#define NYyz 301
#define NYxy 302
#define NYzx 303
#define NYzy 304
#define NYyx 305


/* Ortho Young moduluses */
#define Ex 400
#define Ey 401
#define Ez 402

/* thermal dilatation coeffs */
#define tAlphax 403
#define tAlphay 404
#define tAlphaz 405
#define tAlpha 406 // 1d.

/* shear modulus */
#define Gyz 407
#define Gxz 408
#define Gxy 409

/* heat capacity */
#define HeatCapaCoeff 450
#define Mass1CapaCoeff 451

/* Viscosity */
#define Viscosity 500
#define YieldStress 501
