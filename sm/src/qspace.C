/* $Header: /home/cvs/bp/oofem/sm/src/qspace.C,v 1.3.4.1 2004/04/05 15:19:47 bp Exp $ */
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

//   file QSPACE.CC

#include "qspace.h"
#include "node.h"
#include "material.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

QSpace :: QSpace (int n, Domain* aDomain)
      : StructuralElement (n,aDomain)
   // Constructor.
{
   numberOfDofMans  = 20 ;
   numberOfGaussPoints = 27;
   // this -> computeGaussPoints() ; => moved to instanciateYourself();
}


void
QSpace :: computeBmatrixAt (GaussPoint *aGaussPoint, FloatMatrix& answer, int li, int ui)
   // Returns the [6x60] strain-displacement matrix {B} of the receiver, eva-
   // luated at aGaussPoint.
{
   FloatMatrix jacMtrx, inv; // ,*answer ;
   FloatArray  *coord ;
   double   u,v,w,
            x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,
            y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,
            z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16,z17,z18,z19,z20 ;
   double j11,j12,j13,j21,j22,j23,j31,j32,j33 ;

   // natural coordinates ksi, eta and dzeta :

   coord = aGaussPoint -> giveCoordinates() ;
   u = coord -> at(1) ;
   v = coord -> at(2) ;
   w = coord -> at(3) ;
 
   // partial derivatives of ksi and eta, with respect to x and y :

   this -> computeJacobianMatrixAt(jacMtrx, coord) ;
   inv.beInverseOf (jacMtrx);

   j11 = inv . at(1,1) ;
   j12 = inv . at(1,2) ;
   j13 = inv . at(1,3) ;
   j21 = inv . at(2,1) ;
   j22 = inv . at(2,2) ;
   j23 = inv . at(2,3) ;
   j31 = inv . at(3,1) ;
   j32 = inv . at(3,2) ;
   j33 = inv . at(3,3) ;

   // partial derivatives of the shape functions N_i, with respect to ksi
  
   x1  =  0.125*(1.+v)*(1.+w)*(2.*u+v+w-1.);
   x2  = -0.5*u*(1.+v)*(1.+w);
   x3  = -0.125*(1.+v)*(1.+w)*(-2.*u+v+w-1.);
   x4  = -0.25*(1.-v*v)*(1.+w);
   x5  = -0.125*(1.-v)*(1.+w)*(-2.*u-v+w-1.);
   x6  = -0.5*u*(1.-v)*(1.+w);
   x7  =  0.125*(1.-v)*(1.+w)*(2.*u-v+w-1.);
   x8  =  0.25*(1.-v*v)*(1.+w);
   x9  =  0.25*(1.+v)*(1.-w*w);
   x10 = -0.25*(1.+v)*(1.-w*w);
   x11 = -0.25*(1.-v)*(1.-w*w);
   x12 =  0.25*(1.-v)*(1.-w*w);
   x13 =  0.125*(1.+v)*(1.-w)*(2.*u+v-w-1.);
   x14 = -0.5*u*(1.+v)*(1.-w);
   x15 = -0.125*(1.+v)*(1.-w)*(-2.*u+v-w-1.);
   x16 = -0.25*(1.-v*v)*(1.-w);
   x17 = -0.125*(1.-v)*(1.-w)*(-2.*u-v-w-1.);
   x18 = -0.5*u*(1.-v)*(1.-w);
   x19 =  0.125*(1.-v)*(1.-w)*(2.*u-v-w-1.);
   x20 =  0.25*(1.-v*v)*(1.-w);

   // partial derivatives of the shape functions N_i, with respect to eta

   y1  =  0.125*(1.+u)*(1.+w)*(u+2.*v+w-1.);
   y2  =  0.25*(1.-u*u)*(1.+w);
   y3  =  0.125*(1.-u)*(1.+w)*(-u+2.*v+w-1.);
   y4  = -0.5*v*(1.-u)*(1.+w);
   y5  = -0.125*(1.-u)*(1.+w)*(-u-2.*v+w-1.);
   y6  = -0.25*(1.-u*u)*(1.+w);
   y7  = -0.125*(1.+u)*(1.+w)*(u-2.*v+w-1.);
   y8  = -0.5*v*(1.+u)*(1.+w);
   y9  =  0.25*(1.+u)*(1.-w*w);
   y10 =  0.25*(1.-u)*(1.-w*w);
   y11 = -0.25*(1.-u)*(1.-w*w);
   y12 = -0.25*(1.+u)*(1.-w*w);
   y13 =  0.125*(1.+u)*(1.-w)*(u+2.*v-w-1.);
   y14 =  0.25*(1.-u*u)*(1.-w);
   y15 =  0.125*(1.-u)*(1.-w)*(-u+2.*v-w-1.);
   y16 = -0.5*v*(1.-u)*(1.-w);
   y17 = -0.125*(1.-u)*(1.-w)*(-u-2.*v-w-1.);
   y18 = -0.25*(1.-u*u)*(1.-w);
   y19 = -0.125*(1.+u)*(1.-w)*(u-2.*v-w-1.);
   y20 = -0.5*v*(1.+u)*(1.-w);

   // partial derivatives of the shape functions N_i, with respect to dzeta

   z1  =  0.125*(1.+u)*(1.+v)*(u+v+2.*w-1.);
   z2  =  0.25*(1.-u*u)*(1.+v);
   z3  =  0.125*(1.-u)*(1.+v)*(-u+v+2.*w-1.);
   z4  =  0.25*(1.-u)*(1.-v*v);
   z5  =  0.125*(1.-u)*(1.-v)*(-u-v+2.*w-1.);
   z6  =  0.25*(1.-u*u)*(1.-v);
   z7  =  0.125*(1.+u)*(1.-v)*(u-v+2.*w-1.);
   z8  =  0.25*(1.+u)*(1.-v*v);
   z9  = -0.5*w*(1.+u)*(1.+v);
   z10 = -0.5*w*(1.-u)*(1.+v);
   z11 = -0.5*w*(1.-u)*(1.-v);
   z12 = -0.5*w*(1.+u)*(1.-v);
   z13 = -0.125*(1.+u)*(1.+v)*(u+v-2.*w-1.);
   z14 = -0.25*(1.-u*u)*(1.+v);
   z15 = -0.125*(1.-u)*(1.+v)*(-u+v-2.*w-1.);
   z16 = -0.25*(1.-u)*(1.-v*v);
   z17 = -0.125*(1.-u)*(1.-v)*(-u-v-2.*w-1.);
   z18 = -0.25*(1.-u*u)*(1.-v);
   z19 = -0.125*(1.+u)*(1.-v)*(u-v-2.*w-1.);
   z20 = -0.25*(1.+u)*(1.-v*v);

   // B matrix  -  6 rows : epsilon-X, epsilon-Y, epsilon-Z, gamma-YZ, gamma-ZX, gamma-XY  :

   // answer = new FloatMatrix(6,60) ;
  answer.resize (6, 60);

   answer.at(1,1)  = x1*j11+y1*j12+z1*j13 ;
   answer.at(1,4)  = x2*j11+y2*j12+z2*j13 ;
   answer.at(1,7)  = x3*j11+y3*j12+z3*j13 ;
   answer.at(1,10) = x4*j11+y4*j12+z4*j13 ;
   answer.at(1,13) = x5*j11+y5*j12+z5*j13 ;
   answer.at(1,16) = x6*j11+y6*j12+z6*j13 ;
   answer.at(1,19) = x7*j11+y7*j12+z7*j13 ;
   answer.at(1,22) = x8*j11+y8*j12+z8*j13 ;
   answer.at(1,25) = x9*j11+y9*j12+z9*j13 ;
   answer.at(1,28) = x10*j11+y10*j12+z10*j13 ;
   answer.at(1,31) = x11*j11+y11*j12+z11*j13 ;
   answer.at(1,34) = x12*j11+y12*j12+z12*j13 ;
   answer.at(1,37) = x13*j11+y13*j12+z13*j13 ; 
   answer.at(1,40) = x14*j11+y14*j12+z14*j13 ;
   answer.at(1,43) = x15*j11+y15*j12+z15*j13 ;
   answer.at(1,46) = x16*j11+y16*j12+z16*j13 ; 
   answer.at(1,49) = x17*j11+y17*j12+z17*j13 ;
   answer.at(1,52) = x18*j11+y18*j12+z18*j13 ;
   answer.at(1,55) = x19*j11+y19*j12+z19*j13 ;
   answer.at(1,58) = x20*j11+y20*j12+z20*j13 ; 

   answer.at(2,2)  = x1*j21+y1*j22+z1*j23 ;
   answer.at(2,5)  = x2*j21+y2*j22+z2*j23 ;
   answer.at(2,8)  = x3*j21+y3*j22+z3*j23 ;
   answer.at(2,11) = x4*j21+y4*j22+z4*j23 ;
   answer.at(2,14) = x5*j21+y5*j22+z5*j23 ;
   answer.at(2,17) = x6*j21+y6*j22+z6*j23 ;
   answer.at(2,20) = x7*j21+y7*j22+z7*j23 ;
   answer.at(2,23) = x8*j21+y8*j22+z8*j23 ;
   answer.at(2,26) = x9*j21+y9*j22+z9*j23 ;
   answer.at(2,29) = x10*j21+y10*j22+z10*j23 ;
   answer.at(2,32) = x11*j21+y11*j22+z11*j23 ;
   answer.at(2,35) = x12*j21+y12*j22+z12*j23 ;
   answer.at(2,38) = x13*j21+y13*j22+z13*j23 ; 
   answer.at(2,41) = x14*j21+y14*j22+z14*j23 ;
   answer.at(2,44) = x15*j21+y15*j22+z15*j23 ;
   answer.at(2,47) = x16*j21+y16*j22+z16*j23 ; 
   answer.at(2,50) = x17*j21+y17*j22+z17*j23 ;
   answer.at(2,53) = x18*j21+y18*j22+z18*j23 ;
   answer.at(2,56) = x19*j21+y19*j22+z19*j23 ;
   answer.at(2,59) = x20*j21+y20*j22+z20*j23 ; 

   answer.at(3,3)  = x1*j31+y1*j32+z1*j33 ;
   answer.at(3,6)  = x2*j31+y2*j32+z2*j33 ;
   answer.at(3,9)  = x3*j31+y3*j32+z3*j33 ;
   answer.at(3,12) = x4*j31+y4*j32+z4*j33 ;
   answer.at(3,15) = x5*j31+y5*j32+z5*j33 ;
   answer.at(3,18) = x6*j31+y6*j32+z6*j33 ;
   answer.at(3,21) = x7*j31+y7*j32+z7*j33 ;
   answer.at(3,24) = x8*j31+y8*j32+z8*j33 ;
   answer.at(3,27) = x9*j31+y9*j32+z9*j33 ;
   answer.at(3,30) = x10*j31+y10*j32+z10*j33 ;
   answer.at(3,33) = x11*j31+y11*j32+z11*j33 ;
   answer.at(3,36) = x12*j31+y12*j32+z12*j33 ;
   answer.at(3,39) = x13*j31+y13*j32+z13*j33 ; 
   answer.at(3,42) = x14*j31+y14*j32+z14*j33 ;
   answer.at(3,45) = x15*j31+y15*j32+z15*j33 ;
   answer.at(3,48) = x16*j31+y16*j32+z16*j33 ; 
   answer.at(3,51) = x17*j31+y17*j32+z17*j33 ;
   answer.at(3,54) = x18*j31+y18*j32+z18*j33 ;
   answer.at(3,57) = x19*j31+y19*j32+z19*j33 ;
   answer.at(3,60) = x20*j21+y20*j32+z20*j33 ; 


   answer.at(4,2)  = x1*j31+y1*j32+z1*j33 ;
   answer.at(4,5)  = x2*j31+y2*j32+z2*j33 ;
   answer.at(4,8)  = x3*j31+y3*j32+z3*j33 ;
   answer.at(4,11) = x4*j31+y4*j32+z4*j33 ;
   answer.at(4,14) = x5*j31+y5*j32+z5*j33 ;
   answer.at(4,17) = x6*j31+y6*j32+z6*j33 ;
   answer.at(4,20) = x7*j31+y7*j32+z7*j33 ;
   answer.at(4,23) = x8*j31+y8*j32+z8*j33 ;
   answer.at(4,26) = x9*j31+y9*j32+z9*j33 ;
   answer.at(4,29) = x10*j31+y10*j32+z10*j33 ;
   answer.at(4,32) = x11*j31+y11*j32+z11*j33 ;
   answer.at(4,35) = x12*j31+y12*j32+z12*j33 ;
   answer.at(4,38) = x13*j31+y13*j32+z13*j33 ; 
   answer.at(4,41) = x14*j31+y14*j32+z14*j33 ;
   answer.at(4,44) = x15*j31+y15*j32+z15*j33 ;
   answer.at(4,47) = x16*j31+y16*j32+z16*j33 ; 
   answer.at(4,50) = x17*j31+y17*j32+z17*j33 ;
   answer.at(4,53) = x18*j31+y18*j32+z18*j33 ;
   answer.at(4,56) = x19*j31+y19*j32+z19*j33 ;
   answer.at(4,59) = x20*j31+y20*j32+z20*j33 ; 
   answer.at(4,3)  = x1*j21+y1*j22+z1*j23 ;
   answer.at(4,6)  = x2*j21+y2*j22+z2*j23 ;
   answer.at(4,9)  = x3*j21+y3*j22+z3*j23 ;
   answer.at(4,12) = x4*j21+y4*j22+z4*j23 ;
   answer.at(4,15) = x5*j21+y5*j22+z5*j23 ;
   answer.at(4,18) = x6*j21+y6*j22+z6*j23 ;
   answer.at(4,21) = x7*j21+y7*j22+z7*j23 ;
   answer.at(4,24) = x8*j21+y8*j22+z8*j23 ;
   answer.at(4,27) = x9*j21+y9*j22+z9*j23 ;
   answer.at(4,30) = x10*j21+y10*j22+z10*j23 ;
   answer.at(4,33) = x11*j21+y11*j22+z11*j23 ;
   answer.at(4,36) = x12*j21+y12*j22+z12*j23 ;
   answer.at(4,39) = x13*j21+y13*j22+z13*j23 ; 
   answer.at(4,42) = x14*j21+y14*j22+z14*j23 ;
   answer.at(4,45) = x15*j21+y15*j22+z15*j23 ;
   answer.at(4,48) = x16*j21+y16*j22+z16*j23 ; 
   answer.at(4,51) = x17*j21+y17*j22+z17*j23 ;
   answer.at(4,54) = x18*j21+y18*j22+z18*j23 ;
   answer.at(4,57) = x19*j21+y19*j22+z19*j23 ;
   answer.at(4,60) = x20*j21+y20*j22+z20*j23 ; 

   answer.at(5,1)  = x1*j31+y1*j32+z1*j33 ;
   answer.at(5,4)  = x2*j31+y2*j32+z2*j33 ;
   answer.at(5,7)  = x3*j31+y3*j32+z3*j33 ;
   answer.at(5,10) = x4*j31+y4*j32+z4*j33 ;
   answer.at(5,13) = x5*j31+y5*j32+z5*j33 ;
   answer.at(5,16) = x6*j31+y6*j32+z6*j33 ;
   answer.at(5,19) = x7*j31+y7*j32+z7*j33 ;
   answer.at(5,22) = x8*j31+y8*j32+z8*j33 ;
   answer.at(5,25) = x9*j31+y9*j32+z9*j33 ;
   answer.at(5,28) = x10*j31+y10*j32+z10*j33 ;
   answer.at(5,31) = x11*j31+y11*j32+z11*j33 ;
   answer.at(5,34) = x12*j31+y12*j32+z12*j33 ;
   answer.at(5,37) = x13*j31+y13*j32+z13*j33 ; 
   answer.at(5,40) = x14*j31+y14*j32+z14*j33 ;
   answer.at(5,43) = x15*j31+y15*j32+z15*j33 ;
   answer.at(5,46) = x16*j31+y16*j32+z16*j33 ; 
   answer.at(5,49) = x17*j31+y17*j32+z17*j33 ;
   answer.at(5,52) = x18*j31+y18*j32+z18*j33 ;
   answer.at(5,55) = x19*j31+y19*j32+z19*j33 ;
   answer.at(5,58) = x20*j31+y20*j32+z20*j33 ; 
   answer.at(5,3)  = x1*j11+y1*j12+z1*j13 ;
   answer.at(5,6)  = x2*j11+y2*j12+z2*j13 ;
   answer.at(5,9)  = x3*j11+y3*j12+z3*j13 ;
   answer.at(5,12) = x4*j11+y4*j12+z4*j13 ;
   answer.at(5,15) = x5*j11+y5*j12+z5*j13 ;
   answer.at(5,18) = x6*j11+y6*j12+z6*j13 ;
   answer.at(5,21) = x7*j11+y7*j12+z7*j13 ;
   answer.at(5,24) = x8*j11+y8*j12+z8*j13 ;
   answer.at(5,27) = x9*j11+y9*j12+z9*j13 ;
   answer.at(5,30) = x10*j11+y10*j12+z10*j13 ;
   answer.at(5,33) = x11*j11+y11*j12+z11*j13 ;
   answer.at(5,36) = x12*j11+y12*j12+z12*j13 ;
   answer.at(5,39) = x13*j11+y13*j12+z13*j13 ; 
   answer.at(5,42) = x14*j11+y14*j12+z14*j13 ;
   answer.at(5,45) = x15*j11+y15*j12+z15*j13 ;
   answer.at(5,48) = x16*j11+y16*j12+z16*j13 ; 
   answer.at(5,51) = x17*j11+y17*j12+z17*j13 ;
   answer.at(5,54) = x18*j11+y18*j12+z18*j13 ;
   answer.at(5,57) = x19*j11+y19*j12+z19*j13 ;
   answer.at(5,60) = x20*j11+y20*j12+z20*j13 ; 


   answer.at(6,1)  = x1*j21+y1*j22+z1*j23 ;
   answer.at(6,4)  = x2*j21+y2*j22+z2*j23 ;
   answer.at(6,7)  = x3*j21+y3*j22+z3*j23 ;
   answer.at(6,10) = x4*j21+y4*j22+z4*j23 ;
   answer.at(6,13) = x5*j21+y5*j22+z5*j23 ;
   answer.at(6,16) = x6*j21+y6*j22+z6*j23 ;
   answer.at(6,19) = x7*j21+y7*j22+z7*j23 ;
   answer.at(6,22) = x8*j21+y8*j22+z8*j23 ;
   answer.at(6,25) = x9*j21+y9*j22+z9*j23 ;
   answer.at(6,28) = x10*j21+y10*j22+z10*j23 ;
   answer.at(6,31) = x11*j21+y11*j22+z11*j23 ;
   answer.at(6,34) = x12*j21+y12*j22+z12*j23 ;
   answer.at(6,37) = x13*j21+y13*j22+z13*j23 ; 
   answer.at(6,40) = x14*j21+y14*j22+z14*j23 ;
   answer.at(6,43) = x15*j21+y15*j22+z15*j23 ;
   answer.at(6,46) = x16*j21+y16*j22+z16*j23 ; 
   answer.at(6,49) = x17*j21+y17*j22+z17*j23 ;
   answer.at(6,52) = x18*j21+y18*j22+z18*j23 ;
   answer.at(6,55) = x19*j21+y19*j22+z19*j23 ;
   answer.at(6,58) = x20*j21+y20*j22+z20*j23 ; 
   answer.at(6,2)  = x1*j11+y1*j12+z1*j13 ;
   answer.at(6,5)  = x2*j11+y2*j12+z2*j13 ;
   answer.at(6,8)  = x3*j11+y3*j12+z3*j13 ;
   answer.at(6,11) = x4*j11+y4*j12+z4*j13 ;
   answer.at(6,14) = x5*j11+y5*j12+z5*j13 ;
   answer.at(6,17) = x6*j11+y6*j12+z6*j13 ;
   answer.at(6,20) = x7*j11+y7*j12+z7*j13 ;
   answer.at(6,23) = x8*j11+y8*j12+z8*j13 ;
   answer.at(6,26) = x9*j11+y9*j12+z9*j13 ;
   answer.at(6,29) = x10*j11+y10*j12+z10*j13 ;
   answer.at(6,32) = x11*j11+y11*j12+z11*j13 ;
   answer.at(6,35) = x12*j11+y12*j12+z12*j13 ;
   answer.at(6,38) = x13*j11+y13*j12+z13*j13 ; 
   answer.at(6,41) = x14*j11+y14*j12+z14*j13 ;
   answer.at(6,44) = x15*j11+y15*j12+z15*j13 ;
   answer.at(6,47) = x16*j11+y16*j12+z16*j13 ; 
   answer.at(6,50) = x17*j11+y17*j12+z17*j13 ;
   answer.at(6,53) = x18*j11+y18*j12+z18*j13 ;
   answer.at(6,56) = x19*j11+y19*j12+z19*j13 ;
   answer.at(6,59) = x20*j11+y20*j12+z20*j13 ; 

   return  ;
}

   
void
QSpace :: computeNmatrixAt (GaussPoint* aGaussPoint, FloatMatrix& answer) 
   // Returns the displacement interpolation matrix {N} of the receiver, eva-
   // luated at aGaussPoint.
{
   double    x,y,z,
             n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
             n16,n17,n18,n19,n20;
   // FloatMatrix* answer ;

   x = aGaussPoint -> giveCoordinate(1) ;
   y = aGaussPoint -> giveCoordinate(2) ;
   z = aGaussPoint -> giveCoordinate(3) ;

   n1  = 0.125*(1.+x)*(1.+y)*(1.+z)*(x+y+z-2.);
   n2  = 0.25*(1.-x*x)*(1.+y)*(1.+z);
   n3  = 0.125*(1.-x)*(1.+y)*(1.+z)*(-x+y+z-2.);
   n4  = 0.25*(1.-x)*(1.-y*y)*(1.+z);
   n5  = 0.125*(1.-x)*(1.-y)*(1.+z)*(-x-y+z-2.);
   n6  = 0.25*(1.-x*x)*(1.-y)*(1.+z);
   n7  = 0.125*(1.+x)*(1.-y)*(1.+z)*(x-y-z-2.);
   n8  = 0.25*(1.+x)*(1.-y*y)*(1.+z);
   n9  = 0.25*(1.+x)*(1.+y)*(1.-z*z);
   n10 = 0.25*(1.-x)*(1.+y)*(1.-z*z);
   n11 = 0.25*(1.-x)*(1.-y)*(1.-z*z);
   n12 = 0.25*(1.+x)*(1.-y)*(1.-z*z);
   n13 = 0.125*(1.+x)*(1.+y)*(1.-z)*(x+y-z-2.);
   n14 = 0.25*(1.-x*x)*(1.+y)*(1.-z);
   n15 = 0.125*(1.-x)*(1.+y)*(1.-z)*(-x+y-z-2.);
   n16 = 0.25*(1.-x)*(1.-y*y)*(1.-z);
   n17 = 0.125*(1.-x)*(1.-y)*(1.-z)*(-x-y-z-2.);
   n18 = 0.25*(1.-x*x)*(1.-y)*(1.-z);
   n19 = 0.125*(1.+x)*(1.-y)*(1.-z)*(x-y-z-2.);
   n20 = 0.25*(1.+x)*(1.-y*y)*(1.-z);

   //answer = new FloatMatrix(3,60) ;
  answer.resize (3,60);
  answer.zero();

   answer.at(1,1)  = n1  ;
   answer.at(1,4)  = n2  ;
   answer.at(1,7)  = n3  ; 
   answer.at(1,10) = n4  ;
   answer.at(1,13) = n5  ;
   answer.at(1,16) = n6  ;
   answer.at(1,19) = n7  ;
   answer.at(1,22) = n8  ;
   answer.at(1,25) = n9  ;
   answer.at(1,28) = n10 ;
   answer.at(1,31) = n11 ; 
   answer.at(1,34) = n12 ;
   answer.at(1,37) = n13 ;
   answer.at(1,40) = n14 ;
   answer.at(1,43) = n15 ;
   answer.at(1,46) = n16 ;
   answer.at(1,49) = n17 ;
   answer.at(1,52) = n18 ;
   answer.at(1,55) = n19 ; 
   answer.at(1,58) = n20 ;

   answer.at(2,2)  = n1 ;
   answer.at(2,5)  = n2 ;
   answer.at(2,8)  = n3 ;
   answer.at(2,11) = n4 ;
   answer.at(2,14) = n5 ;
   answer.at(2,17) = n6 ;
   answer.at(2,20) = n7 ;
   answer.at(2,23) = n8 ;
   answer.at(2,26) = n9 ;
   answer.at(2,29) = n10 ;
   answer.at(2,32) = n11 ;
   answer.at(2,35) = n12 ;
   answer.at(2,38) = n13 ;
   answer.at(2,41) = n14 ;
   answer.at(2,44) = n15 ;
   answer.at(2,47) = n16 ;
   answer.at(2,50) = n17 ;
   answer.at(2,53) = n18 ;
   answer.at(2,56) = n19 ;
   answer.at(2,59) = n20 ;

   answer.at(3,3)  = n1 ;
   answer.at(3,6)  = n2 ;
   answer.at(3,9)  = n3 ;
   answer.at(3,12) = n4 ;
   answer.at(3,15) = n5 ;
   answer.at(3,18) = n6 ;
   answer.at(3,21) = n7 ;
   answer.at(3,24) = n8 ;
   answer.at(3,27) = n9 ;
   answer.at(3,30) = n10 ;
   answer.at(3,33) = n11 ;
   answer.at(3,36) = n12 ;
   answer.at(3,39) = n13 ;
   answer.at(3,42) = n14 ;
   answer.at(3,45) = n15 ;
   answer.at(3,48) = n16 ;
   answer.at(3,51) = n17 ;
   answer.at(3,54) = n18 ;
   answer.at(3,57) = n19 ;
   answer.at(3,60) = n20 ;

   return  ;
}

void  QSpace :: computeGaussPoints ()
   // Sets up the array containing the four Gauss points of the receiver.
{
   numberOfIntegrationRules = 1 ;
  integrationRulesArray = new IntegrationRule*;
  integrationRulesArray[0] = new GaussIntegrationRule (1,domain,1, 6);
  integrationRulesArray[0]->setUpIntegrationPoints (_Cube, numberOfGaussPoints, this,  _3dMat);

}

double  
QSpace :: computeVolumeAround (GaussPoint* aGaussPoint)
   // Returns the portion of the receiver which is attached to aGaussPoint.
{
   FloatMatrix  jacMtrx ;
   FloatArray*  coord ;
   double       determinant,weight,volume ;

   coord       = aGaussPoint -> giveCoordinates() ;
   this -> computeJacobianMatrixAt(jacMtrx, coord) ;
   determinant = fabs (jacMtrx.giveDeterminant()) ;
   weight      = aGaussPoint -> giveWeight() ;

   volume      = determinant * weight ;


   return volume ;
}

void
QSpace :: computeJacobianMatrixAt (FloatMatrix& answer, FloatArray* coord) 
   // Returns the jacobian matrix  J (x,y,z)/(ksi,eta,dzeta)  of the receiver.
   // Computes it if it does not exist yet.
{
   Node    *n1,*n2,*n3,*n4,*n5,*n6,*n7,*n8,*n9,*n10,*n11,*n12,*n13,*n14,*n15,
           *n16,*n17,*n18,*n19,*n20 ;

   double u,v,w ;

  answer.resize(3,3); answer.zero();

  // nodes
   n1  = this -> giveNode(1)  ;
   n2  = this -> giveNode(2)  ;
   n3  = this -> giveNode(3)  ;
   n4  = this -> giveNode(4)  ;
   n5  = this -> giveNode(5)  ;
   n6  = this -> giveNode(6)  ;
   n7  = this -> giveNode(7)  ;
   n8  = this -> giveNode(8)  ;
   n9  = this -> giveNode(9)  ;
   n10 = this -> giveNode(10) ;
   n11 = this -> giveNode(11) ;
   n12 = this -> giveNode(12) ;
   n13 = this -> giveNode(13) ;
   n14 = this -> giveNode(14) ;
   n15 = this -> giveNode(15) ;
   n16 = this -> giveNode(16) ;
   n17 = this -> giveNode(17) ;
   n18 = this -> giveNode(18) ;
   n19 = this -> giveNode(19) ;
   n20 = this -> giveNode(20) ;

// coordinates
 FloatArray x(numberOfDofMans) , y(numberOfDofMans) , z(numberOfDofMans) ;

   x.at(1)  = n1  -> giveCoordinate(1) ;
   x.at(2)  = n2  -> giveCoordinate(1) ;
   x.at(3)  = n3  -> giveCoordinate(1) ;
   x.at(4)  = n4  -> giveCoordinate(1) ;
   x.at(5)  = n5  -> giveCoordinate(1) ;
   x.at(6)  = n6  -> giveCoordinate(1) ;
   x.at(7)  = n7  -> giveCoordinate(1) ;
   x.at(8)  = n8  -> giveCoordinate(1) ;
   x.at(9)  = n9  -> giveCoordinate(1) ;
   x.at(10) = n10 -> giveCoordinate(1) ;
   x.at(11) = n11 -> giveCoordinate(1) ;
   x.at(12) = n12 -> giveCoordinate(1) ;
   x.at(13) = n13 -> giveCoordinate(1) ;
   x.at(14) = n14 -> giveCoordinate(1) ;
   x.at(15) = n15 -> giveCoordinate(1) ;
   x.at(16) = n16 -> giveCoordinate(1) ;
   x.at(17) = n17 -> giveCoordinate(1) ;
   x.at(18) = n18 -> giveCoordinate(1) ;
   x.at(19) = n19 -> giveCoordinate(1) ;
   x.at(20) = n20 -> giveCoordinate(1) ;

   y.at(1)  = n1  -> giveCoordinate(2) ;
   y.at(2)  = n2  -> giveCoordinate(2) ;
   y.at(3)  = n3  -> giveCoordinate(2) ;
   y.at(4)  = n4  -> giveCoordinate(2) ;
   y.at(5)  = n5  -> giveCoordinate(2) ;
   y.at(6)  = n6  -> giveCoordinate(2) ;
   y.at(7)  = n7  -> giveCoordinate(2) ;
   y.at(8)  = n8  -> giveCoordinate(2) ;
   y.at(9)  = n9  -> giveCoordinate(2) ;
   y.at(10) = n10 -> giveCoordinate(2) ;
   y.at(11) = n11 -> giveCoordinate(2) ;
   y.at(12) = n12 -> giveCoordinate(2) ;
   y.at(13) = n13 -> giveCoordinate(2) ;
   y.at(14) = n14 -> giveCoordinate(2) ;
   y.at(15) = n15 -> giveCoordinate(2) ;
   y.at(16) = n16 -> giveCoordinate(2) ;
   y.at(17) = n17 -> giveCoordinate(2) ;
   y.at(18) = n18 -> giveCoordinate(2) ;
   y.at(19) = n19 -> giveCoordinate(2) ;
   y.at(20) = n20 -> giveCoordinate(2) ;

   z.at(1)  = n1  -> giveCoordinate(3) ;
   z.at(2)  = n2  -> giveCoordinate(3) ;
   z.at(3)  = n3  -> giveCoordinate(3) ;
   z.at(4)  = n4  -> giveCoordinate(3) ;
   z.at(5)  = n5  -> giveCoordinate(3) ;
   z.at(6)  = n6  -> giveCoordinate(3) ;
   z.at(7)  = n7  -> giveCoordinate(3) ;
   z.at(8)  = n8  -> giveCoordinate(3) ;
   z.at(9)  = n9  -> giveCoordinate(3) ;
   z.at(10) = n10 -> giveCoordinate(3) ;
   z.at(11) = n11 -> giveCoordinate(3) ;
   z.at(12) = n12 -> giveCoordinate(3) ;
   z.at(13) = n13 -> giveCoordinate(3) ;
   z.at(14) = n14 -> giveCoordinate(3) ;
   z.at(15) = n15 -> giveCoordinate(3) ;
   z.at(16) = n16 -> giveCoordinate(3) ;
   z.at(17) = n17 -> giveCoordinate(3) ;
   z.at(18) = n18 -> giveCoordinate(3) ;
   z.at(19) = n19 -> giveCoordinate(3) ;
   z.at(20) = n20 -> giveCoordinate(3) ;

   u = coord->at(1);
   v = coord->at(2);
   w = coord->at(3);

   FloatArray dx(numberOfDofMans) , dy(numberOfDofMans) , dz(numberOfDofMans) ;

   dx.at(1)  =  0.125*(1.+v)*(1.+w)*(2.*u+v+w-1.);
   dx.at(2)  = -0.5*u*(1.+v)*(1.+w);
   dx.at(3)  = -0.125*(1.+v)*(1.+w)*(-2.*u+v+w-1.);
   dx.at(4)  = -0.25*(1.-v*v)*(1.+w);
   dx.at(5)  = -0.125*(1.-v)*(1.+w)*(-2.*u-v+w-1.);
   dx.at(6)  = -0.5*u*(1.-v)*(1.+w);
   dx.at(7)  =  0.125*(1.-v)*(1.+w)*(2.*u-v+w-1.);
   dx.at(8)  =  0.25*(1.-v*v)*(1.+w);
   dx.at(9)  =  0.25*(1.+v)*(1.-w*w);
   dx.at(10) = -0.25*(1.+v)*(1.-w*w);
   dx.at(11) = -0.25*(1.-v)*(1.-w*w);
   dx.at(12) =  0.25*(1.-v)*(1.-w*w);
   dx.at(13) =  0.125*(1.+v)*(1.-w)*(2.*u+v-w-1.);
   dx.at(14) = -0.5*u*(1.+v)*(1.-w);
   dx.at(15) = -0.125*(1.+v)*(1.-w)*(-2.*u+v-w-1.);
   dx.at(16) = -0.25*(1.-v*v)*(1.-w);
   dx.at(17) = -0.125*(1.-v)*(1.-w)*(-2.*u-v-w-1.);
   dx.at(18) = -0.5*u*(1.-v)*(1.-w);
   dx.at(19) =  0.125*(1.-v)*(1.-w)*(2.*u-v-w-1.);
   dx.at(20) =  0.25*(1.-v*v)*(1.-w);

   dy.at(1)  =  0.125*(1.+u)*(1.+w)*(u+2.*v+w-1.);
   dy.at(2)  =  0.25*(1.-u*u)*(1.+w);
   dy.at(3)  =  0.125*(1.-u)*(1.+w)*(-u+2.*v+w-1.);
   dy.at(4)  = -0.5*v*(1.-u)*(1.+w);
   dy.at(5)  = -0.125*(1.-u)*(1.+w)*(-u-2.*v+w-1.);
   dy.at(6)  = -0.25*(1.-u*u)*(1.+w);
   dy.at(7)  = -0.125*(1.+u)*(1.+w)*(u-2.*v+w-1.);
   dy.at(8)  = -0.5*v*(1.+u)*(1.+w);
   dy.at(9)  =  0.25*(1.+u)*(1.-w*w);
   dy.at(10) =  0.25*(1.-u)*(1.-w*w);
   dy.at(11) = -0.25*(1.-u)*(1.-w*w);
   dy.at(12) = -0.25*(1.+u)*(1.-w*w);
   dy.at(13) =  0.125*(1.+u)*(1.-w)*(u+2.*v-w-1.);
   dy.at(14) =  0.25*(1.-u*u)*(1.-w);
   dy.at(15) =  0.125*(1.-u)*(1.-w)*(-u+2.*v-w-1.);
   dy.at(16) = -0.5*v*(1.-u)*(1.-w);
   dy.at(17) = -0.125*(1.-u)*(1.-w)*(-u-2.*v-w-1.);
   dy.at(18) = -0.25*(1.-u*u)*(1.-w);
   dy.at(19) = -0.125*(1.+u)*(1.-w)*(u-2.*v-w-1.);
   dy.at(20) = -0.5*v*(1.+u)*(1.-w);

   dz.at(1)  =  0.125*(1.+u)*(1.+v)*(u+v+2.*w-1.);
   dz.at(2)  =  0.25*(1.-u*u)*(1.+v);
   dz.at(3)  =  0.125*(1.-u)*(1.+v)*(-u+v+2.*w-1.);
   dz.at(4)  =  0.25*(1.-u)*(1.-v*v);
   dz.at(5)  =  0.125*(1.-u)*(1.-v)*(-u-v+2.*w-1.);
   dz.at(6)  =  0.25*(1.-u*u)*(1.-v);
   dz.at(7)  =  0.125*(1.+u)*(1.-v)*(u-v+2.*w-1.);
   dz.at(8)  =  0.25*(1.+u)*(1.-v*v);
   dz.at(9)  = -0.5*w*(1.+u)*(1.+v);
   dz.at(10) = -0.5*w*(1.-u)*(1.+v);
   dz.at(11) = -0.5*w*(1.-u)*(1.-v);
   dz.at(12) = -0.5*w*(1.+u)*(1.-v);
   dz.at(13) = -0.125*(1.+u)*(1.+v)*(u+v-2.*w-1.);
   dz.at(14) = -0.25*(1.-u*u)*(1.+v);
   dz.at(15) = -0.125*(1.-u)*(1.+v)*(-u+v-2.*w-1.);
   dz.at(16) = -0.25*(1.-u)*(1.-v*v);
   dz.at(17) = -0.125*(1.-u)*(1.-v)*(-u-v-2.*w-1.);
   dz.at(18) = -0.25*(1.-u*u)*(1.-v);
   dz.at(19) = -0.125*(1.+u)*(1.-v)*(u-v-2.*w-1.);
   dz.at(20) = -0.25*(1.+u)*(1.-v*v);


   for (int i = 1; i <= numberOfDofMans ; i++)
   {
     answer.at(1,1) += dx.at(i)*x.at(i);
     answer.at(1,2) += dx.at(i)*y.at(i);
     answer.at(1,3) += dx.at(i)*z.at(i);
     answer.at(2,1) += dy.at(i)*x.at(i);
     answer.at(2,2) += dy.at(i)*y.at(i);
     answer.at(2,3) += dy.at(i)*z.at(i);
     answer.at(3,1) += dz.at(i)*x.at(i);
     answer.at(3,2) += dz.at(i)*y.at(i);
     answer.at(3,3) += dz.at(i)*z.at(i);
   }
}


IRResultType
QSpace :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

  this->Element :: initializeFrom (ir);
 numberOfGaussPoints = 8;
 IR_GIVE_OPTIONAL_FIELD (ir, numberOfGaussPoints, IFT_QSpace_nip, "nip"); // Macro

  if (!((numberOfGaussPoints == 8) || (numberOfGaussPoints == 27))) numberOfGaussPoints = 8; 
  // set - up Gaussian integration points
  this -> computeGaussPoints();

  return IRRT_OK;
} 


void
QSpace ::   giveDofManDofIDMask  (int inode, EquationID, IntArray& answer) const {
// returns DofId mask array for inode element node.
// DofId mask array determines the dof ordering requsted from node.
// DofId mask array contains the DofID constants (defined in cltypes.h)
// describing physical meaning of particular DOFs.
 //IntArray* answer = new IntArray (3);
 answer.resize (3);

 answer.at(1) = D_u;
 answer.at(2) = D_v;
 answer.at(3) = D_w;

 return ;
}

double
QSpace::giveCharacteristicLenght(GaussPoint* gp, const FloatArray &normalToCrackPlane) 
{
 double factor = __OOFEM_POW((double)this->numberOfGaussPoints, 1./3.);
 return this -> giveLenghtInDir (normalToCrackPlane)/factor;
}

int
QSpace :: computeGlobalCoordinates (FloatArray& answer, const FloatArray& lcoords) 
{
 int i;
 double    x,y,z;
 FloatArray n(20);

 x = lcoords.at(1);
 y = lcoords.at(2);
 z = lcoords.at(3);
 
 n.at(1)  = 0.125*(1.+x)*(1.+y)*(1.+z)*(x+y+z-2.);
 n.at(2)  = 0.25*(1.-x*x)*(1.+y)*(1.+z);
 n.at(3)  = 0.125*(1.-x)*(1.+y)*(1.+z)*(-x+y+z-2.);
 n.at(4)  = 0.25*(1.-x)*(1.-y*y)*(1.+z);
 n.at(5)  = 0.125*(1.-x)*(1.-y)*(1.+z)*(-x-y+z-2.);
 n.at(6)  = 0.25*(1.-x*x)*(1.-y)*(1.+z);
 n.at(7)  = 0.125*(1.+x)*(1.-y)*(1.+z)*(x-y-z-2.);
 n.at(8)  = 0.25*(1.+x)*(1.-y*y)*(1.+z);
 n.at(9)  = 0.25*(1.+x)*(1.+y)*(1.-z*z);
 n.at(10) = 0.25*(1.-x)*(1.+y)*(1.-z*z);
 n.at(11) = 0.25*(1.-x)*(1.-y)*(1.-z*z);
 n.at(12) = 0.25*(1.+x)*(1.-y)*(1.-z*z);
 n.at(13) = 0.125*(1.+x)*(1.+y)*(1.-z)*(x+y-z-2.);
 n.at(14) = 0.25*(1.-x*x)*(1.+y)*(1.-z);
 n.at(15) = 0.125*(1.-x)*(1.+y)*(1.-z)*(-x+y-z-2.);
 n.at(16) = 0.25*(1.-x)*(1.-y*y)*(1.-z);
 n.at(17) = 0.125*(1.-x)*(1.-y)*(1.-z)*(-x-y-z-2.);
 n.at(18) = 0.25*(1.-x*x)*(1.-y)*(1.-z);
 n.at(19) = 0.125*(1.+x)*(1.-y)*(1.-z)*(x-y-z-2.);
 n.at(20) = 0.25*(1.+x)*(1.-y*y)*(1.-z);

 answer.resize (3); answer.zero();
 for (i=1; i<=20; i++) {
  answer.at(1) += n.at(i)*this->giveNode(i)->giveCoordinate(1);
  answer.at(2) += n.at(i)*this->giveNode(i)->giveCoordinate(2);
  answer.at(3) += n.at(i)*this->giveNode(i)->giveCoordinate(3);
 }

 return 1;
}
