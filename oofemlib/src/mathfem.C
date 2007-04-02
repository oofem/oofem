/* $Header: /home/cvs/bp/oofem/oofemlib/src/mathfem.C,v 1.9 2003/04/06 14:08:25 bp Exp $ */
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


#ifndef __MAKEDEPEND
#include <math.h>
#endif
#include "mathfem.h"

// measure dependent constant
#define CUBIC_ZERO 1.0e-100


void cubic (double a, double b, double c, double d, double*r1, double *r2, double *r3,int *num)
//
//   solves cubic equation for real roots
//
//   input:
//     a,b,c,d - coefficients of equation in form:
//     ax^3 + bx^2 + cx + d = 0
//
//   output:
//     r1,r2,r3 - roots (only first num roots is valid)
//     num      - number of roots resolved
//
{
 double  aa, p, q, D, pq, u, v, phi;
 double help;
 
 if(fabs(a) < CUBIC_ZERO){
  if(fabs(b) < CUBIC_ZERO){
   *r1 = -d / c;
   *num = 1;
  }
  else{
   if((D = c * c - 4.0 * b * d) < 0.0) {
    *num = 0;
    return;
   }else{
     //*r1 = (-c + sqrt(D)) / 2.0 / b;
     //*r2 = (-c - sqrt(D)) / 2.0 / b;
     help = -(c+sgn(c)*sqrt(D))/2.0;
     *r1 = help/b;
     *r2 = d/help;
     *num = 2;
   }
  }
 }
 else {
  aa = a;
  a = b/aa;
  b = c/aa;
  c = d/aa;
  p = b - a * a / 3.0;
  q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;
  pq = p * q;
  D = q * q / 4.0 + p * p * p / 27.0;
  if(fabs(D) < CUBIC_ZERO){
   if(fabs(pq) < CUBIC_ZERO){
    *r1 = 0.0 - a / 3.0;
    *r2 = *r1;
    *r3 = *r1;
    *num = 3;
   }
   else {
    if(q < 0.0)
     *r2 = -exp(log(-q / 2.0) / 3.0);
    else
     *r2 = exp(log(q / 2.0) / 3.0);
    *r1 = -2.0 * *r2 - a / 3.0;
    *r2 -= a / 3.0;
    *num = 2;
   }
  }
  else{
   if(D > 0.0){
    u = -q / 2.0 + sqrt(D);
    v = -q / 2.0 - sqrt(D);
    if(u < 0.0)
     u = -exp(log(-u) / 3.0);
    else
     u = exp(log(u) / 3.0);
    if(v < 0.0)
     v = -exp(log(-v) / 3.0);
    else
     v = exp(log(v) / 3.0);
    *r1 = u + v - a / 3.0;
    *num = 1;
   }
   else {
    p = sqrt(fabs(p) / 3.0);
    help = (-q / (2.0 * p * p * p));
    if (fabs (help) > 1.0) help = sgn(help); // prevent rounding errors
    
    phi = acos(help) / 3.0;
    *r1 = 2.0 * p * cos(phi) - a / 3.0;
    *r2 = -2.0 * p * cos(phi - M_PI / 3.0) - a / 3.0;
    *r3 = -2.0 * p * cos(phi + M_PI / 3.0) - a / 3.0;
    *num = 3;
   }
  }
 }
 return ;
} 


void cubic3r (double a, double b, double c, double d, double*r1, double *r2, double *r3,int *num)
//
//   solves cubic equation for real roots
//
//   input:
//     a,b,c,d - coefficients of equation in form:
//     ax^3 + bx^2 + cx + d = 0
//
//   output:
//     r1,r2,r3 - roots (only first num roots is valid)
//     num      - number of roots resolved
//
{
 double  aa, r, q, p, D, phi;
 double help;
 
 if(fabs(a) < CUBIC_ZERO){
  if(fabs(b) < CUBIC_ZERO){
   *r1 = -d / c;
   *num = 1;
  }
  else{
   if((D = c * c - 4.0 * b * d) < 0.0) {
    *num = 0;
    return;
   }else{
     //*r1 = (-c + sqrt(D)) / 2.0 / b;
     //*r2 = (-c - sqrt(D)) / 2.0 / b;
     help=-(c+sgn(c)*sqrt(D))/2.0;
     *r1 = help/b;
     *r2 = d/help;
     *num = 2;
    return;
   }
  }
 } else {
  aa = a;
  a = b/aa;
  b = c/aa;
  c = d/aa;
  
  q = (a*a-3.0*b)/9.0;
  r = (2.0*a*a*a-9.0*a*b+27.0*c)/54.0;
  
  D = q*q*q-r*r;
//  if (D > 0.) {
   // three real roots
    help = r/sqrt(q*q*q);
    if (fabs (help) > 1.0) help = sgn(help); // prevent rounding errors
   phi = acos (help);
   p = sqrt(q);
   
   *r1 = -2.0*p*cos(phi/3.0)-a/3.0;
   *r2 = -2.0*p*cos((phi+2.0*M_PI)/3.0)-a/3.0;
   *r3 = -2.0*p*cos((phi-2.0*M_PI)/3.0)-a/3.0;
   *num = 3;
/*  } else {
   
   help = fabs(r) + sqrt(-D);
   A = -sgn(r)*__OOFEM_POW(D, 1./3.);
   if (fabs(A) > CUBIC_ZERO) 
    B = q/A;
   else
    B = 0.0;
   
   *r1 = (A+B) - a/3.0;
   *num = 1;
  }
*/
  return ;
 } 
}



int iperm (int val, int rank)
//
// returns iperm of val, in specific rank
//
// rank 3 : val is valid from 1,2,3 and returns for val 1,2,3 2,3,1
// rank 2 : -||-
{
  return ( val+1 > rank ?  1 : val+1);
}

/*
#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#endif

// matherr - error-handling function
// When an math error occurs, a pointer to the exception
// structure x will be passed to the user-supplied matherr function.  This
// structure, which is defined in the math.h header file.
int matherr(struct __exception *x )
{
 char *sterr;
 
 switch (x->type) {
  
 case DOMAIN:
 case OVERFLOW:
 case UNDERFLOW:
  
  sterr = strerror(EDOM);
  fprintf(stderr, "%s failed: %s\n", x->name, sterr);
  abort();

 default:
  sterr = strerror(EDOM);
  fprintf(stderr, "%s failed: %s\n", x->name, sterr);

  return (0);    // libm prints error message and sets errno 
 }
}
*/
