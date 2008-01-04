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
   Author: Richard Vondracek, <richard.vondracek@seznam.cz>
*/


// DSSAfx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
//
#if !defined(_DSSAFX_H__)
#define _DSSAFX_H__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifdef __WATCOMC__
 #ifndef BOOL_DEF
  enum bool {false = 0, true = 1};
  #define BOOL_DEF
 #endif 
#endif

// This is just for compilation with MFC
#ifdef BUILDING_CCFEModel
   #include <afx.h>
   #ifdef _DEBUG
      #define new DEBUG_NEW
      #undef THIS_FILE
      static char THIS_FILE[] = __FILE__;
   #endif
#endif


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>


#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif

typedef int                 BOOL;
typedef unsigned char       BYTE;
typedef unsigned long ULONG;

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef PURE 
#define PURE =0
#endif

#define DSS_NAMESPASE_BEGIN    
#define DSS_NAMESPASE_END 
//#define DSS_NAMESPASE_BEGIN namespace DSS {
//#define DSS_NAMESPASE_END }

#endif // !defined(_DSSAFX_H__)

