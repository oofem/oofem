/* $Header: /home/cvs/bp/oofem/oofemlib/src/clock.C,v 1.5 2003/05/19 13:03:57 bp Exp $ */
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

/*  file CLOCK.C  */

/*
  This file customizes function 'timeNow' to your platform.
  Function 'timeNow' is used for measuring elapsed or CPU time.
*/

#include "clock.h"
#include "compiler.h"
#ifndef _MSC_VER

void getUtime (oofem_timeval& answer)
{
 struct rusage rsg;
 getrusage (RUSAGE_SELF, &rsg);
 answer.tv_sec  = rsg.ru_utime.tv_sec;
 answer.tv_usec = rsg.ru_utime.tv_usec;
}

void getRelativeUtime (oofem_timeval& answer, oofem_timeval& from)
{
 struct rusage rsg;
 getrusage (RUSAGE_SELF, &rsg);

 if (rsg.ru_utime.tv_usec < from.tv_usec) {

  answer.tv_usec = (OOFEM_USEC_LIM-from.tv_usec) + rsg.ru_utime.tv_usec;
  answer.tv_sec  = rsg.ru_utime.tv_sec-from.tv_sec-1;

 } else {

  answer.tv_usec = rsg.ru_utime.tv_usec-from.tv_usec;
  answer.tv_sec  = rsg.ru_utime.tv_sec-from.tv_sec;

 }

}


time_t getTime ()
{
  time_t t;
  t = time(NULL);
  return t;
}
 
#else // #ifndef _MSC_VER

void getUtime (oofem_timeval& answer)
{
	clock_t	utime = clock();
	answer.tv_sec = utime/CLOCKS_PER_SEC;
	answer.tv_usec = 0;
}

void getRelativeUtime (oofem_timeval& answer, oofem_timeval& from)
{
	clock_t	utime = clock();

	answer.tv_sec = utime/CLOCKS_PER_SEC - from.tv_sec;
	answer.tv_usec = 0;
}


time_t getTime ()
{
  time_t t;
  t = time(NULL);
  return t;
}

#endif 
