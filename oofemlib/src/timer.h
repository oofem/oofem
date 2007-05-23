/* $Header: /home/cvs/bp/oofem/oofemlib/src/combuff.h,v 1.5 2003/04/06 14:08:23 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2006   Borek Patzak                                       



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
#ifndef timer_h
#include "clock.h"

/**
   Class implementing single timer, providing wall clock and user time capabilities.
*/
class Timer {
  // wall clock time structures
  time_t start_wtime, end_wtime;
  // user time struct
  oofem_timeval start_utime, end_utime; 

 public:
  void startTimer () {start_wtime = ::getTime(); ::getUtime(start_utime);}
  void stopTimer  () {end_wtime = ::getTime(); ::getUtime(end_utime);}

  /// Returns total user time elapsed in seconds
  long int getUtime  () const {return (long int) (end_utime.tv_sec-start_utime.tv_sec);}
  /// Returns total elapsed wall clock time in seconds
  long int getWtime ()  const {return end_wtime-start_wtime;}

  // converts total seconds into hours, mins, and seconds
  void convert2HMS (int &nhrs, int &nmin, int &nsec, long int tsec) const {::convertTS2HMS (nhrs, nmin, nsec, tsec);}
  
  /// short prints receiver state into a string
  void sprintYourselfShort (char* buff) const {sprintf (buff, "ut: %lds, wt: %lds", getUtime(), getWtime());}
};
  
  

/**
   Timer class, assumed to be an attribute of engng model, serving stop-watch facility for engng model.
   It can handle several timers independently, each corresponding to different solution stage, etc.
   Each timer is capable to track elapsed wall clock time as well as user time.
 */
class EngngModelTimer
{
 public:
  enum EngngModelTimerType {EMTT_AnalysisTimer, EMTT_SolutionStepTimer, EMTT_LoadBallancingTimer, EMTT_DataTransferTimer, EMTT_LastTimer};
 protected:

  /// Array of Timer classes.
  Timer timers[EMTT_LastTimer];

 public:

  EngngModelTimer () {}
  virtual ~EngngModelTimer () {}
  
  

  /**@name Profiling routines */
  //@{
  void startTimer (EngngModelTimerType t) {timers[t].startTimer();}
  void stopTimer  (EngngModelTimerType t) {timers[t].stopTimer();}
  //@}
  
  /** Reporting routines */
  //@{
  /// Returns total user time elapsed 
  long int getUtime  (EngngModelTimerType t) const {return timers[t].getUtime();}
  /// Returns elapsed wall clock time
  long int getWtime (EngngModelTimerType t) const {return timers[t].getWtime();}
  /// Returns pointer to timer determined by EngngModelTimerType
  const Timer* getTimer (EngngModelTimerType t) const {return timers+t;}
  /// converts total seconds into hours, mins, and seconds
  void convert2HMS (int &nhrs, int &nmin, int &nsec, long int tsec) const {::convertTS2HMS (nhrs, nmin, nsec, tsec);}
  //@}
};

#define timer_h
#endif
