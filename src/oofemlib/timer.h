/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef timer_h
#define timer_h

#include "oofemcfg.h"

#include <chrono>

namespace oofem {
/**
 * Class implementing single timer, providing wall clock and user time capabilities.
 */
class OOFEM_EXPORT Timer
{
    /// Wall clock time markers.
    std :: chrono :: time_point< std :: chrono :: high_resolution_clock >start_wtime, end_wtime;
    /// User time.
    std :: chrono :: duration< double >start_utime, end_utime;
    /// Accumulated wtime and utime (in seconds) from start.
    std :: chrono :: duration< double >elapsedWTime, elapsedUTime;
    /// Flag indicating whether timer is running.
    bool running;

public:
    Timer();

    void startTimer();
    void stopTimer();
    void pauseTimer();
    void resumeTimer();
    void initTimer();
    bool isRunning() { return running; }

    /// Returns total user time elapsed in seconds
    double getUtime();
    /// Returns total elapsed wall clock time in seconds
    double getWtime();

    /**
     * Converts total seconds into hours, minutes and remaining seconds
     * @param[out] nhrs Number of hours.
     * @param[out] nmin Number of minutes.
     * @param[out] nsec Number of seconds.
     * @param[in] tsec Total time in seconds.
     */
    static void convert2HMS(int &nhrs, int &nmin, int &nsec, double tsec);

    /// Prints receiver state into a string.
    void toString(char *buff);

    void updateElapsedTime();

private:
    /// Platform independent wrapper for user time
    void getUtime(std :: chrono :: duration< double > &answer);
    /// Platform independent wrapper for wall time
    void getTime(std :: chrono :: time_point< std :: chrono :: high_resolution_clock > &answer);
};

/**
 * Timer class, assumed to be an attribute of engineering model, serving stop-watch facility for engineering model.
 * It can handle several timers independently, each corresponding to different solution stage, etc.
 * Each timer is capable to track elapsed wall clock time as well as user time.
 */
class OOFEM_EXPORT EngngModelTimer
{
public:
    /**
     * Enumeration to distinguish different type of timers.
     *
     * EMTT_NetComputationalStepTimer timer (and particularly its wall clock time) should measure only computation itself,
     * no communication, therefore it should be measure of workload (in terms of wall clock time) on particular processors.
     * It also typically not include time needed to solve the system of equations, since this has to be done in parallel,
     * so solution takes the same time on all processors and include unwanted synchronization.
     */
    enum EngngModelTimerType {
        EMTT_AnalysisTimer,
        EMTT_SolutionStepTimer,
        EMTT_NetComputationalStepTimer,
        EMTT_LoadBalancingTimer,
        EMTT_DataTransferTimer,
        EMTT_LastTimer
    };

protected:
    /// Array of Timer classes.
    Timer timers [ EMTT_LastTimer ];

public:
    EngngModelTimer() { }
    ~EngngModelTimer() { }

    /**@name Profiling routines. */
    //@{
    void startTimer(EngngModelTimerType t) { timers [ t ].startTimer(); }
    void stopTimer(EngngModelTimerType t) { timers [ t ].stopTimer(); }
    void pauseTimer(EngngModelTimerType t) { timers [ t ].pauseTimer(); }
    void resumeTimer(EngngModelTimerType t) { timers [ t ].resumeTimer(); }
    void initTimer(EngngModelTimerType t) { timers [ t ].initTimer(); }
    //@}

    /**@name Reporting routines. */
    //@{
    /// Returns total user time elapsed.
    double getUtime(EngngModelTimerType t);
    /// Returns elapsed wall clock time.
    double getWtime(EngngModelTimerType t);
    /// Returns pointer to timer determined by EngngModelTimerType.
    const Timer *getTimer(EngngModelTimerType t)  { return timers + t; }
    /// Converts total seconds into hours, mins, and seconds.
    static void convert2HMS(int &nhrs, int &nmin, int &nsec, double tsec);
    /// Printing & formatting.
    void toString(EngngModelTimerType t, char *buff);
    //@}
};
} // end namespace oofem
#endif // timer_h
