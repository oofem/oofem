/* $Header: /home/cvs/bp/oofem/sm/src/piecewisper.h,v 1.4 2003/04/06 14:08:31 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

//   ************************************************
//   *** CLASS PERIODIC PIECEWISE LINEAR FUNCTION ***
//   ************************************************
#ifndef piecewisper_h
#define piecewisper_h

#include "flotarry.h"
#include "piecewis.h"

namespace oofem {
class PeriodicPiecewiseLinFunction : public PiecewiseLinFunction
{
    /*
     * This class implements an enhanced piecewise linear function with periodicity
     * and possibility to add another arbitrary time function.
     * DESCRIPTION
     * The function is defined by 'numberOfPoints' points. 'dates' and 'values'
     * store respectively the abscissas (t) and the values (f(t)) of the points
     * The values are repeated after 'period'. 'AddTF' parameter specifies number
     * of function to add.
     */

private:
    /** If nonzero, the value of time function specified by addTF is added to computed value */
    int addTF;
    /** If less than zero no periodicity, if >=0 date time is computed as
    * given time%period.
    * If points span more than period, span of LAST period is repeated */
    double period;

public:
    PeriodicPiecewiseLinFunction(int i, Domain *d) : PiecewiseLinFunction(i, d)
    {
        period = -1.0;
        addTF = 0;
    }
    ~PeriodicPiecewiseLinFunction()             { }

    //      void    getPoints () ;
    IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     * @param str string to be filled by input record
     * @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    classType    giveClassID() const { return PeriodicPiecewiseClass; }
    const char *giveClassName() const { return "PeriodicPiecewiseClass"; }
    /// Returns input record name of the receiver.
    const char *giveInputRecordName() const { return "PeriodicPiecewiseLinFunction"; }
    /**
     * Returns the value of load time function at given time. Abstract service.
     * Must be implemented by derived classes.
     * @param t time
     * @return load time function value
     */
    virtual double     __at(double);
    /**
     * Returns the first time derivative of load time function at given time. Abstract service.
     * Must be implemented by derived classes.
     * @param t time
     * @return load time function value
     */
    virtual double    __derAt(double);
};
} // end namespace oofem
#endif // piecewisper_h
