/* $Header: /home/cvs/bp/oofem/sm/src/piecewis.h,v 1.4 2003/04/06 14:08:31 bp Exp $ */
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

//   ***************************************
//   *** CLASS PIECEWISE LINEAR FUNCTION ***
//   ***************************************


#ifndef piecewis_h
#define piecewis_h



#include "flotarry.h"
#include "loadtime.h"

class PiecewiseLinFunction : public LoadTimeFunction
{
    /*
     * This class implements a piecewise linear function.
     * DESCRIPTION
     * The function is defined by 'numberOfPoints' points. 'dates' and 'values'
     * store respectively the abscissas (t) and the values (f(t)) of the points
     */

protected:
    int numberOfPoints;
    FloatArray dates;
    FloatArray values;

public:
    PiecewiseLinFunction(int i, Domain *d) : LoadTimeFunction(i, d), dates(), values()
    { numberOfPoints = 0; }
    ~PiecewiseLinFunction()             { }

    //      void    getPoints () ;
    IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     * @param str string to be filled by input record
     * @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    classType    giveClassID() const { return PiecewiceClass; }
    const char *giveClassName() const { return "PiecewiceClass"; }
    /// Returns input record name of the receiver.
    const char *giveInputRecordName() const { return "PiecewiseLinFunction"; }

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
    /**
     * Returns the second time derivative of load time function at given time. Abstract service.
     * Must be implemented by derived classes.
     * @param t time
     * @return load time function value
     */
    virtual double    __accelAt(double) { _error("accelAt: not supported");
                                          return 0.; }
};

#endif // piecewis_h

