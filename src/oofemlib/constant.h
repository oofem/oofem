/* $Header: /home/cvs/bp/oofem/oofemlib/src/constant.h,v 1.8 2003/04/06 14:08:23 bp Exp $ */
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

//   *******************************
//   *** CLASS CONSTANT FUNCTION ***
//   *******************************

#ifndef constant_h
#define constant_h

#include "loadtime.h"

/**
 * Class implementing time function y=f(t) that is constant in time.
 */
class ConstantFunction : public LoadTimeFunction
{
    /*
     * This class implement a function  y = f(t)  that is constant in time.
     * DESCRIPTION
     * 'value' is the constant value of the function. It is a pointer, rather
     * than a number, so that its state (initialized or not) can be checked.
     */
private:
    /// Value of receiver
    double value;

public:
    /**
     * Constructor. Creates constant load time function with given number, belonging to given domain.
     * @param n load time function number
     * @param d domain to which new object will belongs.
     */
    ConstantFunction(int i, Domain *d) : LoadTimeFunction(i, d) { value = 0; }
    /// Destructor.
    ~ConstantFunction()                       { }

    /**
     * Returns value member of receiver.
     */
    double  giveValue();
    /**
     * Initializes receiver acording to object description stored in input record.
     */
    IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     *  @param str string to be filled by input record
     *  @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    /**
     * Returns classType id of receiver.
     * @return ConstantFunctionClass value.
     */
    classType   giveClassID() const { return ConstantFunctionClass; }
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "ConstantFunction"; }
    /// Returns input record name of the receiver.
    const char *giveInputRecordName() const { return "ConstantFunction"; }

    /**
     * Returns the value of load time function at given time. Abstract service.
     * Must be implemented by derived classes.
     * @param t time
     * @return load time function value
     */
    virtual double     __at(double)            { return this->giveValue(); }
};

#endif // constant_h
