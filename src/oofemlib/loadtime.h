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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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
/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#ifndef loadtime_h
#define loadtime_h

#include "femcmpnn.h"
#include "domain.h"
#include "classtype.h"
#include "valuemodetype.h"

namespace oofem {
/**
 * Abstract base class representing load time function. Classes derived from Load class typically
 * describe load from spatial point of view. The purpose of introducing load time function is to express
 * variation of some components in time. Load time function typically belongs to domain and is
 * attribute of one or more loads. Generally load time function is real function of time, @f$y=f(t)@f$.
 *
 * See TJR Hughes, "The Finite Element Method", p 677.
 */
class LoadTimeFunction : public FEMComponent
{
protected:
    /**
     * By default, the increment of receiver is computed as a difference between values evaluated at given solution step and in previous step.
     * However, if the solution step is the first step, the difference is typically set to the total value of receiver at the first step.
     * This is quite natural, as a loading with constant time function is expected to be applied at first step.
     * In certain cases, this default behavior has to be changed. The initial value (set by default to zero)
     * allows to set initial value of receiver. This initial value is used only when the increment of receiver is evaluated at first step,
     * when result is defined as value of receiver at given step minus the initial value.
     * This allows to correctly handle temperature loading, that is specified with respect to some reference temperature.
     * In this case, the initial value should be set to the reference temperature, allowing to obtain correct temperate increment in first step.
     */
    double initialValue;
public:

    /**
     * Constructor. Creates load time function with given number, belonging to given domain.
     * @param n Load time function number.
     * @param d Domain to which new object will belongs.
     */
    LoadTimeFunction(int n, Domain *d) : FEMComponent(n, d) { initialValue = 0.0; }
    /// Destructor
    virtual ~LoadTimeFunction() { }

    /**
     * Returns the value of load time function at given time. Abstract service.
     * Must be implemented by derived classes.
     * @param tStep Time. Incremental and total mode uses intrinsic time from giveIntrinsicTime.
     * @param mode Determines the mode of the requested value.
     * @return Load time function value.
     */
    double evaluate(TimeStep *tStep, ValueModeType mode);

    /**
     * Returns the value of load time function at given time.
     * @param t Time.
     * @return @f$ f(t) @f$.
     */
    virtual double  __at(double t) { return 0.; }
    /**
     * Returns the first time derivative of load time function at given time.
     * @param t Time.
     * @return @f$ f'(t) @f$.
     */
    virtual double __derAt(double t) { return 0.; }
    /**
     * Returns the second time derivative of load time function at given time.
     * @param t Time.
     * @return @f$ f''(t) @f$.
     */
    virtual double __accelAt(double t) { return 0.; }

    // Overloaded methods:
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual IRResultType initializeFrom(InputRecord *ir, DataReader *dr);
    
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    virtual classType giveClassID() const { return LoadTimeFunctionClass; }
    virtual const char *giveClassName() const { return "LoadTimeFunction"; }
};
} // end namespace oofem
#endif // loadtime_h
