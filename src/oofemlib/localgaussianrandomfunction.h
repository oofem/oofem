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

#ifndef localgaussianrandomfunction_h
#define localgaussianrandomfunction_h

#include "function.h"

///@name Input fields for LocalGaussianRandomFunction
//@{
#define _IFT_LocalGaussianRandomFunction_Name "localgaussrandomfunction"
#define _IFT_LocalGaussianRandomFunction_mean "mean"
#define _IFT_LocalGaussianRandomFunction_variance "variance"
#define _IFT_LocalGaussianRandomFunction_seed "seed"
//@}

namespace oofem {
/**
 * This class implements a local (no spatial correlation) random function using Gaussian distribution.
 * @author Peter Grassl
 */
class OOFEM_EXPORT LocalGaussianRandomFunction : public Function
{
protected:
    /// Integer which is the input of the pseudo-random number generator.
    long randomInteger;
    /// Gauss distribution parameters.
    double mean, variance;

public:
    /// Constructor.
    LocalGaussianRandomFunction(int n, Domain * d);
    /// Destructor
    virtual ~LocalGaussianRandomFunction();

    virtual void evaluate(FloatArray &answer, const std :: map< std :: string, FunctionArgument > &valDict);
    virtual double evaluateAtTime(double t);
    virtual double evaluateVelocityAtTime(double t);
    virtual double evaluateAccelerationAtTime(double t);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "LocalGaussianRandomFunction"; }
    virtual const char *giveInputRecordName() const { return _IFT_LocalGaussianRandomFunction_Name; }

protected:
    /**
     * Computes pseudo-random numbers.
     * @param idum Pointer to start integer (must be negative).
     * @return Random number between 0 and 1.
     */
    double ran1(long *idum);

    /**
     * Computes the inverse of the Gaussian CDF
     * @param cdf Input probability.
     * @param a Mean.
     * @param b Standard deviation.
     * @return Inverse.
     */
    double normalCdfInverse(double cdf, double a, double b);

    /**
     * Computes the inverse of the normal distribution.
     * @param p Input probability.
     * @return Inverse.
     */
    double normal01CdfInverse(double p);

    double dpolyValue(int n, double a[], double x);
};
} // end namespace oofem
#endif
