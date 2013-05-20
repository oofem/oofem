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

#ifndef localrandomgenerator_h
#define localrandomgenerator_h

#include "randomfieldgenerator.h"
#include "gausspoint.h"

///@name Input fields for LocalGaussianRandomGenerator
//@{
#define _IFT_LocalGaussianRandomGenerator_Name "localgaussrandomgenerator"
#define _IFT_LocalGaussianRandomGenerator_mean "mean"
#define _IFT_LocalGaussianRandomGenerator_variance "variance"
#define _IFT_LocalGaussianRandomGenerator_seed "seed"
//@}

namespace oofem {
/**
 * This class implements a local (no spatial correlation) random generator using Gaussian distribution.
 * @author Peter Grassl
 */
class LocalGaussianRandomGenerator : public RandomFieldGenerator
{
protected:
    /// Integer which is the input of the pseudo-random number generator.
    long randomInteger;
    /// Gauss distribution parameters.
    double mean, variance;

public:
    /// Constructor. Creates empty RandomFieldGenerator
    LocalGaussianRandomGenerator(int n, Domain *d);
    /// Destructor
    virtual ~LocalGaussianRandomGenerator();

    void generateRandomValue(double &value, FloatArray *position);
    void generateRandomValueAt(double &value, GaussPoint *gp) {
        this->generateRandomValue(value, NULL);
    }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "LocalGaussianRandomGenerator"; }
    virtual const char *giveInputRecordName() const { return _IFT_LocalGaussianRandomGenerator_Name; }
    virtual classType giveClassID() const { return LocalGaussianRandomGeneratorClass; }

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

