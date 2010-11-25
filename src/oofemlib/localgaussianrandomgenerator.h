/* $Header: $ */
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

/*
 * Author: Peter Grassl
 */

// file: localrandomgenerator.h

// **********************************************************
// CLASS LOCAL RANDOM GENERATOR
// **********************************************************


#ifndef localrandomgenerator_h
#define localrandomgenerator_h

#include "randomfieldgenerator.h"
#include "gausspnt.h"

namespace oofem {
/**
 *  This class implements a local (no spatial correlation) random generator using Guassian distribution.
 */
class LocalGaussianRandomGenerator : public RandomFieldGenerator
{
protected:
    /// integer which is the input of the pseudo-random number generator
    long randomInteger;
    /// gauss distribution parameters
    double mean, variance;
public:

    /// Constructor. Creates empty RandomFieldGenerator
    LocalGaussianRandomGenerator(int n, Domain *d);

    /// Destructor
    virtual ~LocalGaussianRandomGenerator();

    /**
     * Computes the random value.
     */
    void generateRandomValue(double &value, FloatArray *position);
    void generateRandomValueAt(double &value, GaussPoint *gp) {
        this->generateRandomValue(value, NULL);
    }

    virtual IRResultType initializeFrom(InputRecord *ir);
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const { return "LocalGaussianRandomGenerator"; }


protected:
    /**
     * Computes pseudo-random numbers.
     * @param idum Pointer to start integer (must be negative)
     * @return Random number between 0 and 1
     */
    double ran1(long *idum);

    /**
     * Computes the inverse of the Gaussian CDF
     * @param x Input probability
     * @param a Mean
     * @param b Standard deviation
     * @returns Inverse
     */
    double normalCdfInverse(double cdf, double a, double b);

    /**
     * Computes the inverse of the normal distribution
     * @param p Input probability
     * @returns Inverse
     */
    double normal01CdfInverse(double p);

    double dpolyValue(int n, double a[], double x);
};
} // end namespace oofem
#endif

