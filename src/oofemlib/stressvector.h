/* $Header: /home/cvs/bp/oofem/oofemlib/src/stressvector.h,v 1.1.4.1 2004/04/05 15:19:44 bp Exp $ */
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

#ifndef stressvector_h
#define stressvector_h

#include "stressstrainbasevector.h"

namespace oofem {

class StrainVector;
class FloatArray;
class FloatMatrix;

/**
 */
class StressVector : public StressStrainBaseVector
{
protected:
public:
    /// Constructor. Creates zero value stress/strain vector for given material mode.
    StressVector(MaterialMode);
    /// Constructor. Creates stress/strain vector, values taken from given vector, mode is parameter.
    StressVector(const FloatArray &, MaterialMode);
    /// Destructor
    ~StressVector() { }
    /**
     *      Member function that computes principal values of receiver (stress vector).
     *      @param answer computed principal values (sorted)
     *      @param s stress/strain vector which eigenvalues are computed
     */
    void computePrincipalValues(FloatArray &answer) const;
    /**
     * Computes principal values and directions of receiver vector.
     * @param answer computed principal values (sorted)
     * @param dir principal directions (stored columwise)
     * @param s stress/strain vector
     */
    void computePrincipalValDir(FloatArray &answer, FloatMatrix &dir) const;
    /**
     * Computes split of receiver into deviatoric and volumetric part
     */
    void computeDeviatoricVolumetricSplit(StressVector &dev, double &vol) const;

    /**
     * Computes sum of deviatoric and volumetric part.
     */
    void  computeDeviatoricVolumetricSum(StressVector &answer,
                                         double &vol) const;
    /**
     * Prints receiver on stdout, usefull for debugging
     */
    void printYourself() const;
    /**
     * Computes the first invariant I1 of the stress.
     */
    double computeFirstInvariant() const;
    /**
     * Computes the second invariant J2 of the deviatoric stress.
     */
    double computeSecondInvariant() const;
    /**
     * Computes the third invariant J3 of the deviatoric stress state s.
     */
    double computeThirdInvariant() const;

    /**
     * Computes all three Haigh-Westergaard coordinate of the stress.
     * @param xsi first HW-coordinate
     *       @param rho second HW-coordinate
     *       @param theta third HW-coordinate
     */
    void computeAllThreeHWCoordinates(double &xsi, double &rho, double &theta) const;

    /**
     * Computes the first Haigh-Westergaard coordinate of the  stress.
     * xsi = I1 / sqrt(3.)
     */
    double computeFirstCoordinate() const;
    /**
     * Computes the second Haigh-Westergaard coordinate of the deviatoric stress.
     * rho = sqrt(2.*J2)
     */
    double computeSecondCoordinate() const;
    /**
     * Computes the third Haigh-Westergaard coordinate of the deviatoric stress.
     * theta = 3.*sqrt(3.)/2.*J3/pow(J2,3./2.)
     */
    double computeThirdCoordinate() const;

    /**
     * Applies the elastic compliance to the stress.
     * @param strain computed strain
     * @param Emodulus Emodulus of the material
     * @param nu Poisson's ratio of the material
     */
    void applyElasticCompliance(StrainVector &strain,
                                const double EModulus,
                                const double nu) const;
    /**
     * Applies the elastic stiffness to the deviatoric stress.
     * @param strain computed strain
     * @param Emodulus Emodulus of the material
     * @param nu Poisson's ratio of the material
     */
    void applyDeviatoricElasticCompliance(StrainVector &strain,
                                          const double EModulus,
                                          const double nu) const;

    /**
     * Applies the elastic stiffness to the deviatoric stress.
     * @param strain computed strain
     * @param GModulus Gmodulus of the material
     */
    void applyDeviatoricElasticCompliance(StrainVector &strain,
                                          const double GModulus) const;
    /**
     * Computes the norm of the stress tensor using engineering notation.
     */
    double computeStressNorm() const;
protected:

    void giveTranformationMtrx(FloatMatrix &answer, const FloatMatrix &base,
                               int transpose = 0) const;
};

} // end namespace oofem
#endif // stressvector_h
