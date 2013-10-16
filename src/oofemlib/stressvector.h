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

#ifndef stressvector_h
#define stressvector_h

#include "stressstrainbasevector.h"

namespace oofem {
class StrainVector;
class FloatArray;
class FloatMatrix;

/**
 * Specialization of a floating point array for representing a stress state.
 */
class OOFEM_EXPORT StressVector : public StressStrainBaseVector
{
public:
    /// Constructor. Creates zero value stress/strain vector for given material mode.
    StressVector(MaterialMode);
    /// Constructor. Creates stress/strain vector, values taken from given vector, mode is parameter.
    StressVector(const FloatArray &, MaterialMode);
    /// Destructor.
    ~StressVector() { }
    /**
     * Member function that computes principal values of receiver (stress vector).
     * @param answer computed principal values (sorted)
     */
    void computePrincipalValues(FloatArray &answer) const;
    /**
     * Computes principal values and directions of receiver vector.
     * @param answer Computed principal values (sorted).
     * @param dir Principal directions (stored column wise).
     */
    void computePrincipalValDir(FloatArray &answer, FloatMatrix &dir) const;
    /**
     * Computes split of receiver into deviatoric and volumetric part.
     * @param dev Computed deviatoric part.
     * @param vol Computed volumetric part.
     */
    void computeDeviatoricVolumetricSplit(StressVector &dev, double &vol) const;

    /**
     * Sums volumetric part to receiver.
     * @param answer Computed full stress vector.
     * @param vol Volumetric part to add.
     */
    void computeDeviatoricVolumetricSum(StressVector &answer, double vol) const;
    /// Prints receiver on stdout, useful for debugging.
    void printYourself() const;
    /**
     * Computes the first invariant I1 of the stress.
     * @return First invariant.
     */
    double computeFirstInvariant() const;
    /**
     * Computes the second invariant J2 of the deviatoric stress.
     * @return Second invariant.
     */
    double computeSecondInvariant() const;
    /**
     * Computes the third invariant J3 of the deviatoric stress state.
     * @return Third invariant.
     */
    double computeThirdInvariant() const;

    /**
     * Computes all three Haigh-Westergaard coordinate of the stress.
     * @param xsi First HW-coordinate.
     * @param rho Second HW-coordinate.
     * @param theta third HW-coordinate.
     */
    void computeAllThreeHWCoordinates(double &xsi, double &rho, double &theta) const;

    /**
     * Computes the first Haigh-Westergaard coordinate of the  stress.
     * @return @f$ \displaystyle \xi = \frac{I_1}{\sqrt{3}} @f$.
     */
    double computeFirstCoordinate() const;
    /**
     * Computes the second Haigh-Westergaard coordinate of the deviatoric stress.
     * @return @f$ \rho = \sqrt{2 J_2} @f$.
     */
    double computeSecondCoordinate() const;
    /**
     * Computes the third Haigh-Westergaard coordinate of the deviatoric stress.
     * @return @f$ \displaystyle \theta = 3\frac{\sqrt{3}}{2} \frac{J_3}{J_2^{\frac 32}} @f$.
     */
    double computeThirdCoordinate() const;

    /**
     * Applies the isotropic elastic compliance to the stress.
     * @param strain Computed strain.
     * @param EModulus Elasticity modulus of the material.
     * @param nu Poisson's ratio of the material.
     */
    void applyElasticCompliance(StrainVector &strain,
                                const double EModulus,
                                const double nu) const;
    /**
     * Applies the isotropic elastic stiffness to the deviatoric stress.
     * @param strain Computed strain
     * @param EModulus Elasticity modulus of the material.
     * @param nu Poisson's ratio of the material.
     */
    void applyDeviatoricElasticCompliance(StrainVector &strain,
                                          const double EModulus,
                                          const double nu) const;
    /**
     * Applies the isotropic elastic stiffness to the deviatoric stress.
     * @param strain Computed strain.
     * @param GModulus Shear modulus of the material.
     */
    void applyDeviatoricElasticCompliance(StrainVector &strain,
                                          const double GModulus) const;
    /**
     * Computes the norm of the stress tensor using engineering notation.
     * @return Norm of the deviatoric part.
     */
    double computeStressNorm() const;

protected:
    void giveTranformationMtrx(FloatMatrix &answer, const FloatMatrix &base,
                               int transpose = 0) const;
};
} // end namespace oofem
#endif // stressvector_h
