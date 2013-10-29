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

#ifndef strainvector_h
#define strainvector_h

#include "stressstrainbasevector.h"

namespace oofem {
class StressVector;
/**
 * Specialization of a floating point array for representing a strain state.
 */
class OOFEM_EXPORT StrainVector : public StressStrainBaseVector
{
public:
    /// Constructor. Creates zero value stress/strain vector for given material mode.
    StrainVector(MaterialMode);
    /// Constructor. Creates stress/strain vector, values taken from given vector, mode is parameter.
    StrainVector(const FloatArray &, MaterialMode);
    /// Destructor
    ~StrainVector() { }

    /**
     * Computes the principal values of the receiver.
     * @param answer Principal values.
     */
    void computePrincipalValues(FloatArray &answer) const;
    /**
     * Computes the principal values and principal directions of the receiver.
     * The principal values are ordered from the largest to the smallest.
     * @param answer Principal values.
     * @param dir Principal directions (first index refers to component, second index to eigenvalue number).
     */
    void computePrincipalValDir(FloatArray &answer, FloatMatrix &dir) const;
    /**
     * Computes the principal direction of the receiver
     * associated with the maximum principal value.
     * @param answer Principal direction.
     */
    void computeMaxPrincipalDir(FloatArray &answer) const;
    /**
     * Computes split of receiver into deviatoric and volumetric part.
     * @param answer Deviatoric strain.
     * @param vol Volumetric strain.
     */
    void computeDeviatoricVolumetricSplit(StrainVector &answer, double &vol) const;
    /**
     * Computes sum of deviatoric and volumetric part.
     * @param answer Total strain.
     * @param vol Volumetric strain.
     */
    void computeDeviatoricVolumetricSum(StrainVector &answer, const double vol) const;
    /**
     * Prints receiver on stdout, useful for debugging
     */
    void printYourself() const;
    /**
     * Computes the change of volume.
     * @return @f$ \epsilon_v = \epsilon_1 + \epsilon_2 + \epsilon_3 @f$.
     */
    double computeVolumeChange() const;
    /**
     * Computes the tensorial norm of the strain in engineering notation.
     */
    double computeStrainNorm() const;
    /**
     * Applies the elastic stiffness to the strain.
     * @param stress Computed stress.
     * @param EModulus Elasticity modulus of the material.
     * @param nu Poisson's ratio of the material
     */
    void applyElasticStiffness(StressVector &stress, const double EModulus, const double nu) const;
    /**
     * Applies the elastic stiffness to the deviatoric strain.
     * @param stress Computed stress.
     * @param EModulus Elasticity modulus of the material.
     * @param nu Poisson's ratio of the material.
     */
    void applyDeviatoricElasticStiffness(StressVector &stress, const double EModulus, const double nu) const;
    /**
     * Applies the elastic stiffness to the deviatoric strain.
     * @param stress Computed stress.
     * @param GModulus Shear modulus of the material
     */
    void applyDeviatoricElasticStiffness(StressVector &stress, const double GModulus) const;

protected:

    void giveTranformationMtrx(FloatMatrix &answer, const FloatMatrix &base,
                               int transpose = 0) const;
};
} // end namespace oofem
#endif // strainvector_h
