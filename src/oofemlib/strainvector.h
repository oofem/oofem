/* $Header: /home/cvs/bp/oofem/oofemlib/src/strainvector.h,v 1.1.4.1 2004/04/05 15:19:44 bp Exp $ */
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

#ifndef strainvector_h
#define strainvector_h

#include "stressstrainbasevector.h"

namespace oofem {

class StressVector;
/**
 */
class StrainVector : public StressStrainBaseVector
{
protected:
public:
    /// Constructor. Creates zero value stress/strain vector for given material mode.
    StrainVector(MaterialMode);
    /// Constructor. Creates stress/strain vector, values taken from given vector, mode is parameter.
    StrainVector(const FloatArray &, MaterialMode);
    /// Destructor
    ~StrainVector() { }

    /**
     * Member function that computes principal values of receiver (strain vector).
     * @param answer computed principal values (sorted)
     * @param s stress/strain vector which eigenvalues are computed
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
     * @param answer deviatoric strain
     * @param vol volumetric strain
     */
    void computeDeviatoricVolumetricSplit(StrainVector &dev, double &vol) const;
    /**
     * Computes sum of deviatoric and volumetric part.
     * @param answer total strain
     *        @param vol volumetric strain
     */
    void  computeDeviatoricVolumetricSum(StrainVector &answer, const double vol) const;
    /**
     *      Prints receiver on stdout, usefull for debugging
     */
    void printYourself() const;
    /**
     * Computes the change of volume. epsv = eps1 + eps2 + eps3
     */
    double computeVolumeChange() const;
    /**
     * Computes the tensorial norm of the strain in engineering notation.
     */
    double computeStrainNorm() const;
    /**
     *        Applies the elastic stiffness to the strain.
     * @param stress computed stress
     * @param Emodulus Emodulus of the material
     * @param nu Poisson's ratio of the material
     */
    void applyElasticStiffness(StressVector &stress, const double EModulus, const double nu) const;
    /**
     *        Applies the elastic stiffness to the deviatoric strain.
     * @param stress computed stress
     * @param Emodulus Emodulus of the material
     * @param nu Poisson's ratio of the material
     */
    void applyDeviatoricElasticStiffness(StressVector &stress, const double EModulus, const double nu) const;
    /**
     * Applies the elastic stiffness to the deviatoric strain.
     * @param stress computed stress
     * @param Gmodulus Gmodulus of the material
     */
    void applyDeviatoricElasticStiffness(StressVector &stress, const double GModulus) const;

protected:

    void giveTranformationMtrx(FloatMatrix &answer, const FloatMatrix &base,
                               int transpose = 0) const;
};

} // end namespace oofem
#endif // strainvector_h
