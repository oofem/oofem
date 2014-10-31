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

#ifndef stressstrainbasevector_h
#define stressstrainbasevector_h

#include "floatarray.h"
#include "materialmode.h"

namespace oofem {
typedef MaterialMode StressStrainMatMode;
/**
 * Base class for stress/strain vector representations. It is derived from FloatArray class, which is used
 * to store the stress/strain vector components in reduced form. Additional attribute is introduced
 * to keep the stress/strain mode. This allows to implement various stress/strain vector methods
 * regarding the correct stress/strain mode without doing so in full form.
 * The full form includes all components, even if they are zero due to stress/strain mode nature,
 * but in the reduced format, only relevant nonzero components are stored.
 * If in particular mode particular stress component is zero, the corresponding strain is not computed
 * and not stored in reduced vector, and in full vector there is zero value on corresponding position.
 * On the other hand, if zero strain component is imposed this condition must be taken into account in geometrical
 * relations (at element level), and corresponding component are included stress/strain reduced vectors.
 * Methods for converting vectors between reduced and full format are provided.
 */
class OOFEM_EXPORT StressStrainBaseVector : public FloatArray
{
protected:
    /// Stress strain mode.
    StressStrainMatMode mode;

public:
    /// Constructor. Creates zero value stress/strain vector for given material mode.
    StressStrainBaseVector(MaterialMode);
    /// Constructor. Creates stress/strain vector, values taken from given vector, mode is parameter.
    StressStrainBaseVector(const FloatArray &, MaterialMode);
    /// Destructor
    ~StressStrainBaseVector() { }

    /** Assignment operator. Defines the assignment between two Stress Strain Base Vectors.
     *  The Assignment between FloatArray and StressStrainBaseVector is also possible, in this case
     *  the FloatArray assignment operator is used, where the operand is converted to FloatArray (mode is lost).
     *  The assignment from FloatArray to StressStrainBaseVector is not defined (needs material mode).
     */
    StressStrainBaseVector &operator = ( const StressStrainBaseVector & );

    /// Returns the material mode of receiver.
    MaterialMode giveStressStrainMode() const { return ( MaterialMode ) mode; }

    /**
     * Changes the material mode of receiver.
     * @param newMode New mode.
     */
    void letStressStrainModeBe(const MaterialMode newMode);

    /**
     * Computes the full form of receiver.
     * @param fullform Requested full form of receiver.
     */
    void convertToFullForm(FloatArray &fullform) const;
    /**
     * Assign to receiver the reduced form of given vector.
     * @param reducedform Reduced form of receiver to expand.
     * @param mode Mode of the stress/strain.
     */
    void convertFromFullForm(const FloatArray &reducedform, MaterialMode mode);

    /**
     * Member function that computes principal values of receiver (strain vector).
     * @param answer Computed principal values (sorted).
     */
    virtual void computePrincipalValues(FloatArray &answer) const = 0;
    /**
     * Computes principal values and directions of receiver vector.
     * @param answer Computed principal values (sorted).
     * @param dir Principal directions (stored column wise).
     */
    virtual void computePrincipalValDir(FloatArray &answer, FloatMatrix &dir) const = 0;
    /**
     * Transforms receiver vector into another coordinate system.
     * @param answer transformed strain vector
     * @param base Transformation matrix.  There are on each column stored unit vectors of
     * coordinate system (so called base vectors) to which we do transformation. These vectors must
     * be expressed in the same coordinate system as source strainVector.
     * @param answer Transformed 3d strain.
     * @param base New base to express vector in.
     * @param transpose If transpose == 1 then transpose base matrix before transforming.
     */
    void transformTo(StressStrainBaseVector &answer, const FloatMatrix &base, int transpose = 0) const;

    contextIOResultType storeYourself(DataStream &stream);
    contextIOResultType restoreYourself(DataStream &stream);

    /// Returns the volumetric part of the vector.
    double computeVolumetricPart() const;

protected:
    /**
     * Computes 3d transformation matrix from standard vector transformation matrix.
     * @param answer transformation matrix for strain vector
     * @param base (3,3) matrix, where on each column are stored unit direction vectors of
     * local coordinate axes to which we do transformation.
     * @param transpose If transpose == 1 then transpose base matrix before transforming
     */
    virtual void giveTranformationMtrx(FloatMatrix &answer, const FloatMatrix &base,
                                       int transpose = 0) const = 0;
};
} // end namespace oofem
#endif // stressstrainbasevector_h
