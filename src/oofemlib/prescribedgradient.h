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
 *               Copyright (C) 1993 - 2010   Borek Patzak
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

#ifndef prescribedgradient_h
#define prescribedgradient_h

#include "boundary.h"
#include "dof.h"
#include "bctype.h"
#include "valuemodetype.h"
#include "classtype.h"
#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
/**
 * Prescribes @f$ v_i = d_{ij}(x_j-\bar{x}_j) @f$ or @f$ s = d_{1j}(x_j - \bar{x}_j) @f$
 * where @f$ v_i @f$ are primary unknowns for the subscale.
 * This is typical boundary condition in multiscale analysis where @f$ d = \partial_x s@f$
 * would a macroscopic gradient at the integration point, i.e. this is a boundary condition for prolongation.
 * It is also convenient to use when one wants to test a arbitrary specimen for shear.
 * @date 2010-04-22
 * @author Mikael Öhman
 */
class PrescribedGradient : public BoundaryCondition
{
protected:
    /// Prescribed gradient @f$ d_{ij} @f$
    FloatMatrix gradient;

    /// Center coordinate @f$ \bar{x}_i @f$
    FloatArray centerCoord;

public:
    /** Creates boundary condition with given number, belonging to given domain.
     * @param n boundary condition number
     * @param d domain to which new object will belongs
     */
    PrescribedGradient(int i, Domain *d) : BoundaryCondition(i, d) { }

    /// Destructor
    ~PrescribedGradient() { }

    /** Returns the value of a prescribed unknown, respecting requested mode for given time.
     * Its physical meaning is determined by corresponding DOF.
     * @param dof determines the dof subjected to receiver bc.
     * @param mode unknown char type (if total or incremental value is returned)
     * @return prescribed value of unknown or zero if not prescribed
     */
    virtual double give(Dof *dof, ValueModeType mode, TimeStep *tStep);

    /** Returns receiver load type. It distinguish particular boundary conditions according to
     * their "physical" meaning (like StructuralTemperatureLoadLT, StructuralLoadLT).
     * @return returns BoundaryConditionLT value.
     */
    bcType giveType() const { return DirichletBT; }

    /** Initializes receiver acording to object description stored in input record.
     * The input record contains two fields;
     * - gradient #rows #columns { d_11 d_12 ... ; d_21 ... } (required)
     * - cCoords #columns x_1 x_2 ... (optional, default 0)
     * The prescribed tensor's columns must be equal to the size of the center coordinates.
     * The size of the center coordinates must be equal to the size of the coordinates in the applied nodes.
     */
    IRResultType initializeFrom(InputRecord *ir);

    /** Setups the input record string of receiver.
     * keyword parameter is ignored as I don't know what it does (nor is it documented).
     * @param str string to be filled by input record
     * @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    /** Scales the receiver according to given value.
     * @param s scaling factor for entire tensor
     */
    virtual void scale(double s) { gradient.times(s); }

    /** Set prescribed tensor.
     *  @param s prescribed value
     */
    virtual void setPrescribedGradient(const FloatMatrix &t) { gradient = t; }

    /** Set prescribed value at Expresses the matrix from given voigt notation.
     * Assumes use of double values for off-diagonal, usually the way for strain in Voigt form.
     * @param t Vector in voigt format.
     */
    virtual void setPrescribedGradientVoigt(const FloatArray &t);

    /// @warning Not used. Do not call.
    virtual void setPrescribedValue(double) { OOFEM_ERROR("Scalar value not used for prescribed tensors."); }

    /** Set the center coordinate for the prescribed values to be set for.
     *  @param x center coordinate.
     */
    virtual void setCenterCoordinate(const FloatArray &x) { centerCoord = x; }

    /// Returns the center coordinate
    virtual FloatArray &giveCenterCoordinate() { return centerCoord; }

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "PrescribedGradient"; }

    /// Returns classType id of receiver.
    classType giveClassID() const { return PrescribedGradientClass; }
};
} // end namespace oofem

#endif // prescribedgradient_h

