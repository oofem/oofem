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

#ifndef prescribedgradienthomogenization_h
#define prescribedgradienthomogenization_h

#include "inputrecord.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for PrescribedGradientHomogenization
//@{
#define _IFT_PrescribedGradientHomogenization_Name "prescribedgradient"
#define _IFT_PrescribedGradientHomogenization_centercoords "ccoord"
#define _IFT_PrescribedGradientHomogenization_gradient "gradient"
//@}

namespace oofem {
class TimeStep;
class DynamicInputRecord;
class Domain;

/**
 * Class for homogenization of applied gradients.
 * This is typically applied to a boundary condition in multiscale analysis where @f$ d = \partial_x s@f$
 * 
 * @author Mikael Öhman
 */
class OOFEM_EXPORT PrescribedGradientHomogenization
{
protected:
    /// Prescribed gradient @f$ d_{ij} @f$
    FloatMatrix mGradient;

    /// Center coordinate @f$ \bar{x}_i @f$
    FloatArray mCenterCoord;

    virtual double domainSize(Domain *d, int set);

public:
    PrescribedGradientHomogenization() { }
    virtual ~PrescribedGradientHomogenization() { }

    /**
     * Initializes receiver according to object description stored in input record.
     * The input record contains two fields;
     * - gradient \#rows \#columns { d_11 d_12 ... ; d_21 ... } (required)
     * - cCoords \#columns x_1 x_2 ... (optional, default 0)
     * The prescribed gradients columns must be equal to the size of the center coordinates.
     * The size of the center coordinates must be equal to the size of the coordinates in the applied nodes.
     */
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    /**
     * Computes the homogenized, macroscopic, field (stress).
     * @param sigma Output quantity (typically stress).
     * @param tStep Active time step.
     */
    virtual void computeField(FloatArray &sigma, TimeStep *tStep) = 0;

    /**
     * Computes the macroscopic tangent for homogenization problems through sensitivity analysis.
     * @param tangent Output tangent.
     * @param tStep Active time step.
     */
    virtual void computeTangent(FloatMatrix &tangent, TimeStep *tStep) = 0;

    /**
     * Set prescribed gradient.
     * @param t New prescribed gradient.
     */
    void setPrescribedGradient(const FloatMatrix &t) { mGradient = t; }

    /**
     * Sets the prescribed gradient from the matrix from given voigt notation.
     * Assumes use of double values for off-diagonal, usually the way for strain in Voigt form.
     * @param t Vector in voigt format.
     */
    void setPrescribedGradientVoigt(const FloatArray &t);
    /**
     * Gives back the applied gradient in Voigt form.
     * @param oGradient The applied gradient, in Voigt form.
     */
    void giveGradientVoigt(FloatArray &oGradient) const;

    /**
     * Set the center coordinate for the prescribed values to be set for.
     * @param x Center coordinate.
     */
    void setCenterCoordinate(FloatArray &x) { mCenterCoord = x; }
    /// Returns the center coordinate
    FloatArray &giveCenterCoordinate() { return mCenterCoord; }

    virtual const char *giveClassName() const = 0;
};
} // end namespace oofem

#endif // prescribedgradienthomogenization_h
