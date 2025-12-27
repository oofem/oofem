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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#pragma once

#include "boundarycondition.h"
#include "dof.h"
#include "bctype.h"
#include "valuemodetype.h"
#include "floatarray.h"
#include "floatmatrix.h"

#include <map>

///@name Input fields for TransportGradientDirichlet
//@{
#define _IFT_TransportGradientDirichlet_Name "tmgraddirichlet"
#define _IFT_TransportGradientDirichlet_gradient "gradient"
#define _IFT_TransportGradientDirichlet_centerCoords "centercoords"
#define _IFT_TransportGradientDirichlet_surfSets "surfsets"
#define _IFT_TransportGradientDirichlet_edgeSets "edgesets"
#define _IFT_TransportGradientDirichlet_tractionControl "usexi"
//@}

namespace oofem {

/**
 * Prescribes @f$ T = g_{i}(x_i-\bar{x}_i) @f$ where @f$ T @f$ are primary unknown.
 * This is typical boundary condition in multiscale analysis where @f$ g = \partial_x T @f$
 * would be a macroscopic gradient at the integration point, i.e. this is a boundary condition for prolongation.
 * It is also convenient to use when one wants to test a arbitrary specimen for a given average gradient.
 * 
 * @note This class has optional support for the "psi" parameter, which is computed to improve the 
 * bounds for high permeable inclusions in a matrix.
 * Using this parameter requires additional inputs, edgeSets and surfSets, which should be defined in the particular order.
 * The numbering of the surface should be the same as for the hex cube interpolator as defined in the element reference manual.
 * 
 * @note Experimental code
 * 
 * @author Mikael Öhman
 */
class OOFEM_EXPORT TransportGradientDirichlet : public BoundaryCondition
{
protected:
    FloatArray mGradient;
    FloatArray mCenterCoord;

    bool tractionControl = false;
    /// Stores one "psi" value for each node
    std :: map< int, FloatArray > xis;
    IntArray surfSets;

public:
    TransportGradientDirichlet(int n, Domain *d) : BoundaryCondition(n, d) { }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    void postInitialize() override;

    double give(Dof *dof, ValueModeType mode, double time) override;

    bcType giveType() const override { return DirichletBT; }

    double domainSize();
    /**
     * Constructs a coefficient matrix for all prescribed unknowns.
     * Helper routine for computational homogenization.
     * @todo Perhaps this routine should only give results for the dof it prescribes?
     * @param C Coefficient matrix to fill.
     */
    void computeCoefficientMatrix(FloatMatrix &C);

    /**
     * Computes the homogenized, macroscopic, field (stress).
     * @param sigma Output quantity (typically stress).
     * @param tStep Active time step.
     */
    virtual void computeField(FloatArray &sigma, TimeStep *tStep);

    /**
     * Computes the macroscopic tangent for homogenization problems through sensitivity analysis.
     * @param tangent Output tangent.
     * @param tStep Active time step.
     */
    virtual void computeTangent(FloatMatrix &tangent, TimeStep *tStep);

    /// Computes the offset values for "improved" Dirichlet. See class description.
    void computeXi();

    void scale(double s) override { mGradient.times(s); }

    const char *giveClassName() const override { return "TransportGradientDirichlet"; }
    const char *giveInputRecordName() const override { return _IFT_TransportGradientDirichlet_Name; }
};
} // end namespace oofem
