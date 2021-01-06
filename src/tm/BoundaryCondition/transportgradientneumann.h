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

#pragma once

//#include "prescribedgradienthomogenization.h"
#include "activebc.h"
#include "floatarray.h"

#include <memory>

#define _IFT_TransportGradientNeumann_Name "tmgradneumann"
#define _IFT_TransportGradientNeumann_gradient "gradient"
#define _IFT_TransportGradientNeumann_centerCoords "centercoords"
#define _IFT_TransportGradientNeumann_surfSets "surfsets" /// Numered x+, y+, z+, x-, y-, z-
#define _IFT_TransportGradientNeumann_dispControl "useeta" /// For activating Voigt-Neumann b.c.

namespace oofem {
class Node;
class Element;

/**
 * Homogenization boundary condition that imposes a gradient weakly on the boundary with scaled Neumann boundary condition.
 * The scaling is determined such that the is preserved while retaining the statistical properties.
 * 
 * @note Experimental code
 * 
 * @author Mikael Öhman
 */
class OOFEM_EXPORT TransportGradientNeumann : public ActiveBoundaryCondition //, public Homogenization
{
public:
    TransportGradientNeumann(int n, Domain *d);

    int giveNumberOfInternalDofManagers() override { return 1; }
    DofManager *giveInternalDofManager(int i) override;

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    void postInitialize() override;

    bcType giveType() const override { return UnknownBT; }

    void scale(double s) override;

    double domainSize();

    void assembleVector(FloatArray &answer, TimeStep *tStep,
                        CharType type, ValueModeType mode,
                        const UnknownNumberingScheme &s, FloatArray *eNorm=nullptr, void*lock=nullptr) override;

    void assemble(SparseMtrx &answer, TimeStep *tStep, CharType type, const UnknownNumberingScheme &r_s,
                  const UnknownNumberingScheme &c_s, double scale=1.0, void*lock=nullptr) override;

    void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                            const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) override;

    const char *giveClassName() const override { return "TransportGradientNeumann"; }
    const char *giveInputRecordName() const override { return _IFT_TransportGradientNeumann_Name; }

    virtual void computeField(FloatArray &flux, TimeStep *tStep);
    virtual void computeTangent(FloatMatrix &tangent, TimeStep *tStep);

    void giveFluxLocationArray(IntArray &oCols, const UnknownNumberingScheme &r_s);

protected:
    /// DOF-manager containing the unknown homogenized stress.
    std :: unique_ptr< Node > mpFluxHom;
    IntArray mFluxIds;
    FloatArray mGradient;
    FloatArray mCenterCoord;
    bool dispControl;

    IntArray surfSets;

    /// Scaling factor (one array per edge with one scaling factor per GP)
    std :: vector< std :: vector< FloatArray > > eta;

    /// Help function computes phi by solving a diffusion problem on the RVE-surface.
    void computeEta();

    /// Help function that integrates the tangent contribution from a single element boundary.
    void integrateTangent(FloatMatrix &oTangent, Element *e, int iBndIndex, int surfSet, int pos);
};
} /* namespace oofem */
