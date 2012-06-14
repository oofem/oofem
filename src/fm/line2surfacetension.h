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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef line2surfacetension_h
#define line2surfacetension_h

#include "linesurfacetension.h"

namespace oofem {
class FEI2dLineQuad;

/**
 * 3 node line elements for surface tension.
 * @see LineSurfaceTension
 * @author Mikael Öhman
 */
class Line2SurfaceTension :
    public LineSurfaceTension
{
protected:
    static FEI2dLineQuad fei;

public:
    /**
     * Constructor.
     * @param n Element's number.
     * @param d Pointer to the domain to which element belongs.
     */
    Line2SurfaceTension(int n, Domain *d);
    /// Destructor.
    virtual ~Line2SurfaceTension();

    void computeTangent(FloatMatrix &answer, TimeStep *tStep);
    void computeLoadVector(FloatArray &answer, ValueModeType mode, TimeStep *tStep);

    /**
     * Computes the integral @f$ \int n \cdot x \mathrm{d}s @f$.
     * The normal is defined as left in the direction parameterization.
     * @return Evaluated integral.
     */
    virtual double computeNXIntegral() const;

    virtual void computeN(FloatArray &answer, const FloatArray &lcoords) const;

    virtual double SpatialLocalizerI_giveClosestPoint(FloatArray &lcoords, FloatArray &closest, const FloatArray &gcoords);
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &gcoords) { return false; }

    virtual int EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode,
                                                                 TimeStep *tStep, const FloatArray &gcoords,
                                                                 FloatArray &answer);

    virtual void EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
                                                                 TimeStep *tStep, const FloatArray &lcoords,
                                                                 FloatArray &answer);

    virtual Element_Geometry_Type giveGeometryType() const { return EGT_line_2; }
    virtual int computeNumberOfDofs(EquationID ut) { return 6; }

    virtual const char *giveClassName() const { return "Line2SurfaceTension"; }
    virtual classType giveClassID() const { return Line2SurfaceTensionElementClass; }
};
}

#endif // line2surfacetension_h
