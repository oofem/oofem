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

#define _IFT_Line2SurfaceTension_Name "line2surfacetension"

namespace oofem {
class FEI2dLineQuad;

/**
 * 3 node line elements for surface tension.
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

    virtual FEInterpolation *giveInterpolation() const;

    virtual void computeTangent(FloatMatrix &answer, TimeStep *tStep);

    virtual void computeInternalForcesVector(FloatArray &answer, ValueModeType mode, TimeStep *tStep);

    virtual int SpatialLocalizerI_containsPoint(const FloatArray &gcoords) { return false; }
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &gcoords);

    virtual int EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode,
                                                                 TimeStep *tStep, const FloatArray &gcoords,
                                                                 FloatArray &answer);

    virtual void EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
                                                                 TimeStep *tStep, const FloatArray &lcoords,
                                                                 FloatArray &answer);

    virtual int computeNumberOfDofs() { return 6; }

    virtual const char *giveClassName() const { return "Line2SurfaceTension"; }
    virtual const char *giveInputRecordName() const { return _IFT_Line2SurfaceTension_Name; }
};
}

#endif // line2surfacetension_h
