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

#ifndef line2boundaryelement_h
#define line2boundaryelement_h

#include "fmelement.h"
#include "spatiallocalizer.h"
#include "eleminterpmapperinterface.h"

namespace oofem {
class FEI2dLineQuad;

/**
 * Boundary element used for tracking solutions on arbitrary sections.
 * Also convenient for computing the RVE volume.
 * @author Mikael Öhman
 */
class Line2BoundaryElement :
    public FMElement,
    public SpatialLocalizerInterface,
    public EIPrimaryUnknownMapperInterface
{
protected:
    int boundaryNumber;
    static FEI2dLineQuad fei;

public:
    /**
     * Constructor.
     * @param n Element's number.
     * @param d Pointer to the domain to which element belongs.
     */
    Line2BoundaryElement(int n, Domain *d);
    /// Destructor.
    virtual ~Line2BoundaryElement();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) { answer.resize(0); }
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep) { answer.resize(0,0); }

    virtual void computeN(FloatArray &answer, const FloatArray &lcoords) const;

    /**
     * Gives the boundary number.
     * @todo Use regions instead?
     * @return Boundary number.
     */
    virtual int giveBoundaryNumber() const { return boundaryNumber; }
    /**
     * Computes the integral @f$ \int_S n \cdot x \mathrm{d}s @f$.
     * The normal is defined as left in the direction parameterization.
     * @todo{Move actual computations to FEI class}
     * @return Evaluated integral.
     */
    virtual double computeNXIntegral() const;

    virtual void giveDofManDofIDMask(int i, EquationID eid, IntArray &nodeDofIDMask) const;

    virtual FEInterpolation *giveInterpolation();
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_line_2; }
    virtual int computeNumberOfDofs(EquationID eid) { return 0; }

    virtual const char *giveClassName() const { return "Line2BoundaryElement"; }
    virtual classType giveClassID() const { return Line2BoundaryElementClass; }

    // Interfaces
    virtual Interface* giveInterface(InterfaceType it);

    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords) { return false; }
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &gcoords);
    virtual double SpatialLocalizerI_giveClosestPoint(FloatArray &lcoords, FloatArray &closest, const FloatArray &gcoords);

    virtual int EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode,
        TimeStep *tStep, const FloatArray &gcoords, FloatArray &answer);
    virtual void EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
        TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer);
    virtual void EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer);
};
}

#endif // line2boundaryelement_h
