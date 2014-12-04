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

#ifndef constantedgeload_h
#define constantedgeload_h

#include "boundaryload.h"

#define _IFT_ConstantEdgeLoad_Name "constantedgeload"

namespace oofem {
/**
 * This class implements a boundary load (force, moment,...) that acts
 * directly on a boundary of some finite element (on side, face, ..).
 * A boundary load is usually attribute of one or more elements.
 *
 * The boundary load describes its geometry and values (it is assumed, that user will specify
 * all necessary dofs) on element boundary using isoparametric approximation.
 * Elements can request the order of approximation (for setting up the appropriate
 * integration rule order) and the array of values (for each dof) at specific integration point
 * on the boundary.
 *
 * Elements must take care, on which boundary the load acts on (side number, ...).
 *
 * For some elements it may be better to obtain "vertex values" of boundary load to
 * compute load vector directly using exact formulae.
 *
 * This class is not restricted to structural problems. For example, in thermal
 * analysis, a boundary load load would be a  heat source.
 */
class OOFEM_EXPORT ConstantEdgeLoad : public BoundaryLoad
{
public:
    ConstantEdgeLoad(int i, Domain * d) : BoundaryLoad(i, d) { }

    virtual void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode);
    virtual int giveApproxOrder() { return 0; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual bcGeomType giveBCGeoType() const { return EdgeLoadBGT; }

    virtual const char *giveClassName() const { return "ConstantEdgeLoad"; }
    virtual const char *giveInputRecordName() const { return _IFT_ConstantEdgeLoad_Name; }

private:
    virtual void computeNArray(FloatArray &answer, const FloatArray &coords) const { answer.clear(); }
};
} // end namespace oofem
#endif // constantedgeload_h
