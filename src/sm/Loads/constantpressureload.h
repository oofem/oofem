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

#ifndef constantpressureload_h
#define constantpressureload_h

#include "boundaryload.h"

///@name Input fields for ConstantPressureLoad
//@{
#define _IFT_ConstantPressureLoad_LoadOffset "loadoffset"
#define _IFT_ConstantPressureLoad_Name "constantpressureload"
//@}

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
class ConstantPressureLoad : public BoundaryLoad
{
public:
    //ConstantPressureLoad(int i, Domain *d) : BoundaryLoad(i, d)
    ConstantPressureLoad(int i, Domain * d);

    // Overloaded methods:
    virtual void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode);
    virtual void computeValues(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, const IntArray &dofids, ValueModeType mode);
    virtual int giveApproxOrder() { return 0; }

    /**
     * Sets a new load vector.
     * @param newValue New load.
     */
    void updateLoad(const FloatArray &newValue) { componentArray = newValue; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual bcGeomType giveBCGeoType() const { return SurfaceLoadBGT; }

    virtual const char *giveClassName() const { return "ConstantPressureLoad"; }
    virtual const char *giveInputRecordName() const { return _IFT_ConstantPressureLoad_Name; }
    double giveLoadOffset() { return this->loadOffset; }
private:
    virtual void computeNArray(FloatArray &answer, const FloatArray &coords) const { answer.clear(); }
    double loadOffset;  // xi-coord offset of load. xi=-1 -> bottom, xi=0 -> midsurface (default), xi=1 -> top surface
};
} // end namespace oofem
#endif
