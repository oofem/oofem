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

#ifndef pointload_h
#define pointload_h

#include "bodyload.h"

///@name Input fields for PointLoad
//@{
#define _IFT_PointLoad_Name "pointload"
#define _IFT_PointLoad_ndofs "ndofs"
#define _IFT_PointLoad_coords "coords"
#define _IFT_PointLoad_loadtype "loadtype"
#define _IFT_PointLoad_cstype "cstype"
//@}

namespace oofem {
class TimeStep;

/**
 * Abstract base class representing a point load (force, momentum, ...) that acts
 * directly on or inside of some finite element.
 *
 * Methods for returning values of load components and returning the position of
 * load are provided.
 *
 * @note
 * This class is not restricted to structural problems. For example, in thermal
 * analysis, a point load load could be a point heat source.
 */
class OOFEM_EXPORT PointLoad : public BodyLoad
{
protected:
    /// Load type (its physical meaning).
    bcType lType;
    /// Load coordinate system.
    CoordSystType coordSystemType;
    /// Additional properties (coordinates, point of application).
    FloatArray coords;

public:
    /**
     * Constructor. Creates Boundary Load object with given number, belonging to given domain.
     * @param n Load number.
     * @param d Domain to which new object will belongs.
     */
    PointLoad(int n, Domain * d) : BodyLoad(n, d), coordSystemType(CST_Global) { }

    void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode) override;
    /**
     * Gives coordinates of the receiver
     */
    const FloatArray & giveCoordinates() const { return coords; }

    CoordSystType giveCoordSystMode() override { return coordSystemType; }
    //FormulationType giveFormulationType () override { return FT_Entity; }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    bcType giveType() const override { return lType; }
    bcGeomType giveBCGeoType() const override { return PointLoadBGT; }
    const char *giveClassName() const override { return "PointLoad"; }
    const char *giveInputRecordName() const override { return _IFT_PointLoad_Name; }
};
} // end namespace oofem
#endif // pointload_h
