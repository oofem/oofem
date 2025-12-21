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

#ifndef usrdefboundaryload_h
#define usrdefboundaryload_h

#include "boundaryload.h"

///@name Input fields for ConstantSurfaceLoad
//@{
#define _IFT_UsrDefBoundaryLoad_intensityfunction "intensityfunction"
#define _IFT_UsrDefBoundaryLoad_Name "usrdefboundaryload"
#define _IFT_UsrDefBoundaryLoad_GeomType "geomtype"
#define _IFT_UsrDefBoundaryLoad_approxorder "approxorder"
//@}


namespace oofem {
/**
 * This class implements a boundary load (force, moment,...) that acts
 * directly on a boundary of some finite element (on side, face, ..).
 * A boundary load is usually attribute of one or more elements.
 *
 * The boundary load describes its intensity as a function of spatial coordinates.
 * Elements can request the order of approximation (for setting up the appropriate
 * integration rule order).
 *
 * Elements take care, on which boundary the load acts on (side number, ...).
 *
 * This class is not restricted to structural problems. For example, in thermal
 * analysis, a boundary load load would be a  heat source.
 */
class OOFEM_EXPORT UsrDefBoundaryLoad : public BoundaryLoad
{
    protected:
    /// Reference to user defined function that computes the load intensity
    int intensityFunction;
    /// bcgeotype of the boundary load
    bcGeomType myGeomType = SurfaceLoadBGT; // default value
    /// approximation order to set up the appropriate integration rule 
    int approxOrder = 0;
public:
    UsrDefBoundaryLoad(int i, Domain * d);

    // Overloaded methods:
    void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode) override;
    virtual void computeComponentArrayAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode) override;
    int giveApproxOrder() override { return approxOrder; }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    bcGeomType giveBCGeoType() const override { return myGeomType; }
    FormulationType giveFormulationType() override { return FT_Global; }


    const char *giveClassName() const override { return "UsrDefBoundaryLoad"; }
    const char *giveInputRecordName() const override { return _IFT_UsrDefBoundaryLoad_Name; }

private:
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
    void computeNArray(FloatArray &answer, const FloatArray &coords) const override { answer.clear(); }

};
} // end namespace oofem
#endif // usrdefboundaryload_h
