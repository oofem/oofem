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

#ifndef interactionload_h
#define interactionload_h

#include "linearedgeload.h"
#include "gausspoint.h"
#include "floatarray.h"

///@name Input fields for LinearEdgeLoad
//@{
#define _IFT_InteractionLoad_Name "interactionload"
#define _IFT_InteractionLoad_CoupledParticles "coupledparticles"
//#define _IFT_LinearEdgeLoad_formulation "formulation"
//#define _IFT_LinearEdgeLoad_startcoord "sc"
//#define _IFT_LinearEdgeLoad_endcoord "ec"
//@}

namespace oofem {
class TimeStep;

/**
 * This class implements a fluid-to-solid interface in the FluidStructureProblem.
 * Its attributes are the corresponding fluid PFEMParticles providing the pressure.
 * The pressure load acts in normal edge direction, so a transformation to the global
 * coordinate system is performed.
 */
class OOFEM_EXPORT InteractionLoad : public LinearEdgeLoad
{
protected:
    /// Coordinates of start and end point
    //FloatArray startCoords, endCoords;
    //FormulationType formulation;
    IntArray coupledParticles;

public:
    InteractionLoad(int i, Domain *d) : LinearEdgeLoad(i, d), coupledParticles(2) { }

    virtual void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode);

    virtual int giveApproxOrder() { return 1; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual bcGeomType giveBCGeoType() const { return EdgeLoadBGT; }
    virtual FormulationType giveFormulationType() { return formulation; }

    virtual const char *giveClassName() const { return "InteractionLoad"; }
    virtual const char *giveInputRecordName() const { return _IFT_InteractionLoad_Name; }

protected:
    virtual void computeNArray(FloatArray &answer, const FloatArray &coords) const;
};
} // end namespace oofem
#endif // interactionload_h
