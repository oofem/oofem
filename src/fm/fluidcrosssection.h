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

#ifndef fluidcrosssection_h
#define fluidcrosssection_h

#include "crosssection.h"

#define _IFT_FluidCrossSection_Name "fluidcs"
#define _IFT_FluidCrossSection_material "mat"

namespace oofem {
class FluidDynamicMaterial;

/**
 * Fluid cross-section. It's functionality is essentially only to keep track of the material used.
 *
 * @author Mikael Öhman
 */
class OOFEM_EXPORT FluidCrossSection : public CrossSection
{
protected:
    int matNumber;

public:
    /**
     * Constructor. Creates cross section with number n belonging to domain d.
     * @param n Cross section number.
     * @param d Domain for cross section.
     */
    FluidCrossSection(int n, Domain * d);
    /// Destructor.
    virtual ~FluidCrossSection();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual int checkConsistency();

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep);

    virtual double giveDensity(GaussPoint *gp);
    FluidDynamicMaterial *giveFluidMaterial();

    virtual int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp);
    virtual int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp);
    virtual int estimatePackSize(DataStream &buff, GaussPoint *gp);

    virtual const char *giveClassName() const { return "FluidCrossSection"; }
    virtual const char *giveInputRecordName() const { return _IFT_FluidCrossSection_Name; }
};
} // end namespace oofem
#endif // fluidcrosssection_h
