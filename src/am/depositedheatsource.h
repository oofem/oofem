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

#ifndef depositedheatsource_h
#define depositedheatsource_h

#include "bodyload.h"
#include "bcgeomtype.h"
#include "valuemodetype.h"

#define _IFT_DepositedHeatSource_Name "depositedheatsource"
#define _IFT_DepositedHeatSource_depositedmaterialid "mat"
#define _IFT_DepositedHeatSource_depositiontemperature "tdep"
#define _IFT_DepositedHeatSource_depositedmassfractionfunction "depmassfractionfunction" // mass fraction of deposited material
#define _IFT_DepositedHeatSource_power "value"

namespace oofem {
/**
 * This class implements a volumetric heat source to account for incrementally deposited material (with specific given temperature) in element
 *
 */
class OOFEM_EXPORT DepositedHeatSource : public BodyLoad
{
protected:
    int depositedMaterialID; /// ID of the deposited material; not used at the moment
    double powervalue; /// power = specificHeat * density
    double depositionTemperature; /// Temperature of the deposited material
    int depositedMassFractionFunction; /// Mass fraction of the deposited material
    
public:
    /// Constructor
    DepositedHeatSource(int i, Domain * d) : BodyLoad(i, d) { }
    DepositedHeatSource(int i, Domain * d, int depositedMaterialID, double powervalue, double depositionTemperature, int depositedMassFractionFunction) :
        BodyLoad(i, d), depositedMaterialID(depositedMaterialID), powervalue(powervalue), depositionTemperature(depositionTemperature), depositedMassFractionFunction(depositedMassFractionFunction) { }

    void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode) override;
    void computeValueAt(FloatArray &answer, TimeStep *tStep, GaussPoint *gp, ValueModeType mode) override;

    bcValType giveBCValType() const override { return ForceLoadBVT; }
    bcGeomType giveBCGeoType() const override { return BodyLoadBGT; }

    void setComponents(int depositedMaterialID, double powervalue, double depositionTemperature, double depositedMassFractionFunction)
    {
        this->depositedMaterialID = depositedMaterialID;
        this->depositionTemperature = depositionTemperature;
        this->depositedMassFractionFunction = depositedMassFractionFunction;
        this->powervalue = powervalue;
    }   
    FormulationType giveFormulationType() override  { return FT_Global; }
    void initializeFrom(InputRecord &ir) override;
    
    const char *giveClassName() const override { return "DepositedHeatSource"; }
    const char *giveInputRecordName() const override { return _IFT_DepositedHeatSource_Name; }
};
} // end namespace oofem
#endif // depositedheatsource_h
