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

#include "depositedheatsource.h"
#include "timestep.h"
#include "function.h"
#include "classfactory.h"
#include "material.h"
#include "engngm.h"
#include "field.h"

namespace oofem {
REGISTER_BoundaryCondition(DepositedHeatSource);

void DepositedHeatSource :: computeValueAt (FloatArray &answer, TimeStep *tStep, GaussPoint *gp, ValueModeType mode)
{
    answer.resize(1);

    if (mode == VM_Total) {

        double depositedMassFractionIncrement = this->domain->giveFunction(this->depositedMassFractionFunction)->evaluate(tStep, VM_Incremental);
        if (depositedMassFractionIncrement > 0.0) {

        Element* elem = gp->giveElement();
        FloatArray nodaltemp;
        elem->computeVectorOf({T_f}, VM_TotalIntrinsic, tStep, nodaltemp);
        double temperature = nodaltemp.sum() / nodaltemp.size(); // average element value
        double factor = 1.0; // this->giveTimeFunction()->evaluate(tStep, mode);
            answer.at(1) = depositedMassFractionIncrement*this->powervalue*(this->depositionTemperature-temperature)/tStep->giveTimeIncrement(); // W/m3
            answer.times(factor);
        } else {
            answer.zero();
        }
    } else {
        OOFEM_ERROR("DepositedHeatSource::computeValueAt - unsupported ValueModeType %d\n", mode);
    }
    
}

void DepositedHeatSource :: computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode)
{
    answer.resize(1);

    if (mode == VM_Total) {

        double depositedMassFractionIncrement = this->domain->giveFunction(this->depositedMassFractionFunction)->evaluate(tStep, VM_Incremental);
        if (depositedMassFractionIncrement > 0.0) {
            FloatArray val(1);
            //Material *mat = this->domain->giveMaterial(depositedMaterialID); // just to check material existence
            //double rho = 3.0; //mat ? mat->give('d', nullptr) : 0.0;
            //double specificHeat = 3.0; // mat ? mat->give('c', nullptr) : 0.0;
            this->domain->giveEngngModel()->giveField(FT_Temperature, tStep)->evaluateAt(val, coords, mode, tStep);
            double temperature = val.at(1);
            // get massFraction Increment over the time step
            double factor = 1.0; // this->giveTimeFunction()->evaluate(tStep, mode);
            answer.at(1) = depositedMassFractionIncrement*this->powervalue*(this->depositionTemperature-temperature)/tStep->giveTimeIncrement(); // W/m3
            answer.times(factor);
        } else {
            answer.zero();
        }
    } else {
        OOFEM_ERROR("DepositedHeatSource::computeValueAt - unsupported ValueModeType %d\n", mode);
    }
}

void DepositedHeatSource::initializeFrom(InputRecord &ir)
{
    GeneralBoundaryCondition::initializeFrom(ir);
    IR_GIVE_FIELD (ir, depositedMaterialID, _IFT_DepositedHeatSource_depositedmaterialid);
    IR_GIVE_FIELD (ir, depositionTemperature, _IFT_DepositedHeatSource_depositiontemperature);
    IR_GIVE_FIELD (ir, depositedMassFractionFunction, _IFT_DepositedHeatSource_depositedmassfractionfunction);
    IR_GIVE_FIELD (ir, powervalue, _IFT_DepositedHeatSource_power);
}

} // end namespace oofem
