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

#include "fm/Materials/fluiddynamicmaterial.h"
#include "fm/Materials/nonlinearfluidmaterial.h"
#include "domain.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "engngm.h"

#include <cstdlib>

namespace oofem {
REGISTER_Material(NonlinearFluidMaterial);


void
NonlinearFluidMaterial :: initializeFrom(InputRecord &ir)
{
    FluidDynamicMaterial :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, viscosity, _IFT_NonlinearFluidMaterial_mu);
    IR_GIVE_FIELD(ir, alpha, _IFT_NonlinearFluidMaterial_alpha);
    IR_GIVE_FIELD(ir, c, _IFT_NonlinearFluidMaterial_C);
}


void
NonlinearFluidMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    FluidDynamicMaterial :: giveInputRecord(input);
    input.setField(this->viscosity, _IFT_NonlinearFluidMaterial_mu);
    input.setField(this->alpha, _IFT_NonlinearFluidMaterial_alpha);
    input.setField(this->c, _IFT_NonlinearFluidMaterial_C);
}


double
NonlinearFluidMaterial :: giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep) const
{
    return this->viscosity;
}


double
NonlinearFluidMaterial :: give(int aProperty, GaussPoint *gp) const
{
    if ( aProperty == Viscosity ) {
        return viscosity;
    } else {
        return FluidDynamicMaterial :: give(aProperty, gp);
    }
}


MaterialStatus *
NonlinearFluidMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new NonlinearFluidMaterialStatus(gp);
}


FloatArrayF<6>
NonlinearFluidMaterial :: computeDeviatoricStress3D(const FloatArrayF<6> &eps, GaussPoint *gp, TimeStep *tStep) const
{
    NonlinearFluidMaterialStatus *status = static_cast< NonlinearFluidMaterialStatus * >( this->giveStatus(gp) );

    ///@todo This doesn't look like a strain norm?
    double normeps2 = eps[0] * eps[0] + eps[1] * eps[1] + eps[2] * eps[2] + 0.5 * ( eps[3] * eps[3] + eps[4] * eps[4] +  eps[5] * eps[5] );

    auto answer = eps;
    answer[3] *= 0.5;
    answer[4] *= 0.5;
    answer[5] *= 0.5;
    answer *= 2.0 * viscosity * ( 1.0 + c * pow(normeps2, alpha * 0.5) );

    status->letTempDeviatoricStressVectorBe(answer);
    status->letTempDeviatoricStrainVectorBe(eps);
    status->letTempStrainNorm2Be(normeps2);
    return answer;
}

FloatMatrixF<6,6>
NonlinearFluidMaterial :: computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    NonlinearFluidMaterialStatus *status = static_cast< NonlinearFluidMaterialStatus * >( this->giveStatus(gp) );
    double normeps2 = status->giveTempStrainNorm2();

    auto answer = eye<6>();
    answer.at(4, 4) *= 0.5;
    answer.at(5, 5) *= 0.5;
    answer.at(6, 6) *= 0.5;

    if ( normeps2 != 0 ) {
        auto eps = status->giveTempDeviatoricStrainVector();

        eps.at(4) *= 0.5;
        eps.at(5) *= 0.5;
        eps.at(6) *= 0.5;

        auto op = dyad(eps, eps);
        answer *= 2 * viscosity * ( 1 + c * pow(normeps2, alpha * 0.5) );
        answer += (2 * viscosity * c * alpha * pow(normeps2, alpha * 0.5 - 1)) * op;
    } else {
        answer *= 2 * viscosity;
    }
    return answer;
}

int
NonlinearFluidMaterial :: checkConsistency()
{
    if ( domain->giveEngngModel()->giveEquationScalingFlag() ) {
        double scale = domain->giveEngngModel()->giveVariableScale(VST_Density);
        propertyDictionary.at('d') /= scale;

        scale = domain->giveEngngModel()->giveVariableScale(VST_Viscosity);
        this->viscosity /= scale;
    }

    return 1;
}

NonlinearFluidMaterialStatus :: NonlinearFluidMaterialStatus(GaussPoint *g) :
    FluidDynamicMaterialStatus(g)
{ }

void
NonlinearFluidMaterialStatus :: initTempStatus()
{
    FluidDynamicMaterialStatus :: initTempStatus();

    temp_deviatoricStressVector = deviatoricStressVector;
    temp_deviatoricStrainVector = deviatoricStrainRateVector;
}

void
NonlinearFluidMaterialStatus :: updateYourself(TimeStep *tStep)
{
    FluidDynamicMaterialStatus :: updateYourself(tStep);

    deviatoricStressVector = temp_deviatoricStressVector;
    deviatoricStrainRateVector = temp_deviatoricStrainVector;
}
} // end namespace oofem
