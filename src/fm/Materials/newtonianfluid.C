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

#include "fm/Materials/newtonianfluid.h"
#include "domain.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "engngm.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(NewtonianFluidMaterial);

void
NewtonianFluidMaterial :: initializeFrom(InputRecord &ir)
{
    FluidDynamicMaterial :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, viscosity, _IFT_NewtonianFluidMaterial_mu);

}


void
NewtonianFluidMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    FluidDynamicMaterial :: giveInputRecord(input);
    input.setField(this->viscosity, _IFT_NewtonianFluidMaterial_mu);
}


double
NewtonianFluidMaterial :: giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep) const
{
    return this->viscosity;
}


double
NewtonianFluidMaterial :: give(int aProperty, GaussPoint *gp) const
{
    if ( aProperty == Viscosity ) {
        return viscosity;
    } else if ( aProperty == YieldStress ) {
        return 0.0;
    } else {
        return FluidDynamicMaterial :: give(aProperty, gp);
    }
}


std::unique_ptr<MaterialStatus> 
NewtonianFluidMaterial :: CreateStatus(GaussPoint *gp) const
{
    return std::make_unique<FluidDynamicMaterialStatus>(gp);
}


FloatArrayF<6>
NewtonianFluidMaterial :: computeDeviatoricStress3D(const FloatArrayF<6> &eps, GaussPoint *gp, TimeStep *tStep) const
{
    double ekk = eps.at(1) + eps.at(2) + eps.at(3);

    FloatArrayF<6> stress = {
        2.0 * viscosity * ( eps.at(1) - ekk / 3.0 ),
        2.0 * viscosity * ( eps.at(2) - ekk / 3.0 ),
        2.0 * viscosity * ( eps.at(3) - ekk / 3.0 ),
        eps.at(4) * viscosity,
        eps.at(5) * viscosity,
        eps.at(6) * viscosity
    };

    static_cast< FluidDynamicMaterialStatus * >( this->giveStatus(gp) )->letDeviatoricStressVectorBe(stress);
    static_cast< FluidDynamicMaterialStatus * >( this->giveStatus(gp) )->letDeviatoricStrainRateVectorBe(eps);

    return stress;
}

FloatMatrixF<6,6>
NewtonianFluidMaterial :: computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    return 2 * viscosity * I_dev6;
}


int
NewtonianFluidMaterial :: checkConsistency()
{
    ///@todo Fix this, checkConsistency shouldn't be a replacement for "post-initialization" but should be completely optional.
    if ( domain->giveEngngModel()->giveEquationScalingFlag() ) {
        double scale;
        scale = domain->giveEngngModel()->giveVariableScale(VST_Density);
        propertyDictionary.at('d') /= scale;

        scale = domain->giveEngngModel()->giveVariableScale(VST_Viscosity);
        this->viscosity /= scale;
    }

    return 1;
}
} // end namespace oofem
