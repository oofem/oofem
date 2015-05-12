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

#include "newtonianfluid.h"
#include "domain.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "engngm.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(NewtonianFluidMaterial);

IRResultType
NewtonianFluidMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, viscosity, _IFT_NewtonianFluidMaterial_mu);

    return FluidDynamicMaterial :: initializeFrom(ir);
}


void
NewtonianFluidMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    FluidDynamicMaterial :: giveInputRecord(input);
    input.setField(this->viscosity, _IFT_NewtonianFluidMaterial_mu);
}


double
NewtonianFluidMaterial :: giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep)
{
    return this->viscosity;
}


double
NewtonianFluidMaterial :: give(int aProperty, GaussPoint *gp)
{
    if ( aProperty == Viscosity ) {
        return viscosity;
    } else if ( aProperty == YieldStress ) {
        return 0.0;
    } else {
        return FluidDynamicMaterial :: give(aProperty, gp);
    }
}


MaterialStatus *
NewtonianFluidMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new FluidDynamicMaterialStatus(1, this->giveDomain(), gp);
}


void
NewtonianFluidMaterial :: computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep)
{
    int size = eps.giveSize();
    answer.resize(size);

    if ( gp->giveMaterialMode() == _2dFlow ) {
        double ekk = eps.at(1) + eps.at(2);

        answer = {
            2.0 * viscosity * ( eps.at(1) - ekk / 3.0 ),
            2.0 * viscosity * ( eps.at(2) - ekk / 3.0 ),
            eps.at(3) * viscosity
        };
    } else if ( gp->giveMaterialMode() == _2dAxiFlow ) {
#if 1
        double ekk = eps.at(1) + eps.at(2) + eps.at(3);

        answer = {
            2.0 * viscosity * ( eps.at(1) - ekk / 3.0 ),
            2.0 * viscosity * ( eps.at(2) - ekk / 3.0 ),
            2.0 * viscosity * ( eps.at(3) - ekk / 3.0 ),
            eps.at(4) * viscosity
        };
#else
        answer = {
            2.0 * viscosity * ( eps.at(1) ),
            2.0 * viscosity * ( eps.at(2) ),
            2.0 * viscosity * ( eps.at(3) ),
            eps.at(4) * viscosity
        };
#endif
    } else if ( gp->giveMaterialMode() == _3dFlow ) {
        double ekk = eps.at(1) + eps.at(2) + eps.at(3);

        answer = {
            2.0 * viscosity * ( eps.at(1) - ekk / 3.0 ),
            2.0 * viscosity * ( eps.at(2) - ekk / 3.0 ),
            2.0 * viscosity * ( eps.at(3) - ekk / 3.0 ),
            eps.at(4) * viscosity,
            eps.at(5) * viscosity,
            eps.at(6) * viscosity
        };
    }  else {
        OOFEM_ERROR("unsupported material mode");
    }


    static_cast< FluidDynamicMaterialStatus * >( this->giveStatus(gp) )->letDeviatoricStressVectorBe(answer);
    static_cast< FluidDynamicMaterialStatus * >( this->giveStatus(gp) )->letDeviatoricStrainRateVectorBe(eps);
}

void
NewtonianFluidMaterial :: giveDeviatoricStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp,
                                                        TimeStep *tStep)
{
    if ( ( gp->giveMaterialMode() == _2dFlow ) ) {
        answer.resize(3, 3);
        answer.zero();

        answer.at(1, 1) = answer.at(2, 2) = 2.0 * viscosity * ( 2. / 3. );
        answer.at(1, 2) = answer.at(2, 1) = -2.0 * viscosity * ( 1. / 3. );
        answer.at(3, 3) = viscosity;
    } else if ( gp->giveMaterialMode() == _2dAxiFlow ) {
#if 1
        answer.resize(4, 4);
        answer.zero();

        answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = 2.0 * viscosity * ( 2. / 3. );
        answer.at(1, 2) = answer.at(1, 3) = -2.0 * viscosity * ( 1. / 3. );
        answer.at(2, 1) = answer.at(2, 3) = -2.0 * viscosity * ( 1. / 3. );
        answer.at(3, 1) = answer.at(3, 2) = -2.0 * viscosity * ( 1. / 3. );

        answer.at(4, 4) = viscosity;
#else
        answer.resize(4, 4);
        answer.zero();

        answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = 2.0 * viscosity;
        answer.at(4, 4) = viscosity;
#endif
    } else if ( gp->giveMaterialMode() == _3dFlow ) {
        answer.resize(6, 6);
        answer.zero();

        answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = 2.0 * viscosity * ( 2. / 3. );
        answer.at(1, 2) = answer.at(1, 3) = -2.0 * viscosity * ( 1. / 3. );
        answer.at(2, 1) = answer.at(2, 3) = -2.0 * viscosity * ( 1. / 3. );
        answer.at(3, 1) = answer.at(3, 2) = -2.0 * viscosity * ( 1. / 3. );

        answer.at(4, 4) = answer.at(5, 5) = answer.at(6, 6) = viscosity;
    } else {
        OOFEM_ERROR("unsupportted material mode");
    }
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
