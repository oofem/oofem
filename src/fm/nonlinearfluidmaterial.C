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

#include "fluiddynamicmaterial.h"
#include "nonlinearfluidmaterial.h"
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


IRResultType
NonlinearFluidMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, viscosity, _IFT_NonlinearFluidMaterial_mu);
    IR_GIVE_FIELD(ir, alpha, _IFT_NonlinearFluidMaterial_alpha);
    IR_GIVE_FIELD(ir, c, _IFT_NonlinearFluidMaterial_C);

    return FluidDynamicMaterial :: initializeFrom(ir);
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
NonlinearFluidMaterial :: giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep)
{
    return this->viscosity;
}


double
NonlinearFluidMaterial :: give(int aProperty, GaussPoint *gp)
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
    return new NonlinearFluidMaterialStatus(1, this->giveDomain(), gp);
}


void
NonlinearFluidMaterial :: computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep)
{
    NonlinearFluidMaterialStatus *status = static_cast< NonlinearFluidMaterialStatus * >( this->giveStatus(gp) );

    double normeps2;

    answer = eps;
    if ( eps.giveSize() == 3 ) {
        normeps2 = eps.at(1) * eps.at(1) + eps.at(2) * eps.at(2) + 0.5 * ( eps.at(3) * eps.at(3) );
        answer.at(3) *= 0.5;
    } else if ( eps.giveSize() == 4 ) {
        normeps2 = eps.at(1) * eps.at(1) + eps.at(2) * eps.at(2) + eps.at(3) * eps.at(3) + 0.5 * ( eps.at(4) * eps.at(4) );
        answer.at(4) *= 0.5;
    } else {
        normeps2 = eps.at(1) * eps.at(1) + eps.at(2) * eps.at(2) + eps.at(3) * eps.at(3) + 0.5 * ( eps.at(4) * eps.at(4) + eps.at(5) * eps.at(5) +  eps.at(6) * eps.at(6) );
        answer.at(4) *= 0.5;
        answer.at(5) *= 0.5;
        answer.at(6) *= 0.5;
    }

    answer.times( 2.0 * viscosity * ( 1.0 + c * pow(normeps2, alpha * 0.5) ) );

    status->letTempDeviatoricStressVectorBe(answer);
    status->letTempDeviatoricStrainVectorBe(eps);
    status->letTempStrainNorm2Be(normeps2);
}

void
NonlinearFluidMaterial :: giveDeviatoricStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp,
                                                        TimeStep *tStep)
{
    FloatArray eps;
    double normeps2;

    NonlinearFluidMaterialStatus *status = static_cast< NonlinearFluidMaterialStatus * >( this->giveStatus(gp) );
    eps = status->giveTempDeviatoricStrainVector();
    normeps2 = status->giveTempStrainNorm2();

    answer.resize( eps.giveSize(), eps.giveSize() );
    answer.zero();
    for ( int i = 1; i <= answer.giveNumberOfRows(); i++ ) {
        answer.at(i, i) = 1.;
    }
    if ( eps.giveSize() == 3 ) {
        answer.at(3, 3) *= 0.5;
    } else if ( eps.giveSize() == 4 ) {
        answer.at(4, 4) *= 0.5;
    } else {
        answer.at(4, 4) *= 0.5;
        answer.at(5, 5) *= 0.5;
        answer.at(6, 6) *= 0.5;
    }

    if ( normeps2 != 0 ) {
        FloatMatrix op;
        if ( eps.giveSize() == 3 ) {
            eps.at(3) *= 0.5;
        } else if ( eps.giveSize() == 4 ) {
            eps.at(4) *= 0.5;
        } else {
            eps.at(4) *= 0.5;
            eps.at(5) *= 0.5;
            eps.at(6) *= 0.5;
        }
        op.beDyadicProductOf(eps, eps);
        answer.times( 2 * viscosity * ( 1 + c * pow(normeps2, alpha * 0.5) ) );
        answer.add(2 * viscosity * c * alpha * pow(normeps2, alpha * 0.5 - 1), op);
    } else {
        answer.times(2 * viscosity);
    }
}

int
NonlinearFluidMaterial :: checkConsistency()
{
    if ( domain->giveEngngModel()->giveEquationScalingFlag() ) {
        double scale;
        scale = domain->giveEngngModel()->giveVariableScale(VST_Density);
        propertyDictionary.at('d') /= scale;

        scale = domain->giveEngngModel()->giveVariableScale(VST_Viscosity);
        this->viscosity /= scale;
    }

    return 1;
}

NonlinearFluidMaterialStatus :: NonlinearFluidMaterialStatus(int n, Domain *d, GaussPoint *g) :
    FluidDynamicMaterialStatus(n, d, g),
    temp_deviatoricStressVector(),
    temp_deviatoricStrainVector(),
    temp_norm2(0)
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
