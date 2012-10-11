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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "newtonianfluid.h"
#include "fluiddynamicmaterial.h"
#include "domain.h"
#include "flotmtrx.h"
#include "gausspnt.h"
#include "engngm.h"

namespace oofem {
int
NewtonianFluidMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    if ( ( mode == _2dFlow ) || ( mode == _3dFlow ) ) {
        return 1;
    }

    return 0;
}


IRResultType
NewtonianFluidMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->FluidDynamicMaterial :: initializeFrom(ir);
    // we use rather object's member data than to store data into slow
    // key-val dictionary with lot of memory allocations

    IR_GIVE_FIELD(ir, viscosity, IFT_NewtonianFluidMaterial_mu, "mu"); // Macro

    return IRRT_OK;
}


int
NewtonianFluidMaterial :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    FluidDynamicMaterial :: giveInputRecordString(str, keyword);
    sprintf(buff, " mu %e ", this->viscosity);
    str += buff;

    return 1;
}


double
NewtonianFluidMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                                  GaussPoint *gp,
                                                  TimeStep *atTime)
{
    if ( mode == MRM_Density ) {
        return this->give('d', gp);
    } else if ( mode == MRM_Viscosity ) {
        return this->viscosity;
    } else {
        return FluidDynamicMaterial :: giveCharacteristicValue(mode, gp, atTime);
    }
}


double
NewtonianFluidMaterial :: give(int aProperty, GaussPoint *gp)
//
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
//
{
       if ( ( aProperty == Viscosity ) ) {
        return viscosity;
    } else if ( ( aProperty == YieldStress ) ) {
        return 0.0;
    } else {
        return FluidDynamicMaterial :: give(aProperty, gp);
    }
}


MaterialStatus *
NewtonianFluidMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates new  material status  corresponding to this class
 */
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

        //if (fabs(ekk) > 1.e-3) printf ("Non-zero volumetric strain (ev=%e), element %d)\n", ekk, gp->giveElement()->giveNumber());

        answer.at(1) = 2.0 * viscosity * ( eps.at(1) - ekk / 3.0 );
        answer.at(2) = 2.0 * viscosity * ( eps.at(2) - ekk / 3.0 );
        answer.at(3) = eps.at(3) * viscosity;
    } else if ( gp->giveMaterialMode() == _2dAxiFlow ) {
#if 1
        double ekk = eps.at(1) + eps.at(2) + eps.at(3);

        answer.at(1) = 2.0 * viscosity * ( eps.at(1) - ekk / 3.0 );
        answer.at(2) = 2.0 * viscosity * ( eps.at(2) - ekk / 3.0 );
        answer.at(3) = 2.0 * viscosity * ( eps.at(3) - ekk / 3.0 );
        answer.at(4) = eps.at(4) * viscosity;
#else
        answer.at(1) = 2.0 * viscosity * ( eps.at(1) );
        answer.at(2) = 2.0 * viscosity * ( eps.at(2) );
        answer.at(3) = 2.0 * viscosity * ( eps.at(3) );
        answer.at(4) = eps.at(4) * viscosity;
#endif
    } else if ( gp->giveMaterialMode() == _3dFlow ) {
        double ekk = eps.at(1) + eps.at(2) + eps.at(3);

        answer.at(1) = 2.0 * viscosity * ( eps.at(1) - ekk / 3.0 );
        answer.at(2) = 2.0 * viscosity * ( eps.at(2) - ekk / 3.0 );
        answer.at(3) = 2.0 * viscosity * ( eps.at(3) - ekk / 3.0 );
        answer.at(4) = eps.at(4) * viscosity;
        answer.at(5) = eps.at(5) * viscosity;
        answer.at(6) = eps.at(6) * viscosity;
    }  else {
        _error("computeDeviatoricStressVector: unsuported material mode");
    }

    ( ( FluidDynamicMaterialStatus * ) this->giveStatus(gp) )->letTempDeviatoricStressVectorBe(answer);
}

void
NewtonianFluidMaterial :: giveDeviatoricStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp,
                                                        TimeStep *atTime)
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
        _error("giveDeviatoricStiffnessMatrix: unsupportted material mode");
    }
}

int
NewtonianFluidMaterial :: checkConsistency()
{
    if ( domain->giveEngngModel()->giveEquationScalingFlag() ) {
        double scale;
        scale = domain->giveEngngModel()->giveVariableScale(VST_Density);
        propertyDictionary->at('d') /= scale;

        scale = domain->giveEngngModel()->giveVariableScale(VST_Viscosity);
        this->viscosity /= scale;
    }

    return 1;
}
} // end namespace oofem
