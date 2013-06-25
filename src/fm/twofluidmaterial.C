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

#include "twofluidmaterial.h"
#include "domain.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "engngm.h"
#include "materialinterface.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
//#include "leplic.h"

namespace oofem {

REGISTER_Material( TwoFluidMaterial );

int
TwoFluidMaterial :: checkConsistency()
{
    return this->giveMaterial(0)->checkConsistency() && this->giveMaterial(1)->checkConsistency();
}

int
TwoFluidMaterial :: hasMaterialModeCapability(MaterialMode mode)
{
    return this->giveMaterial(0)->hasMaterialModeCapability(mode) && 
           this->giveMaterial(1)->hasMaterialModeCapability(mode);
}


IRResultType
TwoFluidMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->slaveMaterial, _IFT_TwoFluidMaterial_mat);
    if ( this->slaveMaterial.giveSize() != 2 ) {
        _error("initializeFrom: mat array should have two values\n");
    }
    return IRRT_OK;
}


void
TwoFluidMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    FluidDynamicMaterial :: giveInputRecord(input);
    input.setField(this->slaveMaterial, _IFT_TwoFluidMaterial_mat);
}


double
TwoFluidMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *atTime)
{
    if ( mode == MRM_Viscosity ) {
        double vof = this->giveTempVOF(gp);
        return ( ( 1.0 - vof ) * giveMaterial(0)->giveCharacteristicValue(MRM_Viscosity, gp, atTime) +
                vof * giveMaterial(1)->giveCharacteristicValue(MRM_Viscosity, gp, atTime) );
    } else {
        _error2("giveCharacteristicValue: sorry, do not know, how to handle mode value %d", ( int ) mode);
    }

    return 0.0;
}


double
TwoFluidMaterial :: give(int aProperty, GaussPoint *gp)
//
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
//
{
    if (  aProperty == Viscosity ) {
        double vof = this->giveTempVOF(gp);
        return ( ( 1.0 - vof ) * giveMaterial(0)->give(Viscosity, gp) +
            vof * giveMaterial(1)->give(Viscosity, gp) );
    } else if ( aProperty == YieldStress ) {
        double vof = this->giveTempVOF(gp);
        return ( ( 1.0 - vof ) * giveMaterial(0)->give(YieldStress, gp) +
            vof * giveMaterial(1)->give(YieldStress, gp) );
    } else {
        _error("give: sorry, do not know, how to return any property for two fluid material");
    }

    return 0.0;
}


MaterialStatus *
TwoFluidMaterial :: CreateStatus(GaussPoint *gp) const
{
    // We can only return one of these, so we return the first one.
    this->giveMaterial(1)->giveStatus(gp);
    return this->giveMaterial(0)->giveStatus(gp);
}


void
TwoFluidMaterial :: computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep)
{
    double vof = this->giveTempVOF(gp);
    FloatArray v0, v1;
    FluidDynamicMaterialStatus *status = static_cast< FluidDynamicMaterialStatus * >( this->giveStatus(gp) );

    this->giveMaterial(0)->computeDeviatoricStressVector(v0, gp, eps, tStep);
    this->giveMaterial(1)->computeDeviatoricStressVector(v1, gp, eps, tStep);

    answer.resize(0);
    answer.add(1.0 - vof, v0);
    answer.add(vof, v1);

    status->letTempDeviatoricStressVectorBe(answer);
}

void
TwoFluidMaterial :: giveDeviatoricStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp,
                                                  TimeStep *atTime)
{
    FloatMatrix a0, a1;
    double vof = this->giveTempVOF(gp);

    this->giveMaterial(0)->giveDeviatoricStiffnessMatrix(a0, mode, gp, atTime);
    this->giveMaterial(1)->giveDeviatoricStiffnessMatrix(a1, mode, gp, atTime);

    answer.beEmptyMtrx();
    answer.add(1.0 - vof, a0);
    answer.add(vof, a1);
}

double
TwoFluidMaterial :: giveTempVOF(GaussPoint *gp)
{
    /*
     * Element* elem = gp->giveElement();
     * LEPlicElementInterface *interface = (LEPlicElementInterface*) elem->giveInterface(LEPlicElementInterfaceType);
     * if (interface) {
     * return interface->giveTempVolumeFraction();
     * } else {
     * return 0.0; // the default
     * }
     */
    FloatArray vof(2);
    MaterialInterface *mi = domain->giveEngngModel()->giveMaterialInterface( domain->giveNumber() );
    if ( mi ) {
        mi->giveElementMaterialMixture( vof, gp->giveElement()->giveNumber() );

        if ( ( vof.at(1) < 0. ) || ( vof.at(1) > 1.0 ) ) {
            _error2( "giveTempVOF: vof value out of range (vof=%lf)", vof.at(1) );
        }

        return vof.at(1);
    } else {
        return 0.0;
    }
}
} // end namespace oofem
