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

#include "twofluidmaterial.h"
#include "domain.h"
#include "flotmtrx.h"
#include "gausspnt.h"
#include "engngm.h"
#include "materialinterface.h"
//#include "leplic.h"

namespace oofem {
int
TwoFluidMaterial :: checkConsistency()
{
    if ( this->giveMaterial(0)->giveClassID() !=
        this->giveMaterial(1)->giveClassID() ) {
        if ( ( this->giveMaterial(0)->giveClassID() == NewtonianFluidMaterialClass ) &&
            ( this->giveMaterial(1)->giveClassID() == BinghamFluidMaterialClass ) ) {
            masterMat = 1;
        } else if ( ( this->giveMaterial(0)->giveClassID() == BinghamFluidMaterialClass ) &&
                   ( this->giveMaterial(1)->giveClassID() == NewtonianFluidMaterialClass ) ) {
            masterMat = 0;
        } else {
            /*
             * This is here necessary, because these two materials can share same material status in gp.
             * It can be more general, but support for multiple statuses in gp should be implemented
             */
            _error("checkConsistency: Materials should be of the same type");
        }
    } else {
        masterMat = 0;
    }

    this->giveMaterial(0)->checkConsistency();
    this->giveMaterial(1)->checkConsistency();
    return 1;
}

int
TwoFluidMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    return this->giveMaterial(masterMat)->hasMaterialModeCapability(mode);
}


IRResultType
TwoFluidMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IntArray mats(2);
    IR_GIVE_FIELD(ir, mats, IFT_TwoFluidMaterial_mat, "mat"); // Macro
    if ( mats.giveSize() != 2 ) {
        _error("initializeFrom: mat array should have two values\n");
    }

    slaveMaterial [ 0 ] = mats.at(1);
    slaveMaterial [ 1 ] = mats.at(2);

    return IRRT_OK;
}


int
TwoFluidMaterial :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    FluidDynamicMaterial :: giveInputRecordString(str, keyword);
    sprintf(buff, " mat 2 %d %d ", slaveMaterial [ 0 ], slaveMaterial [ 1 ]);
    str += buff;

    return 1;
}


double
TwoFluidMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *atTime)
{
    if ( mode == MRM_Density ) {
        double vof = this->giveTempVOF(gp);
        return ( ( 1.0 - vof ) * giveMaterial(0)->giveCharacteristicValue(MRM_Density, gp, atTime) +
                vof * giveMaterial(1)->giveCharacteristicValue(MRM_Density, gp, atTime) );
    } else if ( mode == MRM_Viscosity ) {
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
    _error("give: sorry, do not know, how to return any property for two fluid material");
    return 0.0;
}


MaterialStatus *
TwoFluidMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates new  material status  corresponding to this class
 * NOTE:
 * we do rely on the fact, that the two salve materials are of the same type, so they can share
 * same material status !!! This is checked in checkConsistency() service.
 */
{
    return this->giveMaterial(masterMat)->CreateStatus(gp);
}


void
TwoFluidMaterial :: computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep)
{
    double vof = this->giveTempVOF(gp);
    FloatArray v0, v1;
    int i, size = eps.giveSize();
    FluidDynamicMaterialStatus *status = ( FluidDynamicMaterialStatus * ) this->giveStatus(gp);

    this->giveMaterial(0)->computeDeviatoricStressVector(v0, gp, eps, tStep);
    this->giveMaterial(1)->computeDeviatoricStressVector(v1, gp, eps, tStep);

    answer.resize(size);
    for ( i = 1; i <= size; i++ ) {
        answer.at(i) = ( 1.0 - vof ) * v0.at(i) + vof *v1.at(i);
    }

    status->letTempDeviatoricStressVectorBe(answer);
}

void
TwoFluidMaterial :: giveDeviatoricStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp,
                                                  TimeStep *atTime)
{
    FloatMatrix a0, a1;
    int i, j, size;
    double vof = this->giveTempVOF(gp);
    //FluidDynamicMaterialStatus* status = (FluidDynamicMaterialStatus*) this->giveStatus(gp);

    this->giveMaterial(0)->giveDeviatoricStiffnessMatrix(a0, mode, gp, atTime);
    this->giveMaterial(1)->giveDeviatoricStiffnessMatrix(a1, mode, gp, atTime);

    size = a1.giveNumberOfRows();
    answer.resize(size, size);
    for ( i = 1; i <= size; i++ ) {
        for ( j = 1; j <= size; j++ ) {
            answer.at(i, j) = ( 1.0 - vof ) * a0.at(i, j) + vof *a1.at(i, j);
        }
    }
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
