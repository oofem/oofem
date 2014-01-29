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
REGISTER_Material(TwoFluidMaterial);

int
TwoFluidMaterial :: checkConsistency()
{
    return this->giveMaterial(0)->checkConsistency() &&
           this->giveMaterial(1)->checkConsistency();
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
TwoFluidMaterial :: giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep)
{
    TwoFluidMaterialStatus *status = static_cast< TwoFluidMaterialStatus* >( this->giveStatus(gp) );
    double vof = this->giveTempVOF(gp);
    return ( 1.0 - vof ) * giveMaterial(0)->giveEffectiveViscosity(status->giveSlaveGaussPoint0(), tStep) +
           vof *giveMaterial(1)->giveEffectiveViscosity(status->giveSlaveGaussPoint1(), tStep);
}


double
TwoFluidMaterial :: give(int aProperty, GaussPoint *gp)
{
    if ( aProperty == Viscosity || aProperty == 'd' ) {
        TwoFluidMaterialStatus *status = static_cast< TwoFluidMaterialStatus* >( this->giveStatus(gp) );
        double vof = this->giveTempVOF(gp);
        return ( 1.0 - vof ) * giveMaterial(0)->give(aProperty, status->giveSlaveGaussPoint0()) +
                         vof * giveMaterial(1)->give(aProperty, status->giveSlaveGaussPoint1());
#if 0
    } else if ( aProperty == YieldStress ) {
        double vof = this->giveTempVOF(gp);
        return ( 1.0 - vof ) * giveMaterial(0)->give(YieldStress, status->giveSlaveGaussPoint0()) +
                         vof * giveMaterial(1)->give(YieldStress, status->giveSlaveGaussPoint1());
#endif
    } else {
        _error("give: sorry, do not know, how to return any property for two fluid material");
    }

    return 0.0;
}


MaterialStatus *
TwoFluidMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new TwoFluidMaterialStatus(1, this->giveDomain(), gp, this->slaveMaterial);
}


void
TwoFluidMaterial :: computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep)
{
    double vof = this->giveTempVOF(gp);
    FloatArray v0, v1;
    TwoFluidMaterialStatus *status = static_cast< TwoFluidMaterialStatus* >( this->giveStatus(gp) );

    this->giveMaterial(0)->computeDeviatoricStressVector(v0, status->giveSlaveGaussPoint0(), eps, tStep);
    this->giveMaterial(1)->computeDeviatoricStressVector(v1, status->giveSlaveGaussPoint1(), eps, tStep);

    answer.resize(0);
    answer.add(1.0 - vof, v0);
    answer.add(vof, v1);

    status->letTempDeviatoricStrainRateVectorBe(eps);
    status->letTempDeviatoricStressVectorBe(answer);
}

void
TwoFluidMaterial :: giveDeviatoricStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp,
                                                  TimeStep *tStep)
{
    FloatMatrix a0, a1;
    double vof = this->giveTempVOF(gp);
    TwoFluidMaterialStatus *status = static_cast< TwoFluidMaterialStatus* >( this->giveStatus(gp) );

    this->giveMaterial(0)->giveDeviatoricStiffnessMatrix(a0, mode, status->giveSlaveGaussPoint0(), tStep);
    this->giveMaterial(1)->giveDeviatoricStiffnessMatrix(a1, mode, status->giveSlaveGaussPoint1(), tStep);

    answer.beEmptyMtrx();
    answer.add(1.0 - vof, a0);
    answer.add(vof, a1);
}

double
TwoFluidMaterial :: giveTempVOF(GaussPoint *gp)
{
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




TwoFluidMaterialStatus :: TwoFluidMaterialStatus(int n, Domain *d, GaussPoint *gp, const IntArray &slaveMaterial) :
    FluidDynamicMaterialStatus(n, d, gp)
{
    MaterialMode mmode = gp->giveMaterialMode();
    int _size = 0;

    if ( mmode == _2dFlow ) {
        _size = 3;
    } else if ( mmode == _2dAxiFlow ) {
        _size = 4;
    } else if ( mmode == _3dFlow ) {
        _size = 6;
    }

    deviatoricStressVector.resize(_size);
    deviatoricStressVector.zero();
    deviatoricStrainRateVector.resize(_size);
    deviatoricStrainRateVector.zero();

    this->slaveGp0 = new GaussPoint(NULL, 0, NULL, 0., mmode);
    this->slaveGp1 = new GaussPoint(NULL, 0, NULL, 0., mmode);
    this->slaveGp0->setMaterialStatus( domain->giveMaterial( slaveMaterial(0) )->CreateStatus(this->slaveGp0), this->giveNumber() );
    this->slaveGp1->setMaterialStatus( domain->giveMaterial( slaveMaterial(1) )->CreateStatus(this->slaveGp0), this->giveNumber() );
}


void
TwoFluidMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    this->giveSlaveGaussPoint0()->giveMaterialStatus()->printOutputAt(file, tStep);
    this->giveSlaveGaussPoint1()->giveMaterialStatus()->printOutputAt(file, tStep);
}


void
TwoFluidMaterialStatus :: updateYourself(TimeStep *tStep)
{
    FluidDynamicMaterialStatus :: updateYourself(tStep);
    this->giveSlaveGaussPoint0()->giveMaterialStatus()->updateYourself(tStep);
    this->giveSlaveGaussPoint1()->giveMaterialStatus()->updateYourself(tStep);
}


void
TwoFluidMaterialStatus :: initTempStatus()
{
    FluidDynamicMaterialStatus :: initTempStatus();
    static_cast< MaterialStatus* >( this->giveSlaveGaussPoint0()->giveMaterialStatus() )->initTempStatus();
    static_cast< MaterialStatus* >( this->giveSlaveGaussPoint1()->giveMaterialStatus() )->initTempStatus();
}


contextIOResultType
TwoFluidMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    this->giveSlaveGaussPoint0()->giveMaterialStatus()->saveContext(stream, mode, obj);
    this->giveSlaveGaussPoint1()->giveMaterialStatus()->saveContext(stream, mode, obj);
    return CIO_OK;
}


contextIOResultType
TwoFluidMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    this->giveSlaveGaussPoint0()->giveMaterialStatus()->restoreContext(stream, mode, obj);
    this->giveSlaveGaussPoint1()->giveMaterialStatus()->restoreContext(stream, mode, obj);
    return CIO_OK;
}
} // end namespace oofem
