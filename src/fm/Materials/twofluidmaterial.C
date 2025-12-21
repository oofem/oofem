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

#include "fm/Materials/twofluidmaterial.h"
#include "domain.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "engngm.h"
#include "materialinterface.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "gausspoint.h"

namespace oofem {
REGISTER_Material(TwoFluidMaterial);

int
TwoFluidMaterial :: checkConsistency()
{
    return this->giveMaterial(0)->checkConsistency() &&
           this->giveMaterial(1)->checkConsistency();
}


void
TwoFluidMaterial :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_FIELD(ir, this->slaveMaterial, _IFT_TwoFluidMaterial_mat);
    if ( this->slaveMaterial.giveSize() != 2 ) {
        throw ValueInputException(ir, _IFT_TwoFluidMaterial_mat, "mat array should have two values");
    }
}


void
TwoFluidMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    FluidDynamicMaterial :: giveInputRecord(input);
    input.setField(this->slaveMaterial, _IFT_TwoFluidMaterial_mat);
}


double
TwoFluidMaterial :: giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep) const
{
    TwoFluidMaterialStatus *status = static_cast< TwoFluidMaterialStatus * >( this->giveStatus(gp) );
    double vof = this->giveTempVOF(gp);
    return ( 1.0 - vof ) * giveMaterial(0)->giveEffectiveViscosity(status->giveSlaveGaussPoint0(), tStep) +
           vof *giveMaterial(1)->giveEffectiveViscosity(status->giveSlaveGaussPoint1(), tStep);
}


double
TwoFluidMaterial :: give(int aProperty, GaussPoint *gp) const
{
    TwoFluidMaterialStatus *status = static_cast< TwoFluidMaterialStatus * >( this->giveStatus(gp) );
    double vof = this->giveTempVOF(gp);
    return ( 1.0 - vof ) * giveMaterial(0)->give( aProperty, status->giveSlaveGaussPoint0() ) +
           vof *giveMaterial(1)->give( aProperty, status->giveSlaveGaussPoint1() );
}

int
TwoFluidMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    TwoFluidMaterialStatus *status = static_cast< TwoFluidMaterialStatus * >( this->giveStatus(gp) );
    double vof = this->giveTempVOF(gp);
    FloatArray tmp;
    int ret = giveMaterial(0)->giveIPValue(answer, status->giveSlaveGaussPoint0(), type, tStep);
    answer.times(1.0 - vof);
    ret = ret && giveMaterial(1)->giveIPValue(tmp, status->giveSlaveGaussPoint1(), type, tStep);
    answer.add(vof, tmp);
    return ret;
}


std::unique_ptr<MaterialStatus> 
TwoFluidMaterial :: CreateStatus(GaussPoint *gp) const
{
    std::array<Material*, 2> slaveMaterial({this->giveMaterial(0), this->giveMaterial(1)});
    return std::make_unique<TwoFluidMaterialStatus>(gp, slaveMaterial);
}

FluidDynamicMaterial *
TwoFluidMaterial :: giveMaterial(int i) const
{
    return static_cast< FluidDynamicMaterial * >( domain->giveMaterial( slaveMaterial[i] ) );
}

FloatArrayF<6>
TwoFluidMaterial :: computeDeviatoricStress3D(const FloatArrayF<6> &eps, GaussPoint *gp, TimeStep *tStep) const
{
    double vof = this->giveTempVOF(gp);
    TwoFluidMaterialStatus *status = static_cast< TwoFluidMaterialStatus * >( this->giveStatus(gp) );

    auto v0 = this->giveMaterial(0)->computeDeviatoricStress3D(eps, status->giveSlaveGaussPoint0(), tStep);
    auto v1 = this->giveMaterial(1)->computeDeviatoricStress3D(eps, status->giveSlaveGaussPoint1(), tStep);

    auto stress = (1.0 - vof) * v0 + vof * v1;

    status->letDeviatoricStrainRateVectorBe(eps);
    status->letDeviatoricStressVectorBe(stress);

    return stress;
}

FloatMatrixF<6,6>
TwoFluidMaterial :: computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    double vof = this->giveTempVOF(gp);
    TwoFluidMaterialStatus *status = static_cast< TwoFluidMaterialStatus * >( this->giveStatus(gp) );

    auto a0 = this->giveMaterial(0)->computeTangent3D(mode, status->giveSlaveGaussPoint0(), tStep);
    auto a1 = this->giveMaterial(1)->computeTangent3D(mode, status->giveSlaveGaussPoint1(), tStep);

    return (1.0 - vof) * a0 + vof * a1;
}

double
TwoFluidMaterial :: giveTempVOF(GaussPoint *gp) const
{
    FloatArray vof(2);
    MaterialInterface *mi = domain->giveEngngModel()->giveMaterialInterface( domain->giveNumber() );
    if ( mi ) {
        mi->giveElementMaterialMixture( vof, gp->giveElement()->giveNumber() );

        if ( ( vof.at(1) < 0. ) || ( vof.at(1) > 1.0 ) ) {
            OOFEM_ERROR("vof value out of range (vof=%lf)", vof.at(1));
        }

        return vof.at(1);
    } else {
        return 0.0;
    }
}




TwoFluidMaterialStatus :: TwoFluidMaterialStatus(GaussPoint *gp, const std::array<Material*, 2> &slaveMaterial) :
    FluidDynamicMaterialStatus(gp),
    slaveGps{{{nullptr, 0, 0., gp->giveMaterialMode()}, {nullptr, 0, 0., gp->giveMaterialMode()}}}
{
    for ( int i = 0; i < 2; ++i ) slaveGps[i].setMaterialStatus( slaveMaterial[i]->CreateStatus( &slaveGps[i] ) );
}


void
TwoFluidMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    for ( auto &gp : slaveGps ) gp.giveMaterialStatus()->printOutputAt(file, tStep);
}


void
TwoFluidMaterialStatus :: updateYourself(TimeStep *tStep)
{
    FluidDynamicMaterialStatus :: updateYourself(tStep);
    for ( auto &gp : slaveGps ) gp.giveMaterialStatus()->updateYourself(tStep);
}


void
TwoFluidMaterialStatus :: initTempStatus()
{
    FluidDynamicMaterialStatus :: initTempStatus();
    for ( auto &gp : slaveGps ) static_cast< MaterialStatus * >( gp.giveMaterialStatus() )->initTempStatus();
}


void
TwoFluidMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    for ( auto &gp : slaveGps ) gp.giveMaterialStatus()->saveContext(stream, mode);
}


void
TwoFluidMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    for ( auto &gp : slaveGps ) gp.giveMaterialStatus()->restoreContext(stream, mode);
}
} // end namespace oofem
