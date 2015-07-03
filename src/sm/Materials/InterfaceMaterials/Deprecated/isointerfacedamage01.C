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

#include "isointerfacedamage01.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(IsoInterfaceDamageMaterial);

IsoInterfaceDamageMaterial :: IsoInterfaceDamageMaterial(int n, Domain *d) : StructuralInterfaceMaterial(n, d)
{
    maxOmega = 0.999999;
    beta = 0.;
}


IsoInterfaceDamageMaterial :: ~IsoInterfaceDamageMaterial() { }


void
IsoInterfaceDamageMaterial :: giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
{
    IsoInterfaceDamageMaterialStatus *status = static_cast< IsoInterfaceDamageMaterialStatus * >( this->giveStatus(gp) );
    FloatMatrix de;
    double equivStrain, tempKappa = 0.0, omega = 0.0;

    // compute equivalent strain
    this->computeEquivalentStrain(equivStrain, jump, gp, tStep);

    // compute value of loading function if strainLevel crit apply
    tempKappa = status->giveKappa();

    if ( tempKappa >= equivStrain ) {
        // damage does not grow
        omega = status->giveDamage();
    } else {
        // damage grows
        tempKappa = equivStrain;
        // evaluate damage
        this->computeDamageParam(omega, tempKappa, jump, gp);
    }

    this->give3dStiffnessMatrix_Eng(de, ElasticStiffness, gp, tStep);
    de.times(1.0 - omega);
    answer.beProductOf(de, jump);

    // update gp
    status->letTempJumpBe(jump);
    status->letTempTractionBe(answer);
    status->setTempKappa(tempKappa);
    status->setTempDamage(omega);
}


void
IsoInterfaceDamageMaterial :: give3dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    double om, un;
    IsoInterfaceDamageMaterialStatus *status = static_cast< IsoInterfaceDamageMaterialStatus * >( this->giveStatus(gp) );

    // assemble eleastic stiffness
    answer.resize(3, 3);
    answer.zero();
    answer.at(1, 1) = kn;
    answer.at(2, 2) = ks;
    answer.at(3, 3) = ks;

    if ( rMode == ElasticStiffness ) {
        return;
    }

    if ( rMode == SecantStiffness ) {
        // Secant stiffness
        om = status->giveTempDamage();
        un = status->giveTempJump().at(1);
        om = min(om, maxOmega);
        // damage in tension only
        if ( un >= 0 ) {
            answer.times(1.0 - om);
        }

        return;
    } else {
        // Tangent Stiffness
        FloatArray se;
        const FloatArray &e = status->giveTempJump();
        se.beProductOf(answer, e);

        om = status->giveTempDamage();
        un = status->giveTempJump().at(1);
        om = min(om, maxOmega);
        // damage in tension only
        if ( un >= 0 ) {
            answer.times(1.0 - om);
            return;
#if 0
            // Unreachable code - commented out to supress compiler warnings
            double dom = -( -e0 / un / un * exp( -( ft / gf ) * ( un - e0 ) ) + e0 / un * exp( -( ft / gf ) * ( un - e0 ) ) * ( -( ft / gf ) ) );
            if ( ( om > 0. ) && ( status->giveTempKappa() > status->giveKappa() ) ) {
                answer.at(1, 1) -= se.at(1) * dom;
                answer.at(2, 1) -= se.at(2) * dom;
            }
#endif
        }
    }
}


int
IsoInterfaceDamageMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    IsoInterfaceDamageMaterialStatus *status = static_cast< IsoInterfaceDamageMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_DamageTensor ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveDamage();
        return 1;
    } else if ( type == IST_DamageTensorTemp ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveTempDamage();
        return 1;
    } else if ( type == IST_PrincipalDamageTensor ) {
        answer.resize(3);
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveDamage();
        return 1;
    } else if ( type == IST_PrincipalDamageTempTensor ) {
        answer.resize(3);
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveTempDamage();
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = status->giveKappa();
        return 1;
    } else {
        return StructuralInterfaceMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}


IRResultType
IsoInterfaceDamageMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, kn, _IFT_IsoInterfaceDamageMaterial_kn);
    IR_GIVE_FIELD(ir, ks, _IFT_IsoInterfaceDamageMaterial_ks);

    IR_GIVE_FIELD(ir, ft, _IFT_IsoInterfaceDamageMaterial_ft);
    IR_GIVE_FIELD(ir, gf, _IFT_IsoInterfaceDamageMaterial_gf);
    this->e0 = ft / kn;

    //Set limit on the maximum isotropic damage parameter if needed
    IR_GIVE_OPTIONAL_FIELD(ir, maxOmega, _IFT_IsoInterfaceDamageMaterial_maxOmega);
    maxOmega = min(maxOmega, 0.999999);
    maxOmega = max(maxOmega, 0.0);

    beta = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, beta, _IFT_IsoInterfaceDamageMaterial_beta);

    return StructuralInterfaceMaterial :: initializeFrom(ir);
}


void
IsoInterfaceDamageMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);

    input.setField(this->kn, _IFT_IsoInterfaceDamageMaterial_kn);
    input.setField(this->ks, _IFT_IsoInterfaceDamageMaterial_ks);

    input.setField(this->ft, _IFT_IsoInterfaceDamageMaterial_ft);
    input.setField(this->gf, _IFT_IsoInterfaceDamageMaterial_gf);

    input.setField(this->maxOmega, _IFT_IsoInterfaceDamageMaterial_maxOmega);
}


void
IsoInterfaceDamageMaterial :: computeEquivalentStrain(double &kappa, const FloatArray &jump, GaussPoint *gp, TimeStep *tStep)
{
    double epsNplus = macbra( jump.at(1) );
    double epsT2 = jump.at(2) * jump.at(2);
    if ( jump.giveSize() == 3 ) {
        epsT2 += jump.at(3) * jump.at(3);
    }
    kappa = sqrt(epsNplus * epsNplus + beta * epsT2);
}

void
IsoInterfaceDamageMaterial :: computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp)
{
    if ( kappa > this->e0 ) {
        omega = 1.0 - ( this->e0 / kappa ) * exp( -( ft / gf ) * ( kappa - e0 ) );
    } else {
        omega = 0.0;
    }
}


IsoInterfaceDamageMaterialStatus :: IsoInterfaceDamageMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralInterfaceMaterialStatus(n, d, g)
{
    kappa = tempKappa = 0.0;
    damage = tempDamage = 0.0;
}


IsoInterfaceDamageMaterialStatus :: ~IsoInterfaceDamageMaterialStatus()
{ }


void
IsoInterfaceDamageMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( this->damage > 0.0 ) {
        fprintf(file, "kappa %f, damage %f ", this->kappa, this->damage);
    }

    fprintf(file, "}\n");
}


void
IsoInterfaceDamageMaterialStatus :: initTempStatus()
{
    StructuralInterfaceMaterialStatus :: initTempStatus();
    this->tempKappa = this->kappa;
    this->tempDamage = this->damage;
}

void
IsoInterfaceDamageMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus :: updateYourself(tStep);
    this->kappa = this->tempKappa;
    this->damage = this->tempDamage;
}


contextIOResultType
IsoInterfaceDamageMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralInterfaceMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
IsoInterfaceDamageMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralInterfaceMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}
} // end namespace oofem
