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

#include "intmatisodamage.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(IntMatIsoDamage);

IntMatIsoDamage :: IntMatIsoDamage(int n, Domain *d) : StructuralInterfaceMaterial(n, d)
    //
    // constructor
    //
{
    maxOmega = 0.999999;
}


IntMatIsoDamage :: ~IntMatIsoDamage()
//
// destructor
//
{ }


void
IntMatIsoDamage :: giveEngTraction_3d(FloatArray &answer, GaussPoint *gp,
                                                   const FloatArray &jump,
                                                   TimeStep *tStep)
//
// returns 3d traction vector 
{
    IntMatIsoDamageStatus *status = static_cast< IntMatIsoDamageStatus * >( this->giveStatus(gp) );
    
    double f, equivJump, tempKappa = 0.0, omega = 0.0;

    this->initGpForNewStep(gp);

    // compute equivalent strain
    this->computeEquivalentJump(equivJump, jump);

    // compute value of loading function 
    f = equivJump - status->giveKappa();

    if ( f <= 0.0 ) {
        // no damage evolution
        tempKappa = status->giveKappa();
        omega     = status->giveDamage();
    } else {
        // damage evolution
        tempKappa = equivJump;
        // evaluate damage parameter
        this->computeDamageParam(omega, tempKappa);
    }
    FloatMatrix D;
    this->give3dStiffnessMatrix_Eng(D, ElasticStiffness, gp, tStep);
    // damage in tension only
    if ( equivJump >= 0.0 ) {
        D.times(1.0 - omega);
    }

    answer.beProductOf(D, jump);

    // update gp
    status->letTempJumpBe(jump);
    status->letTempTractionBe(answer);

    status->setTempKappa(tempKappa);
    status->setTempDamage(omega);
}



void
IntMatIsoDamage :: give2dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode,
                                                                     GaussPoint *gp, TimeStep *tStep)
{
    double om, un;
    IntMatIsoDamageStatus *status = static_cast< IntMatIsoDamageStatus * >( this->giveStatus(gp) );


    if ( ( rMode == ElasticStiffness ) || ( rMode == SecantStiffness ) || ( rMode == TangentStiffness ) ) {
        // assemble eleastic stiffness
        answer.resize(2, 2);
        answer.zero();
        answer.at(1, 1) = ks;
        answer.at(2, 2) = kn;

        if ( rMode == ElasticStiffness ) {
            return;
        }

        if ( rMode == SecantStiffness ) {
            // Secant stiffness
            om = status->giveTempDamage();
            un = status->giveTempJump().at(2);
            om = min(om, maxOmega);
            // damage in tension only
            if ( un >= 0 ) {
                answer.times(1.0 - om);
            }

            return;
        } else {
            // Tangent Stiffness
            FloatArray se(2), jump3d, jump2d(2);

            jump3d = status->giveTempJump();
            jump2d.at(1) = jump3d.at(1);
            jump2d.at(2) = jump3d.at(3);
            se.beProductOf(answer, jump2d);

            om = status->giveTempDamage();
            un = jump2d.at(2);
            om = min(om, maxOmega);
            // damage in tension only
            if ( un >= 0 ) {
                answer.times(1.0 - om); ///@todo this is only the secant stiffness - tangent stiffness is broken!

                double omega, omega_plus;
                computeDamageParam(omega, un);
                computeDamageParam(omega_plus, un + 1.0e-8 );
                //double dom = (omega_plus - omega)/ 1.0e-8;
                
                // d( (1-omega)*D*j ) / dj = (1-omega)D - D*j openprod domega/dj
            /*    double fac = ft*(e0 - un)/gf;
                dom = e0*exp(fac) /(un*un + 1.0e-9) + e0*ft*exp(fac) / (gf*un + 1.0e-9);
                
                
                dom = -( -e0 / (un * un+1.0e-9) * exp( -( ft / gf ) * ( un - e0 ) ) + e0 / (un+1.0e-9) * exp( -( ft / gf ) * ( un - e0 ) ) * ( -( ft / gf ) ) );
                if ( ( om > 0. ) && ( status->giveTempKappa() > status->giveKappa() ) ) {
                    answer.at(1, 2) -= se.at(1) * dom;
                    answer.at(2, 2) -= se.at(2) * dom;
                }*/
                return;
            }
        }
    }  else {
        OOFEM_ERROR("Unknown MatResponseMode");
    }
}


void
IntMatIsoDamage :: give3dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode,
                                                                     GaussPoint *gp, TimeStep *tStep)
{
    double om, un;
    IntMatIsoDamageStatus *status = static_cast< IntMatIsoDamageStatus * >( this->giveStatus(gp) );


    if ( ( rMode == ElasticStiffness ) || ( rMode == SecantStiffness ) || ( rMode == TangentStiffness ) ) {
        // assemble eleastic stiffness
        answer.resize(3, 3);
        answer.zero();
        answer.at(1, 1) = ks;
        answer.at(2, 2) = ks;
        answer.at(3, 3) = kn;

        if ( rMode == ElasticStiffness ) {
            return;
        }

        if ( rMode == SecantStiffness ) {
            // Secant stiffness
            om = status->giveTempDamage();
            //un = status->giveTempStrainVector().at(1);
            un = status->giveTempJump().at(3);
            om = min(om, maxOmega);
            // damage in tension only
            if ( un >= 0 ) {
                answer.times(1.0 - om);
            }

            return;
        } else {
            // Tangent Stiffness
            FloatArray se, jump;
            //e = status->giveTempStrainVector();
            jump = status->giveTempJump();
            se.beProductOf(answer, jump);

            om = status->giveTempDamage();
            //un = status->giveTempStrainVector().at(1);
            un = jump.at(3);
            om = min(om, maxOmega);
            // damage in tension only
            if ( un >= 0 ) {
                answer.times(1.0 - om);
                return;
            }
        }
    }  else {
        OOFEM_ERROR("Unknown MatResponseMode");
    }
}


int
IntMatIsoDamage :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    IntMatIsoDamageStatus *status = static_cast< IntMatIsoDamageStatus * >( this->giveStatus(gp) );
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
IntMatIsoDamage :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, kn, _IFT_IntMatIsoDamage_kn);
    IR_GIVE_FIELD(ir, ks, _IFT_IntMatIsoDamage_ks);

    IR_GIVE_FIELD(ir, ft, _IFT_IntMatIsoDamage_ft);
    IR_GIVE_FIELD(ir, gf, _IFT_IntMatIsoDamage_gf);
    this->e0 = ft / kn;

    //Set limit on the maximum isotropic damage parameter if needed
    IR_GIVE_OPTIONAL_FIELD(ir, maxOmega, _IFT_IntMatIsoDamage_maxOmega);
    maxOmega = min(maxOmega, 0.999999);
    maxOmega = max(maxOmega, 0.0);

    return StructuralInterfaceMaterial :: initializeFrom(ir);
}


void
IntMatIsoDamage :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);

    input.setField(this->kn, _IFT_IntMatIsoDamage_kn);
    input.setField(this->ks, _IFT_IntMatIsoDamage_ks);

    input.setField(this->ft, _IFT_IntMatIsoDamage_ft);
    input.setField(this->gf, _IFT_IntMatIsoDamage_gf);

    input.setField(this->maxOmega, _IFT_IntMatIsoDamage_maxOmega);
}


void
IntMatIsoDamage :: computeEquivalentJump(double &kappa, const FloatArray &jump)
{
    // Damage only driven by the normal jump

    if ( jump.giveSize() == 3 ) {
        kappa = macbra( jump.at(3) );
    } else if ( jump.giveSize() == 2 ) {
        kappa = macbra( jump.at(2) );
    } else if( jump.giveSize() == 1 ) {
        kappa = macbra( jump.at(1) );
    } else {
        OOFEM_ERROR("Jump vector is not of dimension 1,2 or 3");
    }
}

void
IntMatIsoDamage :: computeDamageParam(double &omega, double kappa)
{
    if ( kappa > this->e0 ) {
        omega = 1.0 - ( this->e0 / kappa ) * exp( -( ft / gf ) * ( kappa - e0 ) );
    } else {
        omega = 0.0;
    }
}


IntMatIsoDamageStatus :: IntMatIsoDamageStatus(int n, Domain *d, GaussPoint *g) : StructuralInterfaceMaterialStatus(n, d, g)
{
    kappa  = tempKappa  = 0.0;
    damage = tempDamage = 0.0;
}


IntMatIsoDamageStatus :: ~IntMatIsoDamageStatus()
{ }


void
IntMatIsoDamageStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( this->damage > 0.0 ) {
        fprintf(file, "kappa %f, damage %f ", this->kappa, this->damage);
    }

    fprintf(file, "}\n");
}


void
IntMatIsoDamageStatus :: initTempStatus()
{
    StructuralInterfaceMaterialStatus :: initTempStatus();
    this->tempKappa = this->kappa;
    this->tempDamage = this->damage;
}

void
IntMatIsoDamageStatus :: updateYourself(TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus :: updateYourself(tStep);
    this->kappa = this->tempKappa;
    this->damage = this->tempDamage;
}


contextIOResultType
IntMatIsoDamageStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralInterfaceMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream->write(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
IntMatIsoDamageStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralInterfaceMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream->read(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}
} // end namespace oofem
