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

#include "expczmaterial.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
const double tolerance = 1.0e-12; // small number
ExpCZMaterial :: ExpCZMaterial(int n, Domain *d) : StructuralMaterial(n, d)
    //
    // constructor
    //
{ }


ExpCZMaterial :: ~ExpCZMaterial()
//
// destructor
//
{ }

int
ExpCZMaterial :: hasMaterialModeCapability(MaterialMode mode)
{
    // returns whether receiver supports given mode
    if ( mode == _3dInterface ) {
        return 1;
    } else {
        return 0;
    }
}




void
ExpCZMaterial :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                      const FloatArray &jumpVector,
                                      TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    ExpCZMaterialStatus *status = static_cast< ExpCZMaterialStatus * >( this->giveStatus(gp) );

    //this->initGpForNewStep(gp);
    this->initTempStatus(gp);

    answer.resize( jumpVector.giveSize() );
    answer.zero();
    //@todo for now only study normal stress
    // no degradation in shear
    // jumpVector = [shear 1 shear 2 normal]
    double gn  = jumpVector.at(3);
    double gs1 = jumpVector.at(1);
    double gs2 = jumpVector.at(2);
    double gs  = sqrt(gs1 * gs1 + gs2 * gs2);

    double xin  = gn / ( gn0 + tolerance );
    double xit  = gs / ( gs0 + tolerance );
    double tn   = GIc / gn0 *exp(-xin) * ( xin * exp(-xit * xit) + ( 1.0 - q ) / ( r - 1.0 ) * ( 1.0 - exp(-xit * xit) ) * ( r - xin ) );

    answer.at(1) = 0.0;
    answer.at(2) = 0.0;
    answer.at(3) = tn;


    // update gp
    status->letTempStrainVectorBe(jumpVector);
    status->letTempStressVectorBe(answer);
}

void
ExpCZMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                          MatResponseForm form, MatResponseMode rMode,
                                          GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dInterface:
    case _3dMat:
        give3dInterfaceMaterialStiffnessMatrix(answer, form, rMode, gp, tStep);
        break;
    default:
        StructuralMaterial :: giveCharacteristicMatrix(answer, form, rMode, gp, tStep);
    }
}


void
ExpCZMaterial :: give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode,
                                                        GaussPoint *gp, TimeStep *tStep)
{
    ExpCZMaterialStatus *status = static_cast< ExpCZMaterialStatus * >( this->giveStatus(gp) );

    FloatArray jumpVector;
    jumpVector = status->giveTempStrainVector();
    //@todo for now only study normal stress
    // jumpVector = [shear 1 shear 2 normal]
    double gn  = jumpVector.at(3);
    double gs1 = jumpVector.at(1);
    double gs2 = jumpVector.at(2);
    double gs  = sqrt(gs1 * gs1 + gs2 * gs2);

    if ( rMode == TangentStiffness ) {
        answer.resize(3, 3);
        answer.zero();

        double xin  = gn / ( gn0 + tolerance );
        double xit  = gs / ( gs0 + tolerance );
        double dtndgn = GIc / gn0 / gn0 *exp(-xin) *
        ( ( 1.0 - xin ) * exp(-xit * xit) + ( 1.0 - q ) / ( r - 1.0 ) * ( 1.0 - exp(-xit * xit) ) * ( -r + xin - 1.0 ) );

        answer.at(3, 3) = dtndgn;
    }  else {
        OOFEM_ERROR("unknown MatResponseMode");
    }
}



#if 0
int
ExpCZMaterial :: giveSizeOfReducedStressStrainVector(MaterialMode mode)
//
// returns the size of reduced stress-strain vector
// according to mode given by gp.
//
{
    switch ( mode ) {
    case _2dInterface:
        return 2;

    case _3dInterface:
        return 3;

    default:
        return StructuralMaterial :: giveSizeOfReducedStressStrainVector(mode);
    }
}


int
ExpCZMaterial :: giveStressStrainComponentIndOf(MatResponseForm form, MaterialMode mmode, int ind)
//
// this function returns index of reduced(if form == ReducedForm)
// or Full(if form==FullForm) stressStrain component in Full or reduced
// stressStrainVector acording to stressStrain mode of given gp.
//
{
    //MaterialMode mode  = gp -> giveMaterialMode ();

    if ( ( mmode == _2dInterface ) || ( mmode == _3dInterface ) ) {
        return ind;
    } else {
        return StructuralMaterial :: giveStressStrainComponentIndOf(form, mmode, ind);
    }
}


void
ExpCZMaterial :: giveStressStrainMask(IntArray &answer, MatResponseForm form,
                                      MaterialMode mmode) const
//
// this function returns mask of reduced(if form == ReducedForm)
// or Full(if form==FullForm) stressStrain vector in full or
// reduced StressStrainVector
// acording to stressStrain mode of given gp.
//
//
// mask has size of reduced or full StressStrain Vector and  i-th component
// is index to full or reduced StressStrainVector where corresponding
// stressStrain resides.
//
// Reduced form is sub-vector (of stress or strain components),
// where components corresponding to imposed zero stress (plane stress,...)
// are not included. On the other hand, if zero strain component is imposed
// (Plane strain, ..) this condition must be taken into account in geometrical
// relations, and corresponding component is included in reduced vector.
//
{
    if ( mmode == _2dInterface ) {
        answer.resize(2);
        for ( int i = 1; i <= 2; i++ ) {
            answer.at(i) = i;
        }
    } else if ( mmode == _3dInterface ) {
        answer.resize(3);
        for ( int i = 1; i <= 3; i++ ) {
            answer.at(i) = i;
        }
    } else {
        StructuralMaterial :: giveStressStrainMask(answer, form, mmode);
    }
}


void
ExpCZMaterial :: giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                                 const FloatArray &charVector3d)
//
// returns reduced stressVector or strainVector from full 3d vector reduced
// to vector required by gp->giveStressStrainMode()
//
{
    MaterialMode mode = gp->giveMaterialMode();

    if ( ( mode == _2dInterface ) || ( mode == _3dInterface ) ) {
        answer = charVector3d;
        return;
    } else {
        StructuralMaterial :: giveReducedCharacteristicVector(answer, gp, charVector3d);
    }
}


void
ExpCZMaterial :: giveFullCharacteristicVector(FloatArray &answer,
                                              GaussPoint *gp,
                                              const FloatArray &strainVector)
//
// returns full 3d general strain vector from strainVector in reducedMode
// based on StressStrainMode in gp. Included are strains which
// perform nonzero work.
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
// 2) strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
// you must assigng your stress strain mode to one of the folloving modes (or add new)
// FullForm of MaterialStiffnessMatrix must have the same form.
//
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( ( mode == _2dInterface ) || ( mode == _3dInterface ) ) {
        answer = strainVector;
        return;
    } else {
        StructuralMaterial :: giveFullCharacteristicVector(answer, gp, strainVector);
    }
}


void
ExpCZMaterial :: give2dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode,
                                                        GaussPoint *gp, TimeStep *tStep)
{
    double om, un;
    IsoInterfaceDamageMaterialStatus *status = static_cast< IsoInterfaceDamageMaterialStatus * >( this->giveStatus(gp) );


    if ( ( rMode == ElasticStiffness ) || ( rMode == SecantStiffness ) || ( rMode == TangentStiffness ) ) {
        // assemble eleastic stiffness
        answer.resize(2, 2);
        answer.at(1, 1) = kn;
        answer.at(2, 2) = ks;
        answer.at(1, 2) = answer.at(2, 1) = 0.0;

        if ( rMode == ElasticStiffness ) {
            return;
        }

        if ( rMode == SecantStiffness ) {
            // Secant stiffness
            om = status->giveTempDamage();
            un = status->giveTempStrainVector().at(1);
            om = min(om, maxOmega);
            // damage in tension only
            if ( un >= 0 ) {
                answer.times(1.0 - om);
            }

            return;
        } else {
            // Tangent Stiffness
            FloatArray se(2), e(2);
            e = status->giveTempStrainVector();
            se.beProductOf(answer, e);

            om = status->giveTempDamage();
            un = status->giveTempStrainVector().at(1);
            om = min(om, maxOmega);
            // damage in tension only
            if ( un >= 0 ) {
                answer.times(1.0 - om);
                return;

                /* Unreachable code - commented out to supress compiler warnings
                 * double dom = -( -e0 / un / un * exp( -( ft / gf ) * ( un - e0 ) ) + e0 / un * exp( -( ft / gf ) * ( un - e0 ) ) * ( -( ft / gf ) ) );
                 * if ( ( om > 0. ) && ( status->giveTempKappa() > status->giveKappa() ) ) {
                 *  answer.at(1, 1) -= se.at(1) * dom;
                 *  answer.at(2, 1) -= se.at(2) * dom;
                 * }
                 */
            }
        }
    }  else {
        OOFEM_ERROR("unknown MatResponseMode");
    }
}


void
ExpCZMaterial :: give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode,
                                                        GaussPoint *gp, TimeStep *tStep)
{
    double om, un;
    IsoInterfaceDamageMaterialStatus *status = static_cast< IsoInterfaceDamageMaterialStatus * >( this->giveStatus(gp) );


    if ( ( rMode == ElasticStiffness ) || ( rMode == SecantStiffness ) || ( rMode == TangentStiffness ) ) {
        // assemble eleastic stiffness
        answer.resize(3, 3);
        answer.at(1, 1) = kn;
        answer.at(2, 2) = ks;
        answer.at(3, 3) = ks;
        answer.at(1, 2) = answer.at(2, 1) = answer.at(1, 3) = answer.at(3, 1) = answer.at(2, 3) = answer.at(3, 2) = 0.0;

        if ( rMode == ElasticStiffness ) {
            return;
        }

        if ( rMode == SecantStiffness ) {
            // Secant stiffness
            om = status->giveTempDamage();
            un = status->giveTempStrainVector().at(1);
            om = min(om, maxOmega);
            // damage in tension only
            if ( un >= 0 ) {
                answer.times(1.0 - om);
            }

            return;
        } else {
            // Tangent Stiffness
            FloatArray se, e;
            e = status->giveTempStrainVector();
            se.beProductOf(answer, e);

            om = status->giveTempDamage();
            un = status->giveTempStrainVector().at(1);
            om = min(om, maxOmega);
            // damage in tension only
            if ( un >= 0 ) {
                answer.times(1.0 - om);
                return;
                /* Unreachable code - commented out to supress compiler warnings
                 * double dom = -( -e0 / un / un * exp( -( ft / gf ) * ( un - e0 ) ) + e0 / un * exp( -( ft / gf ) * ( un - e0 ) ) * ( -( ft / gf ) ) );
                 * if ( ( om > 0. ) && ( status->giveTempKappa() > status->giveKappa() ) ) {
                 *  answer.at(1, 1) -= se.at(1) * dom;
                 *  answer.at(2, 1) -= se.at(2) * dom;
                 * }
                 */
            }
        }
    }  else {
        OOFEM_ERROR("unknown MatResponseMode");
    }
}


#endif

int
ExpCZMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    ExpCZMaterialStatus *status = static_cast< ExpCZMaterialStatus * >( this->giveStatus(gp) );
    return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
}

int
ExpCZMaterial :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
}


int
ExpCZMaterial :: giveIPValueSize(InternalStateType type, GaussPoint *gp)
{
    return StructuralMaterial :: giveIPValueSize(type, gp);
}




IRResultType
ExpCZMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->GIc, _IFT_ExpCZMaterial_g1c);

    this->sigfs = 0.0;
    IR_GIVE_FIELD(ir, this->sigfn, _IFT_ExpCZMaterial_sigfn);

    GIIc = 0.0;

    this->gn0 = GIc  / ( this->sigfn * exp(1.0) + tolerance );                // normal jump at max stress
    this->gs0 = GIIc / ( sqrt( 0.5 * exp(1.0) ) * this->sigfs + tolerance );      // shear jump at max stress


    this->q = GIIc / GIc;
    this->r = 0.0; // fix

    //this->checkConsistency();                                // check validity of the material paramters
    this->printYourself();
    return IRRT_OK;
}

int
ExpCZMaterial :: checkConsistency()
{
    if ( this->GIc < 0.0 ) {
        OOFEM_ERROR("GIc is negative (%.2e)", this->GIc);
    }
    return 1;
}

void
ExpCZMaterial :: printYourself()
{
    printf("Paramters for ExpCZMaterial: \n");

    printf("-Strength paramters \n");
    printf("  sigfn = %e \n", this->sigfn);
    printf("  GIc   = %e \n", this->GIc);
    printf("  GIIc  = %e \n", this->GIIc);
    //printf("\n");

    printf("-jump limits \n");
    printf("  gn0   = %e \n", this->gn0);
    printf("  gs0   = %e \n", this->gs0);

    printf("-other parameters \n");
    printf("  q     = %e \n", this->q);
    printf("  r     = %e \n", this->r);
}

ExpCZMaterialStatus :: ExpCZMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{ }


ExpCZMaterialStatus :: ~ExpCZMaterialStatus()
{ }


void
ExpCZMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    /*
     * fprintf(file, "status { ");
     * if ( this->damage > 0.0 ) {
     *  fprintf(file, "kappa %f, damage %f ", this->kappa, this->damage);
     * }
     *
     * fprintf(file, "}\n");
     */
}


void
ExpCZMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    //this->tempKappa = this->kappa;
    //this->tempDamage = this->damage;
}

void
ExpCZMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
}

#if 0
contextIOResultType
ExpCZMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
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
ExpCZMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
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
#endif
} // end namespace oofem
