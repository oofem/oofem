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

#include "bilinearczmaterial.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {

BilinearCZMaterial :: BilinearCZMaterial(int n, Domain *d) : StructuralMaterial(n, d)
//
// constructor
//
{
}


BilinearCZMaterial :: ~BilinearCZMaterial()
//
// destructor
//
{ }

int
BilinearCZMaterial :: hasMaterialModeCapability(MaterialMode mode)
{
    // returns whether receiver supports given mode
    if ( mode == _3dInterface ) {
        return 1;
    } else {
        return 0;
    }
}




void
BilinearCZMaterial :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                                   const FloatArray &jumpVector,
                                                   TimeStep *atTime)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    BilinearCZMaterialStatus *status = static_cast< BilinearCZMaterialStatus * >( this->giveStatus(gp) );
    
    this->initGpForNewStep(gp);


    answer.resize( jumpVector.giveSize() );
    answer.zero();
    //@todo for now only study normal stress
    // no degradation in shear
    // jumpVector = [shear 1 shear 2 normal]
    double gn  = jumpVector.at(3);
    double gs1 = jumpVector.at(1);
    double gs2 = jumpVector.at(2);

    if ( gn <= 0.0 ) { // compresssion
        answer.at(1) = this->ks0 * gs1;
        answer.at(2) = this->ks0 * gs2;
        answer.at(3) = this->knc * gn;
    } else if ( gn <= this->gn0 ) { // first linear part
        answer.at(1) = this->ks0 * gs1;
        answer.at(2) = this->ks0 * gs2;
        answer.at(3) = this->kn0 * gn;
    } else if ( gn <= this->gnmax  ) {  // softening branch
        answer.at(1) = this->ks0 * gs1; 
        answer.at(2) = this->ks0 * gs2;
        answer.at(3) = this->sigfn + this->kn1 * (gn - this->gn0 );
        printf("Softening branch...\n");
    } else {
        answer.at(1) = this->ks0 * gs1; 
        answer.at(2) = this->ks0 * gs2;
        answer.at(3) = 0.0; // failed
        printf("Failed...\n");
    }

    //answer.printYourself();
    // update gp
    status->letTempStrainVectorBe(jumpVector);
    status->letTempStressVectorBe(answer);
    
}

void
BilinearCZMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                       MatResponseForm form, MatResponseMode rMode,
                                                       GaussPoint *gp, TimeStep *atTime)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dInterface:
    case _3dMat:
        give3dInterfaceMaterialStiffnessMatrix(answer, form, rMode, gp, atTime);
        break;
    default:
        StructuralMaterial :: giveCharacteristicMatrix(answer, form, rMode, gp, atTime);
    }
}


void
BilinearCZMaterial :: give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode,
                                                                     GaussPoint *gp, TimeStep *atTime)
{
    BilinearCZMaterialStatus *status = static_cast< BilinearCZMaterialStatus * >( this->giveStatus(gp) );

    FloatArray jumpVector;
    jumpVector = status->giveTempStrainVector();
    //@todo for now only study normal stress
    // jumpVector = [shear 1 shear 2 normal]
    double gn  = jumpVector.at(3);
    double gs1 = jumpVector.at(1);
    double gs2 = jumpVector.at(2);

    //jumpVector.printYourself();

    //if ( ( rMode == ElasticStiffness ) || ( rMode == SecantStiffness ) || ( rMode == TangentStiffness ) ) {
    if ( rMode == TangentStiffness ) {
        // assemble eleastic stiffness
        
        answer.resize(3, 3);
        answer.zero();
        if ( gn <= 0.0 ) { // compresssion
            answer.at(1,1) = this->ks0;
            answer.at(2,2) = this->ks0;
            answer.at(3,3) = this->knc;
        } else if ( gn <= this->gn0 ) { // first linear part
            answer.at(1,1) = this->ks0;
            answer.at(2,2) = this->ks0;
            answer.at(3,3) = this->kn0;
            //printf("Linear branch...\n");
        } else if ( gn <= this->gnmax  ) { // softening branch
            answer.at(1,1) = this->ks0; // no degradation in shear
            answer.at(2,2) = this->ks0;
            answer.at(3,3) = this->kn1;
            //printf("Softening branch...\n");
        } else {
            answer.at(1,1) = this->ks0; 
            answer.at(2,2) = this->ks0;
            answer.at(3,3) = 0.0; // failed
            //printf("Failed...\n");
        }
    
    }  else {
        _error("give2dInterfaceMaterialStiffnessMatrix: unknown MatResponseMode");
    }
    
}



#if 0
int
BilinearCZMaterial :: giveSizeOfReducedStressStrainVector(MaterialMode mode)
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
BilinearCZMaterial :: giveStressStrainComponentIndOf(MatResponseForm form, MaterialMode mmode, int ind)
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
BilinearCZMaterial :: giveStressStrainMask(IntArray &answer, MatResponseForm form,
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
BilinearCZMaterial :: giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
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
BilinearCZMaterial :: giveFullCharacteristicVector(FloatArray &answer,
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
BilinearCZMaterial :: give2dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode,
                                                                     GaussPoint *gp, TimeStep *atTime)
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
                double dom = -( -e0 / un / un * exp( -( ft / gf ) * ( un - e0 ) ) + e0 / un * exp( -( ft / gf ) * ( un - e0 ) ) * ( -( ft / gf ) ) );
                if ( ( om > 0. ) && ( status->giveTempKappa() > status->giveKappa() ) ) {
                    answer.at(1, 1) -= se.at(1) * dom;
                    answer.at(2, 1) -= se.at(2) * dom;
                }
                */
            }
        }
    }  else {
        _error("give2dInterfaceMaterialStiffnessMatrix: unknown MatResponseMode");
    }
}


void
BilinearCZMaterial :: give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode,
                                                                     GaussPoint *gp, TimeStep *atTime)
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
                double dom = -( -e0 / un / un * exp( -( ft / gf ) * ( un - e0 ) ) + e0 / un * exp( -( ft / gf ) * ( un - e0 ) ) * ( -( ft / gf ) ) );
                if ( ( om > 0. ) && ( status->giveTempKappa() > status->giveKappa() ) ) {
                    answer.at(1, 1) -= se.at(1) * dom;
                    answer.at(2, 1) -= se.at(2) * dom;
                }
                */
            }
        }
    }  else {
        _error("give2dInterfaceMaterialStiffnessMatrix: unknown MatResponseMode");
    }
}


#endif

int
BilinearCZMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    BilinearCZMaterialStatus *status = static_cast< BilinearCZMaterialStatus * >( this->giveStatus(aGaussPoint) );
    return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    
}


InternalStateValueType
BilinearCZMaterial :: giveIPValueType(InternalStateType type)
{
    return StructuralMaterial :: giveIPValueType(type);
}


int
BilinearCZMaterial :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
        return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
}


int
BilinearCZMaterial :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    return StructuralMaterial :: giveIPValueSize(type, aGaussPoint);
}



const double tolerance = 1.0e-12; // small number
IRResultType
BilinearCZMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";  // Required by IR_GIVE_FIELD macro
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, kn0, _IFT_BilinearCZMaterial_kn);
    this->knc = kn0;                        // Defaults to the same stiffness in compression and tension
    IR_GIVE_OPTIONAL_FIELD(ir, this->knc, _IFT_BilinearCZMaterial_knc);
    
    this->ks0 = 0.0;                        // Defaults to no shear stiffness
    IR_GIVE_OPTIONAL_FIELD(ir, ks0, _IFT_BilinearCZMaterial_ks);

    IR_GIVE_FIELD(ir, GIc, _IFT_BilinearCZMaterial_g1c);

    this->sigfn = 1.0e50; // defaults to infinite strength @todo but then it will not be bilinear only linear
    this->sigfs = 1.0e50;
    IR_GIVE_OPTIONAL_FIELD(ir, sigfn, _IFT_BilinearCZMaterial_sigfn);

    this->gn0 = sigfn / (kn0 + tolerance);                   // normal jump at damage initiation 
    this->gs0 = sigfs / (ks0 + tolerance);                   // shear jump at damage initiation
    this->gnmax = 2.0 * GIc / sigfn;                         // @todo defaults to zero - will this cause problems?
    this->kn1 = - this->sigfn / ( this->gnmax - this->gn0 ); // slope during softening part in normal dir
    double kn0min = 0.5*sigfn*sigfn/GIc;
    this->checkConsistency();                                // check validity of the material paramters
    return IRRT_OK;
}

int
BilinearCZMaterial :: checkConsistency()
{
    double kn0min = 0.5*this->sigfn*this->sigfn/this->GIc;
    if ( this->kn0 < 0.0 ) {
        OOFEM_ERROR2("BilinearCZMaterial :: initializeFrom - stiffness kn0 is negative (%.2e)", this->kn0);
    } else if ( this->ks0 < 0.0 ) {
        OOFEM_ERROR2("BilinearCZMaterial :: initializeFrom - stiffness ks0 is negative (%.2e)", this->ks0);
    } else if ( this->GIc < 0.0 ) {
        OOFEM_ERROR2("BilinearCZMaterial :: initializeFrom - GIc is negative (%.2e)", this->GIc);
    } else if ( this->kn0 < kn0min  ) { // => gn0 > gnmax
        //OOFEM_ERROR3("BilinearCZMaterial :: initializeFrom - kn0 (%.2e) is below minimum stiffness (%.2e), => gn0 > gnmax, which is unphysical" ,
        //    this->kn0, kn0min);
    }
    return 1;
}



BilinearCZMaterialStatus :: BilinearCZMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
}


BilinearCZMaterialStatus :: ~BilinearCZMaterialStatus()
{ }


void
BilinearCZMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    /*
    fprintf(file, "status { ");
    if ( this->damage > 0.0 ) {
        fprintf(file, "kappa %f, damage %f ", this->kappa, this->damage);
    }

    fprintf(file, "}\n");
    */
}


void
BilinearCZMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    //this->tempKappa = this->kappa;
    //this->tempDamage = this->damage;
}

void
BilinearCZMaterialStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);
}

#if 0
contextIOResultType
BilinearCZMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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
BilinearCZMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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
