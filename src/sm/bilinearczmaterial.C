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
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
namespace oofem {

REGISTER_Material( BilinearCZMaterial );

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
BilinearCZMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
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


    answer.resize( jumpVector.giveSize() -10);
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
        //printf("Softening branch...\n");
    } else {
        answer.at(1) = this->ks0 * gs1; 
        answer.at(2) = this->ks0 * gs2;
        answer.at(3) = 0.0; // failed
        //printf("Failed...\n");
    }

    // update gp
    status->letTempStrainVectorBe(jumpVector);
    status->letTempStressVectorBe(answer);
    
}


void
BilinearCZMaterial :: giveStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode rMode,
                                               GaussPoint *gp, 
                                               TimeStep *atTime)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dInterface:
    case _3dMat:
        give3dInterfaceMaterialStiffnessMatrix(answer, rMode, gp, atTime);
        break;
    default:
        //StructuralMaterial :: giveCharacteristicMatrix(answer, rMode, gp, atTime);
        StructuralMaterial ::give3dMaterialStiffnessMatrix(answer, rMode, gp, atTime);
    }
}


void
BilinearCZMaterial :: give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                                                     GaussPoint *gp, TimeStep *atTime)
{
    BilinearCZMaterialStatus *status = static_cast< BilinearCZMaterialStatus * >( this->giveStatus(gp) );

    FloatArray jumpVector;
    jumpVector = status->giveTempStrainVector();
    //@todo for now only study normal stress
    // jumpVector = [shear 1 shear 2 normal]
    double gn  = jumpVector.at(3);
    //double gs1 = jumpVector.at(1);
    //double gs2 = jumpVector.at(2);


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


int
BilinearCZMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    //BilinearCZMaterialStatus *status = static_cast< BilinearCZMaterialStatus * >( this->giveStatus(aGaussPoint) );
    if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.at(1) = 0.0; // no damage
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
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
    
    this->checkConsistency();                                // check validity of the material paramters
    this->printYourself();
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
     //   OOFEM_ERROR3("BilinearCZMaterial :: initializeFrom - kn0 (%.2e) is below minimum stiffness (%.2e), => gn0 > gnmax, which is unphysical" ,
     //       this->kn0, kn0min);
    }
    return 1;
}

void
BilinearCZMaterial :: printYourself()
{
    printf("Paramters for BilinearCZMaterial: \n");

    printf("-Strength paramters \n");
    printf("  sigfn = %e \n", this->sigfn);
    printf("  GIc   = %e \n", this->GIc);
    //printf("\n");

    printf("-Stiffness parameters \n");
    printf("  kn0   = %e \n", this->kn0);
    printf("  kn1   = %e \n", this->kn1);
    printf("  knc   = %e \n", this->knc);
    //printf("\n");

    printf("-jump limits \n");
    printf("  gn0   = %e \n", this->gn0);
    printf("  gnmax = %e \n", this->gnmax);
    
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
