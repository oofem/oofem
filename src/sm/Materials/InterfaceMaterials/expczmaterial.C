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
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
///@todo This unit has a *unit*! There are no small numbers, and we should change the code to avoid the need for this at all cost. / Mikael
const double tolerance = 1.0e-12; // small number

ExpCZMaterial :: ExpCZMaterial(int n, Domain *d) : StructuralInterfaceMaterial(n, d)
{ }


ExpCZMaterial :: ~ExpCZMaterial()
{ }


void
ExpCZMaterial :: giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
{
    ExpCZMaterialStatus *status = static_cast< ExpCZMaterialStatus * >( this->giveStatus(gp) );

    /// @todo only study normal stress for now
    // no degradation in shear
    // jumpVector = [shear 1 shear 2 normal]
    double gn  = jump.at(3);
    double gs1 = jump.at(1);
    double gs2 = jump.at(2);
    double gs  = sqrt(gs1 * gs1 + gs2 * gs2);

    double xin  = gn / ( gn0 + tolerance );
    double xit  = gs / ( gs0 + tolerance );
    double tn   = GIc / gn0 *exp(-xin) * ( xin * exp(-xit * xit) + ( 1.0 - q ) / ( r - 1.0 ) * ( 1.0 - exp(-xit * xit) ) * ( r - xin ) );

    answer.resize(3);
    answer.at(1) = 0.0;
    answer.at(2) = 0.0;
    answer.at(3) = tn;

    // update gp
    status->letTempJumpBe(jump);
    status->letTempTractionBe(answer);
}


void
ExpCZMaterial :: give3dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    ExpCZMaterialStatus *status = static_cast< ExpCZMaterialStatus * >( this->giveStatus(gp) );

    const FloatArray &jumpVector = status->giveTempJump();
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


int
ExpCZMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    ExpCZMaterialStatus *status = static_cast< ExpCZMaterialStatus * >( this->giveStatus(gp) );
    return StructuralInterfaceMaterial :: giveIPValue(answer, gp, type, tStep);
}


IRResultType
ExpCZMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->GIc, _IFT_ExpCZMaterial_g1c);

    this->sigfs = 0.0;
    IR_GIVE_FIELD(ir, this->sigfn, _IFT_ExpCZMaterial_sigfn);

    GIIc = 0.0;

    this->gn0 = GIc  / ( this->sigfn * M_E + tolerance );                // normal jump at max stress
    this->gs0 = GIIc / ( sqrt( 0.5 * M_E ) * this->sigfs + tolerance );      // shear jump at max stress


    this->q = GIIc / GIc;
    this->r = 0.0; // fix

    // check validity of the material paramters
    if ( this->GIc < 0.0 ) {
        OOFEM_WARNING("GIc is negative (%.2e)", this->GIc);
        return IRRT_BAD_FORMAT;
    }
    return IRRT_OK;
}

int
ExpCZMaterial :: checkConsistency()
{
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

ExpCZMaterialStatus :: ExpCZMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralInterfaceMaterialStatus(n, d, g)
{ }


ExpCZMaterialStatus :: ~ExpCZMaterialStatus()
{ }


void
ExpCZMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus :: printOutputAt(file, tStep);
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
    StructuralInterfaceMaterialStatus :: initTempStatus();
    //this->tempKappa = this->kappa;
    //this->tempDamage = this->damage;
}

void
ExpCZMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus :: updateYourself(tStep);
}

#if 0
contextIOResultType
ExpCZMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
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
ExpCZMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
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
#endif
} // end namespace oofem
