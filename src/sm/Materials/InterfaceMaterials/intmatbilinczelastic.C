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

#include "intmatbilinczelastic.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(IntMatBilinearCZElastic);

IntMatBilinearCZElastic :: IntMatBilinearCZElastic(int n, Domain *d) : StructuralInterfaceMaterial(n, d) { }


IntMatBilinearCZElastic :: ~IntMatBilinearCZElastic() { }


void
IntMatBilinearCZElastic :: giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jumpVector,
                                                  const FloatMatrix &F, TimeStep *tStep)
{
    IntMatBilinearCZElasticStatus *status = static_cast< IntMatBilinearCZElasticStatus * >( this->giveStatus(gp) );

    answer.resize( jumpVector.giveSize() );
    answer.zero();
    /// @todo only study normal stress for now 
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
        answer.at(3) = this->sigfn + this->kn1 * ( gn - this->gn0 );
        //printf("Softening branch...\n");
    } else {
        answer.at(1) = this->ks0 * gs1;
        answer.at(2) = this->ks0 * gs2;
        answer.at(3) = 0.0; // failed
        //printf("Failed...\n");
    }

    // update gp
    status->letTempJumpBe(jumpVector);
    status->letTempFirstPKTractionBe(answer);
}


void
IntMatBilinearCZElastic :: give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    IntMatBilinearCZElasticStatus *status = static_cast< IntMatBilinearCZElasticStatus * >( this->giveStatus(gp) );

    const FloatArray &jumpVector = status->giveTempJump();
    double gn = jumpVector.at(3);

    answer.resize(3, 3);
    answer.zero();
    if ( gn <= 0.0 ) { // compresssion
        answer.at(1, 1) = this->ks0;
        answer.at(2, 2) = this->ks0;
        answer.at(3, 3) = this->knc;
    } else if ( gn <= this->gn0 ) { // first linear part
        answer.at(1, 1) = this->ks0;
        answer.at(2, 2) = this->ks0;
        answer.at(3, 3) = this->kn0;
        //printf("Linear branch...\n");
    } else if ( gn <= this->gnmax  ) { // softening branch
        answer.at(1, 1) = this->ks0; // no degradation in shear
        answer.at(2, 2) = this->ks0;
        answer.at(3, 3) = this->kn1;
        //printf("Softening branch...\n");
    } else {
        answer.at(1, 1) = this->ks0;
        answer.at(2, 2) = this->ks0;
        answer.at(3, 3) = 0.0; // failed
        //printf("Failed...\n");
    }
}


int
IntMatBilinearCZElastic :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.at(1) = 0.0; // no damage
        return 1;
    } else {
        return StructuralInterfaceMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}

const double tolerance = 1.0e-12; // small number
IRResultType
IntMatBilinearCZElastic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, kn0, _IFT_IntMatBilinearCZElastic_kn);
    this->knc = kn0;                        // Defaults to the same stiffness in compression and tension
    IR_GIVE_OPTIONAL_FIELD(ir, this->knc, _IFT_IntMatBilinearCZElastic_knc);

    this->ks0 = 0.0;                        // Defaults to no shear stiffness
    IR_GIVE_OPTIONAL_FIELD(ir, ks0, _IFT_IntMatBilinearCZElastic_ks);

    IR_GIVE_FIELD(ir, GIc, _IFT_IntMatBilinearCZElastic_g1c);

    this->sigfn = 1.0e50; // defaults to infinite strength @todo but then it will not be bilinear only linear
    this->sigfs = 1.0e50;
    IR_GIVE_OPTIONAL_FIELD(ir, sigfn, _IFT_IntMatBilinearCZElastic_sigfn);

    this->gn0 = sigfn / ( kn0 + tolerance );                   // normal jump at damage initiation
    this->gs0 = sigfs / ( ks0 + tolerance );                   // shear jump at damage initiation
    this->gnmax = 2.0 * GIc / sigfn;                         // @todo defaults to zero - will this cause problems?
    this->kn1 = -this->sigfn / ( this->gnmax - this->gn0 );  // slope during softening part in normal dir

    result = StructuralInterfaceMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    // check validity of the material paramters
    double kn0min = 0.5 * this->sigfn * this->sigfn / this->GIc;
    if ( this->kn0 < 0.0 ) {
        OOFEM_WARNING("stiffness kn0 is negative (%.2e)", this->kn0);
        return IRRT_BAD_FORMAT;
    } else if ( this->ks0 < 0.0 ) {
        OOFEM_WARNING("stiffness ks0 is negative (%.2e)", this->ks0);
        return IRRT_BAD_FORMAT;
    } else if ( this->GIc < 0.0 ) {
        OOFEM_WARNING("GIc is negative (%.2e)", this->GIc);
        return IRRT_BAD_FORMAT;
    } else if ( this->kn0 < kn0min  ) { // => gn0 > gnmax
        //OOFEM_WARNING("kn0 (%.2e) is below minimum stiffness (%.2e), => gn0 > gnmax, which is unphysical", this->kn0, kn0min);
        //return IRRT_BAD_FORMAT;
    }
    return IRRT_OK;
}

void IntMatBilinearCZElastic :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);

    input.setField(kn0, _IFT_IntMatBilinearCZElastic_kn);
    input.setField(knc, _IFT_IntMatBilinearCZElastic_knc);
    input.setField(ks0, _IFT_IntMatBilinearCZElastic_ks);
    input.setField(GIc, _IFT_IntMatBilinearCZElastic_g1c);
    input.setField(sigfn, _IFT_IntMatBilinearCZElastic_sigfn);
}

int
IntMatBilinearCZElastic :: checkConsistency()
{
    return 1;
}

void
IntMatBilinearCZElastic :: printYourself()
{
    printf("Paramters for IntMatBilinearCZElastic: \n");

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

} // end namespace oofem
