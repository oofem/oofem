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


FloatArrayF<3>
IntMatBilinearCZElastic :: giveFirstPKTraction_3d(const FloatArrayF<3> &jump, const FloatMatrixF<3,3> &F, GaussPoint *gp, TimeStep *tStep) const 
{
    IntMatBilinearCZElasticStatus *status = static_cast< IntMatBilinearCZElasticStatus * >( this->giveStatus(gp) );

    /// @todo only study normal stress for now 
    // no degradation in shear
    // jumpVector = [normal shear1 shear2]
    //auto [gn, gs1, gs2] = jump;
    double gn  = jump.at(1);
    double gs1 = jump.at(2);
    double gs2 = jump.at(3);

    double normal;
    if ( gn <= 0.0 ) { // compresssion
        normal = this->knc * gn;
    } else if ( gn <= this->gn0 ) { // first linear part
        normal = this->kn0 * gn;
    } else if ( gn <= this->gnmax  ) {  // softening branch
        normal = this->sigfn + this->kn1 * ( gn - this->gn0 );
        //printf("Softening branch...\n");
    } else {
        normal = 0.0; // failed
        //printf("Failed...\n");
    }
    FloatArrayF<3> answer = {normal, this->ks0 * gs1, this->ks0 * gs2};

    // update gp
    status->letTempJumpBe(jump);
    status->letTempFirstPKTractionBe(answer);
    status->letTempTractionBe(answer);
    return answer;
}


FloatMatrixF<3,3>
IntMatBilinearCZElastic :: give3dStiffnessMatrix_dTdj(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    IntMatBilinearCZElasticStatus *status = static_cast< IntMatBilinearCZElasticStatus * >( this->giveStatus(gp) );

    const auto &jumpVector = status->giveTempJump();
    double gn = jumpVector.at(1);

    double normal;
    if ( gn <= 0.0 ) { // compresssion
        normal = this->knc;
    } else if ( gn <= this->gn0 ) { // first linear part
        normal = this->kn0;
        //printf("Linear branch...\n");
    } else if ( gn <= this->gnmax  ) { // softening branch
        normal = this->kn1;
        //printf("Softening branch...\n");
    } else {
        normal = 0.0; // failed
        //printf("Failed...\n");
    }
    return diag<3>({normal, this->ks0, this->ks0});
}


const double tolerance = 1.0e-12; // small number
void
IntMatBilinearCZElastic :: initializeFrom(InputRecord &ir)
{
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

    StructuralInterfaceMaterial :: initializeFrom(ir);

    // check validity of the material paramters
    double kn0min = 0.5 * this->sigfn * this->sigfn / this->GIc;
    if ( this->kn0 < 0.0 ) {
        throw ValueInputException(ir, _IFT_IntMatBilinearCZElastic_kn, "stiffness kn0 is negative");
    } else if ( this->ks0 < 0.0 ) {
        throw ValueInputException(ir, _IFT_IntMatBilinearCZElastic_ks, "stiffness ks0 is negative");
    } else if ( this->GIc < 0.0 ) {
        throw ValueInputException(ir, _IFT_IntMatBilinearCZElastic_g1c, "GIc is negative");
    } else if ( this->kn0 < kn0min  ) { // => gn0 > gnmax
        //throw ValueInputException(ir, _IFT_IntMatBilinearCZElastic_kn, "kn0 is below minimum stiffness, which is unphysical");
    }
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
