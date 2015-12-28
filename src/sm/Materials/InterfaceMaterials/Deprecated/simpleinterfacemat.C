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

#include "simpleinterfacemat.h"
#include "../sm/Elements/structuralelement.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(SimpleInterfaceMaterial);

SimpleInterfaceMaterial :: SimpleInterfaceMaterial(int n, Domain *d) : StructuralInterfaceMaterial(n, d) { }


SimpleInterfaceMaterial :: ~SimpleInterfaceMaterial() { }


void
SimpleInterfaceMaterial :: giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
{
    SimpleInterfaceMaterialStatus *status = static_cast< SimpleInterfaceMaterialStatus * >( this->giveStatus(gp) );
    bool shearYieldingFlag = false;
    FloatArray shearJump(2), shearTraction;
    FloatArray tempShearStressShift = status->giveShearStressShift();
    double normalStrain = jump.at(1);
    double normalStress, maxShearStress, dp;
    double shift = -this->kn * this->stiffCoeff * normalClearance;

    if ( normalStrain + normalClearance <= 0. ) {
        normalStress = this->kn * ( normalStrain + normalClearance ) + shift; //in compression and after the clearance gap closed
        maxShearStress = fabs(normalStress) * this->frictCoeff;
    } else {
        normalStress = this->kn * this->stiffCoeff * ( normalStrain + normalClearance ) + shift;
        maxShearStress = 0.;
    }


    shearJump.at(1) = jump.at(2);
    shearJump.at(2) = jump.at(3);
    shearTraction.beScaled(this->ks, shearJump);
    shearTraction.subtract(tempShearStressShift);
    dp = shearTraction.dotProduct(shearTraction, 2);
    if ( dp > maxShearStress * maxShearStress ) {
        shearYieldingFlag = true;
        shearTraction.times( maxShearStress / sqrt(dp) );
    }

    tempShearStressShift.beScaled(this->ks, shearJump);
    tempShearStressShift.subtract(shearTraction);

    double lim = 1.e+50;
    answer.resize(3);
    answer.at(1) = max(min(normalStress, lim), -lim);  //threshold on maximum/minimum
    answer.at(2) = shearTraction.at(1);
    answer.at(3) = shearTraction.at(2);

    // update gp
    status->setShearYieldingFlag(shearYieldingFlag);
    status->setTempShearStressShift(tempShearStressShift);
    status->letTempJumpBe(jump);
    status->letTempTractionBe(answer);
}


void
SimpleInterfaceMaterial :: give3dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    SimpleInterfaceMaterialStatus *status = static_cast< SimpleInterfaceMaterialStatus * >( this->giveStatus(gp) );
    double normalJump = status->giveTempJump().at(1);

    answer.resize(3, 3);
    if ( rMode == SecantStiffness || rMode == TangentStiffness ) {
        if ( normalJump + normalClearance <= 0. ) {
	      answer.at(1, 1) = kn;
	      if ( status->giveShearYieldingFlag() )
		answer.at(2, 2) = answer.at(3, 3) = 0;
	      else
		answer.at(2, 2) = answer.at(3, 3) = ks;//this->kn; //in compression and after the clearance gap closed
        } else {
            answer.at(1, 1) = this->kn * this->stiffCoeff;
            answer.at(2, 2) = answer.at(3, 3) = 0;
        }
    } else {
        answer.at(1, 1) = kn;
        answer.at(2, 2) = answer.at(3, 3) = this->ks;
    }
}


int
SimpleInterfaceMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    return StructuralInterfaceMaterial :: giveIPValue(answer, gp, type, tStep);
}


IRResultType
SimpleInterfaceMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    frictCoeff = 0.;
    stiffCoeff = 0.;
    normalClearance = 0.;
    IR_GIVE_FIELD(ir, kn, _IFT_SimpleInterfaceMaterial_kn);
    ks = kn;
    IR_GIVE_OPTIONAL_FIELD(ir, ks, _IFT_SimpleInterfaceMaterial_ks);
    IR_GIVE_OPTIONAL_FIELD(ir, frictCoeff, _IFT_SimpleInterfaceMaterial_frictCoeff);
    IR_GIVE_OPTIONAL_FIELD(ir, stiffCoeff, _IFT_SimpleInterfaceMaterial_stiffCoeff);
    IR_GIVE_OPTIONAL_FIELD(ir, normalClearance, _IFT_SimpleInterfaceMaterial_normalClearance);

    return StructuralInterfaceMaterial :: initializeFrom(ir);
}


void
SimpleInterfaceMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);
    input.setField(this->kn, _IFT_SimpleInterfaceMaterial_kn);
    input.setField(this->frictCoeff, _IFT_SimpleInterfaceMaterial_frictCoeff);
    input.setField(this->stiffCoeff, _IFT_SimpleInterfaceMaterial_stiffCoeff);
    input.setField(this->normalClearance, _IFT_SimpleInterfaceMaterial_normalClearance);
}


SimpleInterfaceMaterialStatus :: SimpleInterfaceMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralInterfaceMaterialStatus(n, d, g)
{
    shearStressShift.resize(2);
    tempShearStressShift.resize(2);
    shearStressShift.zero();
    tempShearStressShift.zero();
    shearYieldingFlag = false;
}


SimpleInterfaceMaterialStatus :: ~SimpleInterfaceMaterialStatus()
{ }


void
SimpleInterfaceMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf( file, "shearStressShift (%f, %f)", this->shearStressShift.at(1), this->shearStressShift.at(2) );
    fprintf(file, "}\n");
}


void
SimpleInterfaceMaterialStatus :: initTempStatus()
{
    StructuralInterfaceMaterialStatus :: initTempStatus();
    tempShearStressShift = shearStressShift;
}


void
SimpleInterfaceMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus :: updateYourself(tStep);
    shearStressShift = tempShearStressShift;
}


const FloatArray &
SimpleInterfaceMaterialStatus :: giveShearStressShift()
{
    return shearStressShift;
}


contextIOResultType
SimpleInterfaceMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralInterfaceMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    //if ( !stream.write(kappa) ) {
    //THROW_CIOERR(CIO_IOERR);
    //}

    return CIO_OK;
}


contextIOResultType
SimpleInterfaceMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralInterfaceMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    //if ( !stream.read(kappa) ) {
    //THROW_CIOERR(CIO_IOERR);
    //}

    return CIO_OK;
}
} // end namespace oofem
