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

#include "simpleinterfacemat.h"
#include "sm/Elements/structuralelement.h"
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


FloatArrayF<3>
SimpleInterfaceMaterial :: giveEngTraction_3d(const FloatArrayF<3> &jump, GaussPoint *gp, TimeStep *tStep) const
{
    SimpleInterfaceMaterialStatus *status = static_cast< SimpleInterfaceMaterialStatus * >( this->giveStatus(gp) );
    bool shearYieldingFlag = false;
    double normalStrain = jump.at(1);
    double normalStress, maxShearStress;
    double shift = -this->kn * this->stiffCoeff * normalClearance;

    if (regularizedModel) {
        normalStress = 0.5*this->kn*(normalClearance+normalStrain)-
            ((this->kn*0.5)/(m))*log(fabs(cosh(m*(normalClearance+normalStrain))))+
            0.5*this->stiffCoeff*this->kn*(normalClearance+normalStrain)+
            ((this->stiffCoeff*this->kn*0.5)/(m))*log(fabs(cosh(m*(normalClearance+normalStrain))));
            maxShearStress = fabs(normalStress) * this->frictCoeff;
    } else {
      if ( normalStrain + normalClearance <= 0. ) {
        normalStress = this->kn * ( normalStrain + normalClearance ) + shift; //in compression and after the clearance gap closed
        maxShearStress = fabs(normalStress) * this->frictCoeff;
      } else {
        normalStress = this->kn * this->stiffCoeff * ( normalStrain + normalClearance ) + shift;
        maxShearStress = 0.;
      }
    }

    FloatArrayF<2> tempShearStressShift = status->giveShearStressShift();
    FloatArrayF<2> shearJump = {jump.at(2), jump.at(3)};
    auto shearTraction = this->ks * shearJump - tempShearStressShift;
    double dpn = norm(shearTraction);
    if ( dpn > maxShearStress ) {
        shearYieldingFlag = true;
        shearTraction *= maxShearStress / dpn;
    }

    tempShearStressShift = this->ks * shearJump - shearTraction;

    double lim = 1.e+50;
    FloatArrayF<3> answer = {max(min(normalStress, lim), -lim), shearTraction.at(1), shearTraction.at(2)};

    // update gp
    status->setShearYieldingFlag(shearYieldingFlag);
    status->setTempShearStressShift(tempShearStressShift);
    status->letTempJumpBe(jump);
    status->letTempTractionBe(answer);
    return answer;
}


FloatMatrixF<3,3>
SimpleInterfaceMaterial :: give3dStiffnessMatrix_Eng(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    SimpleInterfaceMaterialStatus *status = static_cast< SimpleInterfaceMaterialStatus * >( this->giveStatus(gp) );
    double normalJump = status->giveTempJump().at(1);

    FloatMatrixF<3,3> answer;

    if (regularizedModel) {
      answer.at(1, 1) = this->kn+0.5*tanh(m*(normalJump+this->normalClearance))*(this->kn*this->stiffCoeff-this->kn);
      if ( status->giveShearYieldingFlag() ) {
        FloatArray jump = status->giveTempJump();
        FloatArray traction = status->giveTempTraction();
        answer.at(2, 2) = answer.at(3, 3) = traction.at(2)/jump.at(2);
      } else {
        answer.at(2, 2) = answer.at(3, 3) = ks;//this->kn; //in compression and after the clearance gap closed
      }
    } else {

      if ( rMode == SecantStiffness || rMode == TangentStiffness ) {
        if ( normalJump + normalClearance <= 0. ) {
          answer.at(1, 1) = kn;
          if ( status->giveShearYieldingFlag() ) {
            FloatArray jump = status->giveTempJump();
            FloatArray traction = status->giveTempTraction();
            
            answer.at(2, 2) = answer.at(3, 3) = traction.at(2)/jump.at(2);
          } else {
            answer.at(2, 2) = answer.at(3, 3) = ks;//this->kn; //in compression and after the clearance gap closed
          }
        } else {
          answer.at(1, 1) = this->kn * this->stiffCoeff;
          answer.at(2, 2) = answer.at(3, 3) = 0;
        }
      } else {
        answer.at(1, 1) = kn;
        answer.at(2, 2) = answer.at(3, 3) = this->ks;
      }
    }
    return answer;
}


int
SimpleInterfaceMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    return StructuralInterfaceMaterial :: giveIPValue(answer, gp, type, tStep);
}


void
SimpleInterfaceMaterial :: initializeFrom(InputRecord &ir)
{
    StructuralInterfaceMaterial :: initializeFrom(ir);

    frictCoeff = 0.;
    stiffCoeff = 0.;
    normalClearance = 0.;
    IR_GIVE_FIELD(ir, kn, _IFT_SimpleInterfaceMaterial_kn);
    ks = kn;
    IR_GIVE_OPTIONAL_FIELD(ir, ks, _IFT_SimpleInterfaceMaterial_ks);
    IR_GIVE_OPTIONAL_FIELD(ir, frictCoeff, _IFT_SimpleInterfaceMaterial_frictCoeff);
    IR_GIVE_OPTIONAL_FIELD(ir, stiffCoeff, _IFT_SimpleInterfaceMaterial_stiffCoeff);
    IR_GIVE_OPTIONAL_FIELD(ir, normalClearance, _IFT_SimpleInterfaceMaterial_normalClearance);
    IR_GIVE_OPTIONAL_FIELD(ir, regularizedModel, _IFT_SimpleInterfaceMaterial_regularizedModel);
    IR_GIVE_OPTIONAL_FIELD(ir, m, _IFT_SimpleInterfaceMaterial_regularizationCoeff);

}


void
SimpleInterfaceMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);
    input.setField(this->kn, _IFT_SimpleInterfaceMaterial_kn);
    input.setField(this->frictCoeff, _IFT_SimpleInterfaceMaterial_frictCoeff);
    input.setField(this->stiffCoeff, _IFT_SimpleInterfaceMaterial_stiffCoeff);
    input.setField(this->normalClearance, _IFT_SimpleInterfaceMaterial_normalClearance);

    input.setField(this->regularizedModel, _IFT_SimpleInterfaceMaterial_regularizedModel);
    input.setField(this->m, _IFT_SimpleInterfaceMaterial_regularizationCoeff);
    
}


SimpleInterfaceMaterialStatus :: SimpleInterfaceMaterialStatus(GaussPoint *g) : StructuralInterfaceMaterialStatus(g)
{}


void
SimpleInterfaceMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
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


void
SimpleInterfaceMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralInterfaceMaterialStatus :: saveContext(stream, mode);

    //if ( !stream.write(kappa) ) {
    //THROW_CIOERR(CIO_IOERR);
    //}
}


void
SimpleInterfaceMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralInterfaceMaterialStatus :: restoreContext(stream, mode);

    //if ( !stream.read(kappa) ) {
    //THROW_CIOERR(CIO_IOERR);
    //}
}
} // end namespace oofem
