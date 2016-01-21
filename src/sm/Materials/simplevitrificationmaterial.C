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

#include "simplevitrificationmaterial.h"
#include "gausspoint.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(SimpleVitrificationMaterial);

SimpleVitrificationMaterial :: ~SimpleVitrificationMaterial() { }


IRResultType SimpleVitrificationMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    IR_GIVE_FIELD(ir, this->vitrTime, _IFT_SimpleVitrificationMaterial_vitrificationTime);

    IR_GIVE_FIELD(ir, this->E, _IFT_SimpleVitrificationMaterial_E);
    IR_GIVE_FIELD(ir, this->nu, _IFT_SimpleVitrificationMaterial_nu);
    IR_GIVE_FIELD(ir, this->G, _IFT_SimpleVitrificationMaterial_G);
    IR_GIVE_FIELD(ir, this->alpha, _IFT_SimpleVitrificationMaterial_alpha);

    IR_GIVE_FIELD(ir, this->E_r, _IFT_SimpleVitrificationMaterial_E_r);
    IR_GIVE_FIELD(ir, this->nu_r, _IFT_SimpleVitrificationMaterial_nu_r);
    IR_GIVE_FIELD(ir, this->G_r, _IFT_SimpleVitrificationMaterial_G_r);
    IR_GIVE_FIELD(ir, this->alpha_r, _IFT_SimpleVitrificationMaterial_alpha_r);

    return StructuralMaterial :: initializeFrom(ir);
}


void SimpleVitrificationMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);

    input.setField(this->vitrTime, _IFT_SimpleVitrificationMaterial_vitrificationTime);

    input.setField(this->E, _IFT_SimpleVitrificationMaterial_E);
    input.setField(this->nu, _IFT_SimpleVitrificationMaterial_nu);
    input.setField(this->G, _IFT_SimpleVitrificationMaterial_G);
    input.setField(this->alpha, _IFT_SimpleVitrificationMaterial_alpha);

    input.setField(this->E_r, _IFT_SimpleVitrificationMaterial_E_r);
    input.setField(this->nu_r, _IFT_SimpleVitrificationMaterial_nu_r);
    input.setField(this->G_r, _IFT_SimpleVitrificationMaterial_G_r);
    input.setField(this->alpha_r, _IFT_SimpleVitrificationMaterial_alpha_r);
}


int SimpleVitrificationMaterial :: checkConsistency()
{
    return true;
}



void SimpleVitrificationMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                                  MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    double eksi, nxz, nyz, nxy, nzx, nzy, nyx;
    bool vitr = tStep->giveIntrinsicTime() < this->vitrTime;
    const FloatArray &activeNu = vitr ? this->nu_r : this->nu;
    const FloatArray &activeE = vitr ? this->E_r : this->E;
    const FloatArray &activeG = vitr ? this->G_r : this->G;

    nyz = nzy = activeNu.at(1);
    nxz = nzx = activeNu.at(2);
    nxy = nyx = activeNu.at(3);

    eksi = 1. - ( nxy * nyx + nyz * nzy + nzx * nxz ) - ( nxy * nyz * nzx + nyx * nzy * nxz );

    answer.resize(6, 6);
    answer.zero();
    // switched letters from original oofem -> now produces same material stiffness matrix as Abaqus method
    answer.at(1, 1) = activeE.at(1) * ( 1. - nyz * nzy ) / eksi;
    answer.at(1, 2) = activeE.at(2) * ( nxy + nxz * nzy ) / eksi;
    answer.at(1, 3) = activeE.at(3) * ( nxz + nyz * nxy ) / eksi;
    answer.at(2, 2) = activeE.at(2) * ( 1. - nxz * nzx ) / eksi;
    answer.at(2, 3) = activeE.at(3) * ( nyz + nyx * nxz ) / eksi;
    answer.at(3, 3) = activeE.at(3) * ( 1. - nyx * nxy ) / eksi;

    // define the lower triangle
    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < i; j++ ) {
            answer.at(i, j) = answer.at(j, i);
        }
    }

    answer.at(4, 4) = activeG.at(1);
    answer.at(5, 5) = activeG.at(2);
    answer.at(6, 6) = activeG.at(3);
}


void SimpleVitrificationMaterial :: giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                                            const FloatArray &reducedStrain, TimeStep *tStep)
{
    FloatArray strainVector;
    FloatMatrix d;
    FloatArray deltaStrain;

    StructuralMaterialStatus *status = dynamic_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    this->giveStressDependentPartOfStrainVector(strainVector, gp, reducedStrain, tStep, VM_Total);

    deltaStrain.beDifferenceOf( strainVector, status->giveStrainVector() );

    this->give3dMaterialStiffnessMatrix(d, TangentStiffness, gp, tStep);

    FloatArray deltaStress;
    deltaStress.beProductOf(d, deltaStrain);

    answer = status->giveStressVector();
    answer.add(deltaStress);

    // update gp
    status->letTempStrainVectorBe(reducedStrain);
    status->letTempStressVectorBe(answer);
}


void SimpleVitrificationMaterial :: giveThermalDilatationVector(FloatArray &answer,
                                                                GaussPoint *gp, TimeStep *tStep)
{
    bool vitr = tStep->giveIntrinsicTime() < this->vitrTime;
    answer.resize(6);
    answer.at(1) = vitr ? this->alpha_r.at(1) : this->alpha.at(1);
    answer.at(2) = vitr ? this->alpha_r.at(2) : this->alpha.at(2);
    answer.at(3) = vitr ? this->alpha_r.at(3) : this->alpha.at(3);
}


MaterialStatus *SimpleVitrificationMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(1, domain, gp);
}
} // end namespace oofem
