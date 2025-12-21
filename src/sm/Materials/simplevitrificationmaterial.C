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

#include "simplevitrificationmaterial.h"
#include "gausspoint.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(SimpleVitrificationMaterial);


void SimpleVitrificationMaterial :: initializeFrom(InputRecord &ir)
{
    StructuralMaterial :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, this->vitrTime, _IFT_SimpleVitrificationMaterial_vitrificationTime);

    IR_GIVE_FIELD(ir, this->E, _IFT_SimpleVitrificationMaterial_E);
    IR_GIVE_FIELD(ir, this->nu, _IFT_SimpleVitrificationMaterial_nu);
    IR_GIVE_FIELD(ir, this->G, _IFT_SimpleVitrificationMaterial_G);
    IR_GIVE_FIELD(ir, this->alpha, _IFT_SimpleVitrificationMaterial_alpha);

    IR_GIVE_FIELD(ir, this->E_r, _IFT_SimpleVitrificationMaterial_E_r);
    IR_GIVE_FIELD(ir, this->nu_r, _IFT_SimpleVitrificationMaterial_nu_r);
    IR_GIVE_FIELD(ir, this->G_r, _IFT_SimpleVitrificationMaterial_G_r);
    IR_GIVE_FIELD(ir, this->alpha_r, _IFT_SimpleVitrificationMaterial_alpha_r);
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


FloatMatrixF<6,6>
SimpleVitrificationMaterial :: computeTangent(bool vitr) const
{
    const auto &activeNu = vitr ? this->nu_r : this->nu;
    const auto &activeE = vitr ? this->E_r : this->E;
    const auto &activeG = vitr ? this->G_r : this->G;

    //double [nyz, nxz, nxy] = activeNu; // c++17
    //double [nzy, nzx, nyx] = activeNu;

    double nyz = activeNu.at(1);
    double nxz = activeNu.at(2);
    double nxy = activeNu.at(3);

    double nzy = activeNu.at(1);
    double nzx = activeNu.at(2);
    double nyx = activeNu.at(3);

    double eksi = 1. - ( nxy * nyx + nyz * nzy + nzx * nxz ) - ( nxy * nyz * nzx + nyx * nzy * nxz );

    FloatMatrixF<6,6> tangent;
    // switched letters from original oofem -> now produces same material stiffness matrix as Abaqus method
    tangent.at(1, 1) = activeE.at(1) * ( 1. - nyz * nzy ) / eksi;
    tangent.at(1, 2) = activeE.at(2) * ( nxy + nxz * nzy ) / eksi;
    tangent.at(1, 3) = activeE.at(3) * ( nxz + nyz * nxy ) / eksi;
    tangent.at(2, 2) = activeE.at(2) * ( 1. - nxz * nzx ) / eksi;
    tangent.at(2, 3) = activeE.at(3) * ( nyz + nyx * nxz ) / eksi;
    tangent.at(3, 3) = activeE.at(3) * ( 1. - nyx * nxy ) / eksi;

    // define the lower triangle
    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < i; j++ ) {
            tangent.at(i, j) = tangent.at(j, i);
        }
    }

    tangent.at(4, 4) = activeG.at(1);
    tangent.at(5, 5) = activeG.at(2);
    tangent.at(6, 6) = activeG.at(3);
    
    return tangent;
}



FloatMatrixF<6,6>
SimpleVitrificationMaterial :: give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    bool vitr = tStep->giveIntrinsicTime() < this->vitrTime;
    return this->computeTangent(vitr);
}


FloatArrayF<6>
SimpleVitrificationMaterial :: giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = dynamic_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    
    auto thermalStrain = computeStressIndependentStrainVector_3d(gp, tStep, VM_Total);
    auto deltaStrain = strain - thermalStrain - FloatArrayF<6>(status->giveStrainVector());

    bool vitr = tStep->giveIntrinsicTime() < this->vitrTime;
    auto d = this->computeTangent(vitr);

    auto deltaStress = dot(d, deltaStrain);

    auto stress = status->giveStressVector() + deltaStress;

    // update gp
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(stress);
    
    return stress;
}


FloatArrayF<6> SimpleVitrificationMaterial :: giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const
{
    bool vitr = tStep->giveIntrinsicTime() < this->vitrTime;
    return {
        vitr ? this->alpha_r.at(1) : this->alpha.at(1),
        vitr ? this->alpha_r.at(2) : this->alpha.at(2),
        vitr ? this->alpha_r.at(3) : this->alpha.at(3),
        0., 0., 0.,
    };
}


std::unique_ptr<MaterialStatus> SimpleVitrificationMaterial :: CreateStatus(GaussPoint *gp) const
{
    return std::make_unique<StructuralMaterialStatus>(gp);
}
} // end namespace oofem
