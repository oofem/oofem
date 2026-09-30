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

#include "mpm.h"
#include "material.h"
#include "matstatus.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "internalstatetype.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for StateDependentTangentMaterial
//@{
#define _IFT_StateDependentTangentMaterial_Name "statedeptangentmat"
#define _IFT_StateDependentTangentMaterial_E "e"
#define _IFT_StateDependentTangentMaterial_beta "beta"
//@}

namespace oofem {

/**
 * @brief Status of StateDependentTangentMaterial.
 *
 * Deliberately does not pre-size its vectors in initTempStatus: a quantity that was never pushed
 * stays empty, so a query arriving before any push is detectable rather than reading zeros.
 */
class StateDependentTangentMaterialStatus : public MaterialStatus
{
protected:
    FloatArray strainVector, tempStrainVector;
    FloatArray stressVector, tempStressVector;

public:
    StateDependentTangentMaterialStatus(GaussPoint *g) : MaterialStatus(g) {}

    const FloatArray &giveTempStrainVector() const { return tempStrainVector; }
    const FloatArray &giveTempStressVector() const { return tempStressVector; }
    void letTempStrainVectorBe(const FloatArray &v) { tempStrainVector = v; }
    void letTempStressVectorBe(const FloatArray &v) { tempStressVector = v; }

    void initTempStatus() override {
        MaterialStatus::initTempStatus();
        tempStrainVector = strainVector;
        tempStressVector = stressVector;
    }
    void updateYourself(TimeStep *tStep) override {
        MaterialStatus::updateYourself(tStep);
        strainVector = tempStrainVector;
        stressVector = tempStressVector;
    }
    const char *giveClassName() const override { return "StateDependentTangentMaterialStatus"; }
};


/**
 * @brief Test material whose tangent depends on the pushed state.
 *
 * Exists to make the push/pull contract observable. Every other material in the tree returns a
 * tangent built only from constants, so no in-tree deck can tell whether the tangent was evaluated
 * against current or stale state -- which is precisely why the ordering defect this interface fixes
 * went unnoticed until a python material with a real retention curve hit it.
 *
 * The constitutive law is deliberately trivial, a bar with a strain-softening/stiffening modulus:
 *
 *     sigma = E * eps ,     dsigma/deps = E * (1 + beta * eps)
 *
 * With beta = 0 it is ordinary linear elasticity.
 */
class StateDependentTangentMaterial : public Material
{
protected:
    double e = 1.0;
    double beta = 0.0;

public:
    StateDependentTangentMaterial(int n, Domain *d) : Material(n, d) {}

    StateVariableLayout giveStateVariableIDs(MaterialMode mmode) const override {
        return { { FT_Displacements, SO_SymmetricGradient } };
    }

    void updateTempState(const FloatArray &stateVector, GaussPoint *gp, TimeStep *tStep) override {
        auto status = static_cast< StateDependentTangentMaterialStatus * >( this->giveStatus(gp) );
        FloatArray sig = stateVector;
        sig.times(this->e);
        status->letTempStrainVectorBe(stateVector);
        status->letTempStressVectorBe(sig);
    }

    void giveCharacteristicMatrix(FloatMatrix &answer, MatResponseMode type, GaussPoint *gp, TimeStep *tStep) const override {
        if ( type != TangentStiffness && type != SecantStiffness ) {
            this->Material::giveCharacteristicMatrix(answer, type, gp, tStep);
            return;
        }
        auto status = static_cast< StateDependentTangentMaterialStatus * >( this->giveStatus(gp) );
        const FloatArray &eps = status->giveTempStrainVector();
        if ( eps.isEmpty() ) {
            OOFEM_ERROR("tangent queried on element %d GP %d before any updateTempState established the strain",
                        gp->giveElement()->giveNumber(), gp->giveNumber());
        }
        answer.resize(1, 1);
        answer.at(1, 1) = this->e * ( 1.0 + this->beta * eps.at(1) );
    }

    void giveCharacteristicVector(FloatArray &answer, MatResponseMode type, GaussPoint *gp, TimeStep *tStep) const override {
        if ( type != Stress ) {
            this->Material::giveCharacteristicVector(answer, type, gp, tStep);
            return;
        }
        auto status = static_cast< StateDependentTangentMaterialStatus * >( this->giveStatus(gp) );
        if ( status->giveTempStressVector().isEmpty() ) {
            OOFEM_ERROR("stress queried on element %d GP %d before any updateTempState established it",
                        gp->giveElement()->giveNumber(), gp->giveNumber());
        }
        answer = status->giveTempStressVector();
    }

    bool hasMaterialModeCapability(MaterialMode mode) const override { return mode == _1dMat; }

    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override {
        Material::initializeFrom(ir);
        IR_GIVE_OPTIONAL_FIELD(ir, e, _IFT_StateDependentTangentMaterial_E);
        IR_GIVE_OPTIONAL_FIELD(ir, beta, _IFT_StateDependentTangentMaterial_beta);
    }
    void giveInputRecord(DynamicInputRecord &input) override {}

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override {
        return std::make_unique<StateDependentTangentMaterialStatus>(gp);
    }

    const char *giveClassName() const override { return "StateDependentTangentMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_StateDependentTangentMaterial_Name; }
};

REGISTER_Material(StateDependentTangentMaterial);

} // end namespace oofem
