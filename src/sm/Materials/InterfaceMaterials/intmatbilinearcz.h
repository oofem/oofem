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

#ifndef INTMATBILINEARCZ_H_
#define INTMATBILINEARCZ_H_

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for IntMatBilinearCZ
//@{
#define _IFT_IntMatBilinearCZ_Name "intmatbilinearcz"
#define _IFT_IntMatBilinearCZ_PenaltyStiffness "kn"
#define _IFT_IntMatBilinearCZ_g1c "g1c"
#define _IFT_IntMatBilinearCZ_g2c "g2c"
#define _IFT_IntMatBilinearCZ_mu "mu"
#define _IFT_IntMatBilinearCZ_gamma "gamma"
#define _IFT_IntMatBilinearCZ_sigf "sigf"
#define _IFT_IntMatBilinearCZ_semiexplicit "semiexplicit"
//@}

namespace oofem {
/**
 * This class implements associated Material Status for IntMatBilinearCZFagerstrom
 */
class IntMatBilinearCZStatus : public StructuralInterfaceMaterialStatus
{
public:
    IntMatBilinearCZStatus(GaussPoint * g);

    /// damage variable
    double mDamageNew = 0., mDamageOld = 0.;

    /// Traction
    FloatArrayF<3> mTractionOld, mTractionNew;

    /// Discontinuity
    FloatArrayF<3> mJumpOld, mJumpNew;

    /**
     * Increment of plastic multiplier. Storing this allows
     * semi-explicit update of damage.
     */
    double mPlastMultIncNew = 0., mPlastMultIncOld = 0.;

    const char *giveClassName() const override { return "IntMatBilinearCZStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    /// Functions for MaterialStatusMapperInterface
    void copyStateVariables(const MaterialStatus &iStatus) override;
    void addStateVariables(const MaterialStatus &iStatus) override;
};


/**
 * Bilinear cohesive zone model.
 * Created on: Oct 20, 2013
 * @author Erik Svenning
 */
class IntMatBilinearCZ : public StructuralInterfaceMaterial
{
protected:
    double mPenaltyStiffness = 0.;
    double mGIc = 0.;   // fracture energy, mode 1
    double mGIIc = 0.;  // fracture energy, mode 2
    double mSigmaF = 0.;  // max stress

    double mMu = 0.;    // loading function parameter
    double mGamma = 0.; // loading function parameter

    bool mSemiExplicit = false; // If semi-explicit time integration should be used

    int checkConsistency() override;

public:
    IntMatBilinearCZ(int n, Domain * d);

    const char *giveClassName() const override { return "IntMatBilinearCZ"; }
    const char *giveInputRecordName() const override { return _IFT_IntMatBilinearCZ_Name; }

    FloatArrayF<3> giveFirstPKTraction_3d(const FloatArrayF<3> &jump, const FloatMatrixF<3,3> &F, GaussPoint *gp, TimeStep *tStep) const override;

    // Dummy implementation, we must rely on numerical computation of the tangent.
    FloatMatrixF<3,3> give3dStiffnessMatrix_dTdj(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    bool hasAnalyticalTangentStiffness() const override { return false; }

private:
    // Help functions
    double computeYieldFunction(double iTractionNormal, double iTractionTang) const;
    FloatArrayF<3> computeTraction(const FloatArrayF<3> &iTTrial, double iPlastMultInc) const;

public:
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<IntMatBilinearCZStatus>(gp); }
    void printYourself() override;
};
} /* namespace oofem */
#endif /* INTMATBILINEARCZ_H_ */
