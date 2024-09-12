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

#ifndef RankineMatGrad_h

#include "rankinemat.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "graddamagematerialextensioninterface.h"
#include "cltypes.h"

///@name Input fields for RankineMatGrad
//@{
#define _IFT_RankineMatGrad_Name "rankmatgrad"
#define _IFT_RankineMatGrad_L "l"
#define _IFT_RankineMatGrad_m "m"
#define _IFT_RankineMatGrad_negligibleDamage "negligible_damage"
#define _IFT_RankineMatGrad_formulationType "formtype"
//@}

namespace oofem {
/**
 * Gradient rankine material status.
 */
class RankineMatGradStatus : public RankineMatStatus, public GradientDamageMaterialStatusExtensionInterface
{
protected:

    /**  Type characterizing the dependence of the internal lenght on variable of the state
     *  Note that the assigned numbers to enum values have to correspond to values
     *  used in initializeFrom to resolve internalLenghtDependence. If not, the consistency
     *  between initializeFrom and giveInputRecord methods is lost.
     */

    double kappa_nl = 0.;
    double kappa_hat = 0.;

public:

    RankineMatGradStatus(GaussPoint *g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    // definition
    const char *giveClassName() const override { return "RankineMatGradStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void setKappa_nl(double kap) { kappa_nl = kap; }
    void setKappa_hat(double kap) { kappa_hat = kap; }
    double giveKappa_nl() { return kappa_nl; }
    double giveKappa_hat() { return kappa_hat; }
    virtual double giveNonlocalCumulatedStrain() { return nonlocalDamageDrivingVariable; }
    virtual void setNonlocalCumulatedStrain(double nonlocalCumulatedStrain) { this->nonlocalDamageDrivingVariable = nonlocalCumulatedStrain; }
};


/**
 * Gradient Rankine material.
 */
class RankineMatGrad : public RankineMat, GradientDamageMaterialExtensionInterface
{
protected:
    double L = 0.;
    double mParam = 0.;
    double negligible_damage = 0.;

    enum GradientDamageFormulationType {
        GDFT_Standard = 0,
        GDFT_Eikonal = 2
    };

    GradientDamageFormulationType gradientDamageFormulationType = GDFT_Standard;


public:
    RankineMatGrad(int n, Domain *d);

    const char *giveClassName() const override { return "RankineMatGrad"; }
    const char *giveInputRecordName() const override  { return _IFT_RankineMatGrad_Name; }

    void initializeFrom(InputRecord &ir) override;
    bool hasMaterialModeCapability(MaterialMode mode) const override;
    Interface *giveInterface(InterfaceType t) override {
        if ( t == GradientDamageMaterialExtensionInterfaceType ) {
            return static_cast< GradientDamageMaterialExtensionInterface * >( this );
        } else {
            return nullptr;
        }
    }

    void giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    void giveGradientDamageStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_du_NB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void giveGradientDamageStiffnessMatrix_du_BB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void giveGradientDamageStiffnessMatrix_dd_NN(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_dd_BB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_dd_BN(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;

    void computeLocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep) override;

    void giveNonlocalInternalForces_N_factor(double &answer, double nlddv, GaussPoint *gp, TimeStep *tStep) override;
    void giveNonlocalInternalForces_B_factor(FloatArray &answer, const FloatArray &nlddv, GaussPoint *gp, TimeStep *tStep) override;


    void giveRealStressVectorGradientDamage(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep) override;

    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode, GaussPoint * gp,  TimeStep * tStep) const override;
    void givePlaneStressGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStressKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void giveInternalLength(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    double computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const override;
    double giveNonlocalCumPlasticStrain(GaussPoint *gp);
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);

    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

protected:
    double computeInternalLength(GaussPoint *gp);
    int giveDimension(GaussPoint *gp);

    double computeEikonalInternalLength_a(GaussPoint *gp);
    double computeEikonalInternalLength_b(GaussPoint *gp);
    double computeEikonalInternalLength_aPrime(GaussPoint *gp);
    double computeEikonalInternalLength_bPrime(GaussPoint *gp);

    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new RankineMatGradStatus(gp); }
};
} // end namespace oofem
#define RankineMatGrad_h
#endif
