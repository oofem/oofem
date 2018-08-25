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
#include "graddpmaterialextensioninterface.h"
#include "cltypes.h"

///@name Input fields for RankineMatGrad
//@{
#define _IFT_RankineMatGrad_Name "rankmatgrad"
#define _IFT_RankineMatGrad_L "l"
#define _IFT_RankineMatGrad_m "m"
#define _IFT_RankineMatGrad_negligibleDamage "negligible_damage"
//@}

namespace oofem {
/**
 * Gradient rankine material status.
 */
class RankineMatGradStatus : public RankineMatStatus, public GradDpMaterialStatusExtensionInterface
{
protected:
    double kappa_nl;
    double kappa_hat;

public:
    RankineMatGradStatus(int n, Domain * d, GaussPoint * g);
    virtual ~RankineMatGradStatus() { }

    void printOutputAt(FILE *file, TimeStep *tStep) override;

    const char *giveClassName() const override { return "RankineMatGradStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void setKappa_nl(double kap) { kappa_nl = kap; }
    void setKappa_hat(double kap) { kappa_hat = kap; }
    double giveKappa_nl() { return kappa_nl; }
    double giveKappa_hat() { return kappa_hat; }
    double giveNonlocalCumulatedStrain() override { return nonlocalCumulatedStrain; }
    void setNonlocalCumulatedStrain(double nonlocalCumulatedStrain) override { this->nonlocalCumulatedStrain = nonlocalCumulatedStrain; }
};


/**
 * Gradient Rankine material.
 */
class RankineMatGrad : public RankineMat, GradDpMaterialExtensionInterface
{
protected:
    double L;
    double mParam;
    double negligible_damage;

public:
    RankineMatGrad(int n, Domain * d);
    virtual ~RankineMatGrad() { }

    const char *giveClassName() const override { return "RankineMatGrad"; }
    const char *giveInputRecordName() const override { return _IFT_RankineMatGrad_Name; }

    IRResultType initializeFrom(InputRecord *ir) override;
    int hasMaterialModeCapability(MaterialMode mode) override;
    Interface *giveInterface(InterfaceType t) override {
        if ( t == GradDpMaterialExtensionInterfaceType ) {
            return static_cast< GradDpMaterialExtensionInterface * >(this);
        } else {
            return nullptr;
        }
    }

    void giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;

    void givePDGradMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePDGradMatrix_uk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePDGradMatrix_ku(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePDGradMatrix_kk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePDGradMatrix_LD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveRealStressVectorGrad(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep) override;

    void givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode, GaussPoint *gp,  TimeStep *tStep) override;
    void givePlaneStressGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStressKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void giveInternalLength(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep) override;
    double giveNonlocalCumPlasticStrain(GaussPoint *gp);
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);

    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

protected:
    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new RankineMatGradStatus(1, RankineMat :: domain, gp); }
};
} // end namespace oofem
#define RankineMatGrad_h
#endif
