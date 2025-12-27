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

#ifndef MisesMatGrad_h

#include "Materials/misesmat.h"
#include "graddamagematerialextensioninterface.h"
#include "cltypes.h"

///@name Input fields for MisesMatGrad
//@{
#define _IFT_MisesMatGrad_Name "misesmatgrad"
#define _IFT_MisesMatGrad_l "l"
#define _IFT_MisesMatGrad_m "m"
//@}

namespace oofem {
/**
 * Gradient Mises maaterial status.
 */
class MisesMatGradStatus : public MisesMatStatus, GradientDamageMaterialStatusExtensionInterface
{
protected:
    double localCumPlastStrainForAverage = 0.;

public:
    MisesMatGradStatus(GaussPoint *g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    // definition
    const char *giveClassName() const override { return "MisesMatGradStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    double giveNonlocalCumulatedStrain() const { return nonlocalDamageDrivingVariable; }
    void setNonlocalCumulatedStrain(double nonlocalCumulatedStrain) { this->nonlocalDamageDrivingVariable = nonlocalCumulatedStrain; }
};


/**
 * Gradient Mises material.
 */
class MisesMatGrad : public MisesMat, GradientDamageMaterialExtensionInterface
{
protected:
    double L = 0.;
    double mParam = 0.;

public:
    MisesMatGrad(int n, Domain *d);

    // definition
    const char *giveInputRecordName() const override { return _IFT_MisesMatGrad_Name; }
    const char *giveClassName() const override { return "MisesMatGrad"; }

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

    FloatMatrixF<1,1> give1dStressStiffMtrx(MatResponseMode, GaussPoint * gp, TimeStep * tStep) const override;
    FloatMatrixF<4,4> givePlaneStrainStiffMtrx(MatResponseMode, GaussPoint * gp, TimeStep * tStep) const override;
    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode, GaussPoint * gp, TimeStep * tStep) const override;

    void giveGradientDamageStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_dd_BB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveRealStressVectorGradientDamage(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep) override;


    void give1dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStrainKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void give3dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void give1dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStrainGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void give3dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void giveInternalLength(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);


    void computeLocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep) override;
    void giveNonlocalInternalForces_N_factor(double &answer, double nlddv, GaussPoint *gp, TimeStep *tStep) override;
    void giveNonlocalInternalForces_B_factor(FloatArray &answer, const FloatArray &nlddv, GaussPoint *gp, TimeStep *tStep) override;

    double computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const override;
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);
    //    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

protected:
    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<MisesMatGradStatus>(gp); }
};
} // end namespace oofem
#define MisesMatGrad_h
#endif
