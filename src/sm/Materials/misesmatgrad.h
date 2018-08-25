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

#ifndef MisesMatGrad_h

#include "sm/Materials/misesmat.h"
#include "graddpmaterialextensioninterface.h"
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
class MisesMatGradStatus : public MisesMatStatus, GradDpMaterialStatusExtensionInterface
{
protected:
    double localCumPlastStrainForAverage;

public:
    MisesMatGradStatus(int n, Domain * d, GaussPoint * g);
    virtual ~MisesMatGradStatus();

    void printOutputAt(FILE *file, TimeStep *tStep) override;

    const char *giveClassName() const override { return "MisesMatGradStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    double giveNonlocalCumulatedStrain() override { return nonlocalCumulatedStrain; }
    void setNonlocalCumulatedStrain(double nonlocalCumulatedStrain) override { this->nonlocalCumulatedStrain = nonlocalCumulatedStrain; }
};


/**
 * Gradient Mises material.
 */
class MisesMatGrad : public MisesMat, GradDpMaterialExtensionInterface
{
protected:
    double L;
    double mParam;

public:
    MisesMatGrad(int n, Domain * d);
    virtual ~MisesMatGrad();

    const char *giveInputRecordName() const override { return _IFT_MisesMatGrad_Name; }
    const char *giveClassName() const override { return "MisesMatGrad"; }

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

    void give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode, GaussPoint *gp, TimeStep *tStep) override;
    void givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseMode, GaussPoint *gp, TimeStep *tStep) override;
    void give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode, GaussPoint *gp, TimeStep *tStep) override;

    void givePDGradMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePDGradMatrix_ku(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePDGradMatrix_uk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePDGradMatrix_kk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePDGradMatrix_LD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveRealStressVectorGrad(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep) override;

    void give1dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStrainKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void give3dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void give1dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStrainGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void give3dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void giveInternalLength(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep) override;
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);

protected:
    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new MisesMatGradStatus(1, MisesMat :: domain, gp); }
};
} // end namespace oofem
#define MisesMatGrad_h
#endif
