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

#ifndef TrabBoneGrad3D_h
#define TrabBoneGrad3D_h

#include "trabbone3d.h"
#include "sm/Materials/graddpmaterialextensioninterface.h"
#include "cltypes.h"

#define _IFT_TrabBoneGrad3D_Name "trabbonegrad3d"
#define _IFT_TrabBoneGrad3D_L "l"
#define _IFT_TrabBoneGrad3D_m "mParam"

namespace oofem {
class LinearElasticMaterial;

/**
 * Gradient bone damage-plastic material status.
 */
class TrabBoneGrad3DStatus : public TrabBone3DStatus, GradDpMaterialStatusExtensionInterface
{
protected:
    /// Equivalent strain for avaraging
    double nlKappa;
    /// Reference to the basic elastic material
    LinearElasticMaterial *linearElasticMaterial;

public:
    TrabBoneGrad3DStatus(int n, Domain *d, GaussPoint *g);
    virtual ~TrabBoneGrad3DStatus();

    void printOutputAt(FILE *file, TimeStep *tStep) override;

    const char *giveClassName() const override { return "TrabBoneGrad3DStatus"; }

    double giveNonlocalCumulatedStrain() override { return nonlocalCumulatedStrain; }
    void setNonlocalCumulatedStrain(double nonlocalCumulatedStrain) override { this->nonlocalCumulatedStrain = nonlocalCumulatedStrain; }

    void initTempStatus() override;

    void updateYourself(TimeStep *tStep) override;
};


/**
 * Gradient bone damage-plastic material model.
 */
class TrabBoneGrad3D : public TrabBone3D, GradDpMaterialExtensionInterface
{
protected:
    double L;
    double mParam;

public:
    TrabBoneGrad3D(int n, Domain *d);
    virtual ~TrabBoneGrad3D();

    const char *giveClassName() const override { return "TrabBoneGrad3D"; }
    const char *giveInputRecordName() const override { return _IFT_TrabBoneGrad3D_Name; }

    IRResultType initializeFrom(InputRecord *ir) override;
    int hasMaterialModeCapability(MaterialMode mode) override;
    Interface *giveInterface(InterfaceType t) override {
        if ( t == GradDpMaterialExtensionInterfaceType ) {
            return static_cast< GradDpMaterialExtensionInterface * >( this );
        } else {
            return nullptr;
        }
    }

    void givePDGradMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePDGradMatrix_ku(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePDGradMatrix_uk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePDGradMatrix_kk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePDGradMatrix_LD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;

    void giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void give3dMaterialStiffnessMatrix(FloatMatrix & answer, MatResponseMode, GaussPoint * gp, TimeStep * tStep) override;
    void give3dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void give3dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void giveInternalLength(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void giveRealStressVectorGrad(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep) override;
    void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep) override;
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);

    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new TrabBoneGrad3DStatus(1, TrabBone3D :: domain, gp); }
};
} // end namespace oofem

#endif
