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
#include "../sm/Materials/graddpmaterialextensioninterface.h"
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
    TrabBoneGrad3DStatus(int n, Domain * d, GaussPoint * g);
    virtual ~TrabBoneGrad3DStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual const char *giveClassName() const { return "TrabBoneGrad3DStatus"; }


    virtual double giveNonlocalCumulatedStrain() { return nonlocalCumulatedStrain; }
    virtual void setNonlocalCumulatedStrain(double nonlocalCumulatedStrain) { this->nonlocalCumulatedStrain = nonlocalCumulatedStrain; }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *);
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
    TrabBoneGrad3D(int n, Domain * d);
    virtual ~TrabBoneGrad3D();

    virtual const char *giveClassName() const { return "TrabBoneGrad3D"; }
    virtual const char *giveInputRecordName() const { return _IFT_TrabBoneGrad3D_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual Interface *giveInterface(InterfaceType t) {
        if ( t == GradDpMaterialExtensionInterfaceType ) {
            return static_cast< GradDpMaterialExtensionInterface * >(this);
        } else {
            return NULL;
        }
    }

    virtual void givePDGradMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void givePDGradMatrix_ku(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void givePDGradMatrix_uk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void givePDGradMatrix_kk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void givePDGradMatrix_LD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual void giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode, GaussPoint *gp, TimeStep *tStep);
    void give3dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void give3dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void giveInternalLength(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveRealStressVectorGrad(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep);
    virtual void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep);
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);
    //LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

protected:
    // Creates the corresponding material status
    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new TrabBoneGrad3DStatus(1, TrabBone3D :: domain, gp); }
};
} // end namespace oofem

#endif
