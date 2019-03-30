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
    double localCumPlastStrainForAverage;

public:
    MisesMatGradStatus(GaussPoint *g);
    virtual ~MisesMatGradStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "MisesMatGradStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual double giveNonlocalCumulatedStrain() { return nonlocalDamageDrivingVariable; }
    virtual void setNonlocalCumulatedStrain(double nonlocalCumulatedStrain) { this->nonlocalDamageDrivingVariable = nonlocalCumulatedStrain; }
};


/**
 * Gradient Mises material.
 */
class MisesMatGrad : public MisesMat, GradientDamageMaterialExtensionInterface
{
protected:
    double L;
    double mParam;

public:
    MisesMatGrad(int n, Domain *d);
    virtual ~MisesMatGrad();

    // definition
    virtual const char *giveInputRecordName() const { return _IFT_MisesMatGrad_Name; }
    virtual const char *giveClassName() const { return "MisesMatGrad"; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual Interface *giveInterface(InterfaceType t) {
        if ( t == GradientDamageMaterialExtensionInterfaceType ) {
            return static_cast< GradientDamageMaterialExtensionInterface * >( this );
        } else {
            return NULL;
        }
    }

    virtual void giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual void give1dStressStiffMtrx(FloatMatrix & answer, MatResponseMode, GaussPoint * gp, TimeStep * tStep);
    virtual void givePlaneStrainStiffMtrx(FloatMatrix & answer, MatResponseMode, GaussPoint * gp, TimeStep * tStep);
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer, MatResponseMode, GaussPoint * gp, TimeStep * tStep);

    virtual void giveGradientDamageStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_dd_l(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveRealStressVectorGradientDamage(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep);


    void give1dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStrainKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void give3dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void give1dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStrainGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void give3dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void giveInternalLength(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);


    virtual void computeLocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep);
    virtual void giveNonlocalInternalForces_N_factor(double &answer, double nlddv, GaussPoint *gp, TimeStep *tStep);
    virtual void giveNonlocalInternalForces_B_factor(FloatArray &answer, const FloatArray &nlddv, GaussPoint *gp, TimeStep *tStep);

    virtual void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep);
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);
    //    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

protected:
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new MisesMatGradStatus(gp); }
};
} // end namespace oofem
#define MisesMatGrad_h
#endif
