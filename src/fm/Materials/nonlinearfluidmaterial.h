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

#ifndef nonlinearfluidmaterial_h
#define nonlinearfluidmaterial_h

#include "fm/Materials/fluiddynamicmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matconst.h"
#include "matstatus.h"

///@name Input fields for NonlinearFluidMaterial
//@{
#define _IFT_NonlinearFluidMaterial_Name "nonlinearfluid"
#define _IFT_NonlinearFluidMaterial_mu "mu"
#define _IFT_NonlinearFluidMaterial_alpha "alpha"
#define _IFT_NonlinearFluidMaterial_C "c"
//@}

namespace oofem {
class GaussPoint;

/**
 * Material status class for NonlinearFluidMaterial
 *
 * @author Carl Sandström
 */
class NonlinearFluidMaterialStatus : public FluidDynamicMaterialStatus
{
protected:
    FloatArrayF<6> temp_deviatoricStressVector;
    FloatArrayF<6> temp_deviatoricStrainVector;
    double temp_norm2 = 0.;

public:
    NonlinearFluidMaterialStatus(GaussPoint * g);

    void initTempStatus() override;

    void updateYourself(TimeStep *tStep) override;

    const FloatArrayF<6> &giveTempDeviatoricStressVector() { return temp_deviatoricStressVector; }
    const FloatArrayF<6> &giveTempDeviatoricStrainVector() { return temp_deviatoricStrainVector; }
    double giveTempStrainNorm2() { return temp_norm2; }
    void letTempDeviatoricStressVectorBe(const FloatArrayF<6> &v) { temp_deviatoricStressVector = v; }
    void letTempDeviatoricStrainVectorBe(const FloatArrayF<6> &v) { temp_deviatoricStrainVector = v; }
    void letTempStrainNorm2Be(double v) { temp_norm2 = v; }

    const char *giveClassName() const override { return "NonlinearFluidMaterialStatus"; }
};

/**
 * Constitutive model of a nonlinear fluid material where the deviatoric stress is defined as
 * @f[
 * \boldsymbol{\sigma}_{\text{dev}}=2\mu(1+C \mid\mid \boldsymbol{v} \otimes \boldsymbol{\nabla} \mid\mid^{\alpha})\boldsymbol{v} \otimes \boldsymbol{\nabla}
 * @f]
 * where @f$ C @f$ and @f$ \alpha @f$ are constants and @f$ \mu @f$ the viscosity.
 *
 * @author Carl Sandström
 */
class NonlinearFluidMaterial : public FluidDynamicMaterial
{
protected:
    /// Viscosity @f$ \mu @f$ of material.
    double viscosity;
    /// Material constant @f$ C @f$.
    double c;
    /// Material constant @f$ \alpha @f$.
    double alpha;

public:
    NonlinearFluidMaterial(int n, Domain * d) : FluidDynamicMaterial(n, d) { }

    FloatArrayF<6> computeDeviatoricStress3D(const FloatArrayF<6> &eps, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<6,6> computeTangent3D(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;

    double giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep) const override;
    double give(int aProperty, GaussPoint *) const override;

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    const char *giveClassName() const override { return "NewtonianFluidMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_NonlinearFluidMaterial_Name; }

    int checkConsistency() override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
};
} // end namespace oofem
#endif // nonlinearfluidmaterial_h
