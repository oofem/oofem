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

#ifndef nonlinearfluidmaterial_h
#define nonlinearfluidmaterial_h

#include "fluiddynamicmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "matconst.h"
#include "matstatus.h"

///@name Input fields for NonlinearFluidMaterial
//@{
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

    FloatArray deviatoricStrainVector, temp_deviatoricStrainVector;  // reduced form

public:
    NonlinearFluidMaterialStatus(int n, Domain *d, GaussPoint *g);

    virtual ~NonlinearFluidMaterialStatus() { }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *);

    const FloatArray &giveDeviatoricStrainVector() { return deviatoricStrainVector; }
    const FloatArray &giveTempDeviatoricStrainVector() { return temp_deviatoricStrainVector; }
    void  letTempDeviatoricStrainVectorBe(const FloatArray &v) { temp_deviatoricStrainVector = v; }

    virtual const char *giveClassName() const { return "NonlinearFluidMaterialStatus"; }
    virtual classType giveClassID() const { return NonlinearFluidMaterialStatusClass; }
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
    NonlinearFluidMaterial(int n, Domain *d) : FluidDynamicMaterial(n, d) { }

    virtual ~NonlinearFluidMaterial() { }

    virtual void giveCharacteristicMatrix(FloatMatrix &answer,
                                          MatResponseForm form,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *atTime);

    virtual double giveCharacteristicValue(MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);

    virtual void computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep);

    virtual void giveDeviatoricStiffnessMatrix(FloatMatrix & answer, MatResponseMode, GaussPoint * gp,
                                               TimeStep * atTime);

    virtual double give(int aProperty, GaussPoint *);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual const char *giveClassName() const { return "NewtonianFluidMaterial"; }
    virtual classType giveClassID() const { return NewtonianFluidMaterialClass; }

    virtual int checkConsistency();

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
};
} // end namespace oofem
#endif // nonlinearfluidmaterial_h
