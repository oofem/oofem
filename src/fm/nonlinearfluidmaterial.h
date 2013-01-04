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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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
#include "structuralelement.h"
#include "matstatus.h"

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
    NonlinearFluidMaterialStatus(int n, Domain *d, GaussPoint *g) ;

    ~NonlinearFluidMaterialStatus() { }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *);

    const FloatArray &giveDeviatoricStrainVector()         			  { return deviatoricStrainVector; }
    const FloatArray &giveTempDeviatoricStrainVector()     			  { return temp_deviatoricStrainVector; }
    void         	 letTempDeviatoricStrainVectorBe(const FloatArray &v) { temp_deviatoricStrainVector = v; }

    /// Returns "NonlinearFluidMaterialStatus" string - class name of the receiver.
    const char *giveClassName() const { return "NonlinearFluidMaterialStatus"; }

    /// Returns NonlinearFluidMaterialStatusClass - classType id of receiver.
    classType                giveClassID() const
    { return NonlinearFluidMaterialStatusClass; }
};

/**
 * Constitutive model of a nonlinear fluid material where the deviatoric stress is defined as
 *
 * $\sigma^{\mbox{dev}}=2\mu(1+C \mid\mid \mathbm{v} \otimes \mathbm{\nabla} \mid\mid^{\alpha})\mathbm{v} \otimes \mathbm{\nabla}$
 *
 * where $C$ and $\alpha$ are constants and $\mu$ the viscosity.
 *
 * @author Carl Sandström
 */
class NonlinearFluidMaterial : public FluidDynamicMaterial
{
protected:
	/// Viscosity $\mu$ of o material
    double viscosity;
    /// Material constant c
	double c;
	/// Material constant alpha
	double alpha;
public:

    NonlinearFluidMaterial(int n, Domain *d) : FluidDynamicMaterial(n, d) { }

    ~NonlinearFluidMaterial()                { }

    virtual void  giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);

    virtual double  giveCharacteristicValue(MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *atTime);

    virtual void computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep);

    virtual void giveDeviatoricStiffnessMatrix(FloatMatrix & answer, MatResponseMode, GaussPoint * gp,
                                               TimeStep * atTime);

    virtual double give(int aProperty, GaussPoint*);

    IRResultType initializeFrom(InputRecord *ir);

    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    virtual int hasMaterialModeCapability(MaterialMode mode);

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "NewtonianFluidMaterial"; }
    /// Returns classType id of receiver.
    classType giveClassID()         const { return NewtonianFluidMaterialClass; }

    virtual int    checkConsistency();

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
};

} // end namespace oofem
#endif // nonlinearfluidmaterial_h
