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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef binghamfluid2_h
#define binghamfluid2_h

#include "fluiddynamicmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"

#include "matconst.h"
#include "structuralelement.h"
#include "matstatus.h"

namespace oofem {
class GaussPoint;

#define BINGHAM_DEFAULT_STRESS_GROWTH_RATE 400.0

/**
 * Class representing material status for Bingham material
 */
class BinghamFluidMaterial2Status : public FluidDynamicMaterialStatus
{
protected:
    /// Magnitude of deviatoric strains
    double devStrainMagnitude, temp_devStrainMagnitude;
    /// Magnitude of deviatoric stresses
    double devStressMagnitude, temp_devStressMagnitude;
    /// Deviatoric stresses and strains (reduced form).
    FloatArray deviatoricStrainVector, temp_deviatoricStrainVector;

public:
    /// Constructor - creates new BinghamFluidMaterial2Status with number n, belonging to domain d and IntegrationPoint g.
    BinghamFluidMaterial2Status(int n, Domain *d, GaussPoint *g);
    /// Destructor
    ~BinghamFluidMaterial2Status() { }

    void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    double giveTempDevStressMagnitude() const { return temp_devStressMagnitude; }
    double giveTempDevStrainMagnitude() const { return temp_devStrainMagnitude; }
    double giveDevStressMagnitude() const { return devStressMagnitude; }
    double giveDevStrainMagnitude() const { return devStrainMagnitude; }

    void letTempDevStrainMagnitudeBe(double _val) { temp_devStrainMagnitude = _val; }
    void letTempDevStressMagnitudeBe(double _val) { temp_devStressMagnitude = _val; }

    const FloatArray &giveDeviatoricStrainVector() { return deviatoricStrainVector; }
    const FloatArray &giveTempDeviatoricStrainVector() { return temp_deviatoricStrainVector; }
    void letTempDeviatoricStrainVectorBe(const FloatArray &v) { temp_deviatoricStrainVector = v; }

    const char *giveClassName() const { return "BinghamFluidMaterialStatus"; }
    classType giveClassID() const { return BinghamFluidMaterialStatusClass; }
};



/**
 * Constitutive model of Bingham fluid for concentrated suspensions and pastes.
 * This is the simplest two-constant model, with yield stress and viscosity as parameters.
 */
class BinghamFluidMaterial2 : public FluidDynamicMaterial
{
protected:
    /// Viscosity.
    double mu_0;
    /// Yield stress.
    double tau_0;
    double tau_c;
    double mu_inf;
    // Stress growth rate - parameter controlling the shape of regularized model.
    double stressGrowthRate;
public:

    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    BinghamFluidMaterial2(int n, Domain *d) : FluidDynamicMaterial(n, d) { mu_inf = 1.e6; stressGrowthRate=BINGHAM_DEFAULT_STRESS_GROWTH_RATE;}
    /// Destructor.
    ~BinghamFluidMaterial2() { }

    virtual void giveCharacteristicMatrix(FloatMatrix &answer,
                                          MatResponseForm form,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *atTime) { }
    virtual double giveCharacteristicValue(MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);
    virtual void computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep);

    virtual void giveDeviatoricStiffnessMatrix(FloatMatrix & answer, MatResponseMode, GaussPoint * gp,
                                               TimeStep * tStep);

    virtual double give(int aProperty, GaussPoint *gp);
    IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    virtual int hasMaterialModeCapability(MaterialMode mode);
    const char *giveClassName() const { return "BinghamFluidMaterial2"; }
    classType giveClassID() const { return BinghamFluidMaterialClass; }
    virtual int checkConsistency();
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    double computeActualViscosity(double Tau, double shearRate);
    double computeDevStrainMagnitude(MaterialMode mmode, const FloatArray &epsd);
    double computeDevStressMagnitude(MaterialMode mmode, const FloatArray &sigd);
    void computeDeviatoricStrain(FloatArray &answer, const FloatArray &eps, MaterialMode mmode);
    void computeDeviatoricStress(FloatArray &answer, const FloatArray &deps,
                                    double _nu, MaterialMode mmode);

    void __debug(GaussPoint *gp, TimeStep *atTime);
};
} // end namespace oofem
#endif // binghamfluid2_h
