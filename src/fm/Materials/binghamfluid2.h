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

#ifndef binghamfluid2_h
#define binghamfluid2_h

#include "fm/Materials/fluiddynamicmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matconst.h"
#include "matstatus.h"

///@name Input fields for BinghamFluidMaterial
//@{
#define _IFT_BinghamFluidMaterial2_Name "binghamfluid"
#define _IFT_BinghamFluidMaterial2_mu0 "mu0"
#define _IFT_BinghamFluidMaterial2_tau0 "tau0"
#define _IFT_BinghamFluidMaterial2_muinf "muinf"
#define _IFT_BinghamFluidMaterial2_stressGrowthRate "stressgrowthrate"
//@}

namespace oofem {

#define BINGHAM_DEFAULT_STRESS_GROWTH_RATE 400.0

/**
 * Class representing material status for Bingham material
 */
class BinghamFluidMaterial2Status : public FluidDynamicMaterialStatus
{
protected:
    /// Magnitude of deviatoric strains
    double devStrainMagnitude = 0., temp_devStrainMagnitude = 0.;
    /// Magnitude of deviatoric stresses
    double devStressMagnitude = 0., temp_devStressMagnitude = 0.;
    /// Deviatoric stresses and strains (reduced form).
    FloatArrayF<6> temp_deviatoricStrainVector;

public:
    /// Constructor - creates new BinghamFluidMaterial2Status with number n, belonging to domain d and IntegrationPoint g.
    BinghamFluidMaterial2Status(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    double giveTempDevStressMagnitude() const { return temp_devStressMagnitude; }
    double giveTempDevStrainMagnitude() const { return temp_devStrainMagnitude; }
    double giveDevStressMagnitude() const { return devStressMagnitude; }
    double giveDevStrainMagnitude() const { return devStrainMagnitude; }

    void letTempDevStrainMagnitudeBe(double _val) { temp_devStrainMagnitude = _val; }
    void letTempDevStressMagnitudeBe(double _val) { temp_devStressMagnitude = _val; }

    const FloatArrayF<6> &giveTempDeviatoricStrainVector() const { return temp_deviatoricStrainVector; }
    void letTempDeviatoricStrainVectorBe(const FloatArrayF<6> &v) { temp_deviatoricStrainVector = v; }

    const char *giveClassName() const override { return "BinghamFluidMaterialStatus"; }
};


/**
 * Constitutive model of Bingham fluid for concentrated suspensions and pastes.
 * This is the simplest two-constant model, with yield stress and viscosity as parameters.
 */
class BinghamFluidMaterial2 : public FluidDynamicMaterial
{
protected:
    /// Viscosity.
    double mu_0 = 0.;
    /// Yield stress.
    double tau_0 = 0.;
    double tau_c = 0.;
    double mu_inf = 1.e6;
    /// Stress growth rate - parameter controlling the shape of regularized model.
    double stressGrowthRate = BINGHAM_DEFAULT_STRESS_GROWTH_RATE;

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    BinghamFluidMaterial2(int n, Domain * d);

    FloatArrayF<6> computeDeviatoricStress3D(const FloatArrayF<6> &eps, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<6,6> computeTangent3D(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;

    double giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep) const override;
    double give(int aProperty, GaussPoint *gp) const override;
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    const char *giveClassName() const override { return "BinghamFluidMaterial2"; }
    const char *giveInputRecordName() const override { return _IFT_BinghamFluidMaterial2_Name; }
    int checkConsistency() override;
    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

protected:
    double computeActualViscosity(double tau, double shearRate) const;
    static double computeDevStrainMagnitude(const FloatArrayF<6> &epsd);
    static double computeDevStressMagnitude(const FloatArrayF<6> &sigd);
    static FloatArrayF<6> computeDeviatoricStrain(const FloatArrayF<6> &eps);
    static FloatArrayF<6> computeDeviatoricStress(const FloatArrayF<6> &deps, double nu);

#if 0
    void __debug(GaussPoint *gp, TimeStep *tStep);
#endif
};
} // end namespace oofem
#endif // binghamfluid2_h
