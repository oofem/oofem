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

#ifndef idmgrad_h
#define idmgrad_h

#include "Materials/ConcreteMaterials/idm1.h"
#include "../sm/Materials/graddamagematerialextensioninterface.h"

#define _IFT_IsotropicGradientDamageMaterial_Name "idmgrad"

#define _IFT_IsotropicGradientDamageMaterial_formulationType "formtype"
#define _IFT_IsotropicGradientDamageMaterial_di_rho "di_rho"
#define _IFT_IsotropicGradientDamageMaterial_di_eta "di_eta"

namespace oofem {
/**
 * Gradient-enhanced Isotropic Damage model for concrete in tension,
 */
class IsotropicGradientDamageMaterial : public IsotropicDamageMaterial1, GradientDamageMaterialExtensionInterface
{
protected:

    /**  Type characterizing the dependence of the internal lenght on variable of the state
     *  Note that the assigned numbers to enum values have to correspond to values
     *  used in initializeFrom to resolve internalLenghtDependence. If not, the consistency
     *  between initializeFrom and giveInputRecord methods is lost.
     */
    enum GradientDamageFormulationType {
        GDFT_Standard = 0,
        GDFT_DecreasingInteractions = 1,
        GDFT_Eikonal = 2
    };

    GradientDamageFormulationType gradientDamageFormulationType = GDFT_Standard;

    double di_rho = 0.;
    double di_eta = 0.;


public:
    /// Constructor
    IsotropicGradientDamageMaterial(int n, Domain *d);

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
    // identification and auxiliary functions
    const char *giveClassName() const override { return "IsotropicGradientDamageMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_IsotropicGradientDamageMaterial_Name; }
    void initializeFrom(InputRecord &ir) override;

    Interface *giveInterface(InterfaceType t) override {
        if ( t == GradientDamageMaterialExtensionInterfaceType ) {
            return static_cast< GradientDamageMaterialExtensionInterface * >( this );
        } else {
            return nullptr;
        }
    }
    bool hasMaterialModeCapability(MaterialMode mode) const override;

    void giveGradientDamageStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_dd_NN(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_dd_BB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_dd_BN(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;

    void giveRealStressVectorGradientDamage(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep) override;

    void giveStiffnessMatrix(FloatMatrix &answer,  MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeLocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep) override { }
    void giveNonlocalInternalForces_N_factor(double &answer, double nlddv, GaussPoint *gp, TimeStep *tStep) override;
    void giveNonlocalInternalForces_B_factor(FloatArray &answer, const FloatArray &nlddv, GaussPoint *gp, TimeStep *tStep) override;
protected:
    double computeInternalLength(GaussPoint *gp);
    int giveDimension(GaussPoint *gp);

    double computeEikonalInternalLength_a(GaussPoint *gp);
    double computeEikonalInternalLength_b(GaussPoint *gp);
    double computeEikonalInternalLength_aPrime(GaussPoint *gp);
    double computeEikonalInternalLength_bPrime(GaussPoint *gp);
};





/**
 * Material status for gradient-enhanced isotropic damage model for concrete in tension.
 */
class IsotropicGradientDamageMaterialStatus : public IsotropicDamageMaterial1Status, public GradientDamageMaterialStatusExtensionInterface
{
public:
    IsotropicGradientDamageMaterialStatus(GaussPoint *g);

    const char *giveClassName() const override { return "IsotropicGradientDamageMaterialStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *) override;
};
} // end namespace oofem
#endif // idmgrad_h
