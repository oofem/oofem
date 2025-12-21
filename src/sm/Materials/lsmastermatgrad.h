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

#ifndef lsmastermatgrad_h
#define lsmastermatgrad_h

#include "lsmastermat.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "Materials/linearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "graddamagematerialextensioninterface.h"

namespace oofem {
class GaussPoint;
class Domain;


#define _IFT_LargeStrainMasterMaterialGrad_Name "lsmastermatgrad"


/**
 * This class implements an gradient version of LargeStrainMasterMaterial
 */
class LargeStrainMasterMaterialGrad : public LargeStrainMasterMaterial, GradientDamageMaterialExtensionInterface
{
public:
    LargeStrainMasterMaterialGrad(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;
    bool hasMaterialModeCapability(MaterialMode mode) const override;
    Interface *giveInterface(InterfaceType t) override {
        if ( t == GradientDamageMaterialExtensionInterfaceType ) {
            return static_cast< GradientDamageMaterialExtensionInterface * >( this );
        } else {
            return nullptr;
        }
    }

    void giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    void giveGradientDamageStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_dd_BB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;

    void computeLocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep) override;
    void giveNonlocalInternalForces_N_factor(double &answer, double nlddv, GaussPoint *gp, TimeStep *tStep) override;
    void giveNonlocalInternalForces_B_factor(FloatArray &answer, const FloatArray &nlddv, GaussPoint *gp, TimeStep *tStep) override;

    void give3dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void give3dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void giveInternalLength(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    const char *giveClassName() const override { return "LargeStrainMasterMaterialGrad"; }

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

    virtual void giveFirstPKStressVectorGrad(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalDamageDrivningVariable, TimeStep *tStep);
};
} // end namespace oofem
#endif // misesmat_h
