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

#ifndef newtonianfluid_h
#define newtonianfluid_h

#include "fm/Materials/fluiddynamicmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matconst.h"
#include "matstatus.h"

///@name Input fields for NewtonianFluidMaterial
//@{
#define _IFT_NewtonianFluidMaterial_Name "newtonianfluid"
#define _IFT_NewtonianFluidMaterial_mu "mu"
//@}

namespace oofem {
class GaussPoint;

/**
 * Constitutive model of Newtonian fluid
 */
class NewtonianFluidMaterial : public FluidDynamicMaterial
{
protected:
    double viscosity;

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    NewtonianFluidMaterial(int n, Domain * d) : FluidDynamicMaterial(n, d) { }

    double giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep) const override;

    FloatArrayF<6> computeDeviatoricStress3D(const FloatArrayF<6> &eps, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<6,6> computeTangent3D(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;

    double give(int aProperty, GaussPoint *gp) const override;
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    const char *giveClassName() const override { return "NewtonianFluidMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_NewtonianFluidMaterial_Name; }
    int checkConsistency() override;
    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;
};
} // end namespace oofem
#endif // newtonianfluid_h
