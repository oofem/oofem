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

#ifndef twofluidmaterial_h
#define twofluidmaterial_h

#include "fm/Materials/fluiddynamicmaterial.h"
#include "intarray.h"
#include "matstatus.h"
#include "gausspoint.h"

#include <memory>
#include <array>

///@name Input fields for TwoFluidMaterial
//@{
#define _IFT_TwoFluidMaterial_Name "twofluidmat"
#define _IFT_TwoFluidMaterial_mat "mat"
//@}

namespace oofem {

/**
 * Material coupling the behavior of two particular materials based on
 * rule of mixture. The weighting factor is VOF fraction.
 */
class TwoFluidMaterial : public FluidDynamicMaterial
{
protected:
    IntArray slaveMaterial;

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    TwoFluidMaterial(int n, Domain * d) : FluidDynamicMaterial(n, d) { }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    FloatArrayF<6> computeDeviatoricStress3D(const FloatArrayF<6> &answer, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<6,6> computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    double giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep) const override;
    double give(int aProperty, GaussPoint *gp) const override;
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    const char *giveClassName() const override { return "TwoFluidMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_TwoFluidMaterial_Name; }
    int checkConsistency() override;
    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

protected:
    FluidDynamicMaterial *giveMaterial(int i) const;
    double giveTempVOF(GaussPoint *gp) const;
};


class TwoFluidMaterialStatus : public FluidDynamicMaterialStatus
{
protected:
    std::array<GaussPoint, 2> slaveGps;
    ///@todo This should technically suffice;
    //std::array<std::unique_ptr<MaterialStatus>, 2> slaveStatus;

public:
    /// Constructor
    TwoFluidMaterialStatus(GaussPoint * g, const std::array<Material*, 2> &slaveMaterial);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
    const char *giveClassName() const override { return "TwoFluidMaterialStatus"; }

    GaussPoint *giveSlaveGaussPoint0() { return &this->slaveGps[0]; }
    GaussPoint *giveSlaveGaussPoint1() { return &this->slaveGps[1]; }
};
} // end namespace oofem
#endif // twofluidmaterial_h
