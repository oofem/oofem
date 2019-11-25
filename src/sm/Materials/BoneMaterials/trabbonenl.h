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

#ifndef trabbonenl_h
#define trabbonenl_h

#include "trabbonematerial.h"
#include "sm/Materials/structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "cltypes.h"

///@name Input fields for TrabBoneNL
//@{
#define _IFT_TrabBoneNL_Name "trabbonenl"
#define _IFT_TrabBoneNL_r "r"
#define _IFT_TrabBoneNL_m "m"
//@}

namespace oofem {
/**
 * Trabecular bone nonlocal material status
 */
class TrabBoneNLStatus : public TrabBoneMaterialStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:
    /// Equivalent strain for averaging.
    double localCumPlastStrainForAverage = 0.;

public:
    TrabBoneNLStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    double giveLocalCumPlastStrainForAverage() const { return localCumPlastStrainForAverage; }
    void setLocalCumPlastStrainForAverage(double ls) { localCumPlastStrainForAverage = ls; }

    const char *giveClassName() const override { return "TrabBoneNLStatus"; }

    void initTempStatus() override;

    void updateYourself(TimeStep *tStep) override;

    Interface *giveInterface(InterfaceType) override;
};


/**
 * Trabecular bone nonlocal material.
 */
class TrabBoneNL : public TrabBoneMaterial, public StructuralNonlocalMaterialExtensionInterface
{
protected:
    double R = 0.;
    double mParam = 0.;

public:
    TrabBoneNL(int n, Domain * d);

    const char *giveClassName() const override { return "TrabBoneNL"; }
    const char *giveInputRecordName() const override { return _IFT_TrabBoneNL_Name; }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    Interface *giveInterface(InterfaceType) override;

    double computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const override;

    FloatArrayF<1> giveRealStressVector_1d(const FloatArrayF<1> &strainVector, GaussPoint *gp, TimeStep *tStep) const override;

    double computeLocalCumPlastStrain(const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) const
    {
        return TrabBoneMaterial :: computeCumPlastStrain(gp, tStep);
    }

    void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep) const override;

    double computeWeightFunction(const FloatArray &src, const FloatArray &coord) const override;

    int hasBoundedSupport() const override { return 1; }

    double giveSupportRadius() const { return this->R; }

protected:
    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new TrabBoneNLStatus(gp); }
};
} // end namespace oofem
#endif
