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

#ifndef trabbonenlembed_h
#define trabbonenlembed_h

#include "trabboneembed.h"
#include "sm/Materials/structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "cltypes.h"

///@name Input fields for TrabBoneNLEmbed
//@{
#define _IFT_TrabBoneNLEmbed_Name "trabbonenlembed"
#define _IFT_TrabBoneNLEmbed_r "r"
#define _IFT_TrabBoneNLEmbed_m "m"
//@}

namespace oofem {
class GaussPoint;

/**
 * Trabecular bone nonlocal material status.
 */
class TrabBoneNLEmbedStatus : public TrabBoneEmbedStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:
    /// Equivalent strain for averaging.
    double localCumPlastStrainForAverage;

public:
    TrabBoneNLEmbedStatus(int n, Domain * d, GaussPoint * g);
    virtual ~TrabBoneNLEmbedStatus();

    void printOutputAt(FILE *file, TimeStep *tStep) override;

    /// Gives the local cumulative plastic strain.
    double giveLocalCumPlastStrainForAverage() { return localCumPlastStrainForAverage; }
    /// Sets the local cumulative plastic strain.
    void setLocalCumPlastStrainForAverage(double ls) { localCumPlastStrainForAverage = ls; }

    const char *giveClassName() const override { return "TrabBoneNLEmbedStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    Interface *giveInterface(InterfaceType it) override;
};

/**
 * Trabecular bone nonlocal material.
 */
class TrabBoneNLEmbed : public TrabBoneEmbed, public StructuralNonlocalMaterialExtensionInterface
{
protected:
    double R;
    double mParam;

public:
    TrabBoneNLEmbed(int n, Domain * d);
    virtual ~TrabBoneNLEmbed();

    const char *giveClassName() const override { return "TrabBoneNLEmbed"; }
    const char *giveInputRecordName() const override { return _IFT_TrabBoneNLEmbed_Name; }

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    Interface *giveInterface(InterfaceType) override;

    void computeCumPlastStrain(double &alpha, GaussPoint *gp, TimeStep *tStep) override;

    void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &strainVector, TimeStep *tStep) override;

    void computeLocalCumPlastStrain(double &alpha, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
    {
        TrabBoneEmbed :: computeCumPlastStrain(alpha, gp, tStep);
    }

    void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep) override;
    double computeWeightFunction(const FloatArray &src, const FloatArray &coord) override;

    int hasBoundedSupport() override { return 1; }

    /// Determines the width (radius) of limited support of weighting function.
    virtual void giveSupportRadius(double &radius) { radius = this->R; }

protected:
    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new TrabBoneNLEmbedStatus(1, TrabBoneEmbed :: domain, gp); }
};
} // end namespace oofem
#endif
