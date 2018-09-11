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

#include <memory>

///@name Input fields for TwoFluidMaterial
//@{
#define _IFT_TwoFluidMaterial_Name "twofluidmat"
#define _IFT_TwoFluidMaterial_mat "mat"
//@}

namespace oofem {
class GaussPoint;

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
    /// Destructor.
    virtual ~TwoFluidMaterial() { }

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    void computeDeviatoricStress3D(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep) override;
    void computeTangent3D(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;

    double giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep) override;
    double give(int aProperty, GaussPoint *gp) override;
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    const char *giveClassName() const override { return "TwoFluidMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_TwoFluidMaterial_Name; }
    int checkConsistency() override;
    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

protected:
    FluidDynamicMaterial *giveMaterial(int i) const;
    double giveTempVOF(GaussPoint *gp);
};


class TwoFluidMaterialStatus : public FluidDynamicMaterialStatus
{
protected:
    std :: unique_ptr< GaussPoint >slaveGp0;
    std :: unique_ptr< GaussPoint >slaveGp1;

public:
    /// Constructor
    TwoFluidMaterialStatus(int n, Domain * d, GaussPoint * g, const IntArray & slaveMaterial);
    /// Destructor
    virtual ~TwoFluidMaterialStatus() { }

    void printOutputAt(FILE *file, TimeStep *tStep) override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
    const char *giveClassName() const override { return "TwoFluidMaterialStatus"; }

    GaussPoint *giveSlaveGaussPoint0() { return this->slaveGp0.get(); }
    GaussPoint *giveSlaveGaussPoint1() { return this->slaveGp1.get(); }
};
} // end namespace oofem
#endif // twofluidmaterial_h
