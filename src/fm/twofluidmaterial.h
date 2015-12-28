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

#include "fluiddynamicmaterial.h"
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

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual void computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep);
    virtual void giveDeviatoricStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual double giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep);
    virtual double give(int aProperty, GaussPoint *gp);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual const char *giveClassName() const { return "TwoFluidMaterial"; }
    virtual const char *giveInputRecordName() const { return _IFT_TwoFluidMaterial_Name; }
    virtual int checkConsistency();
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

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

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual const char *giveClassName() const { return "TwoFluidMaterialStatus"; }

    GaussPoint *giveSlaveGaussPoint0() { return this->slaveGp0.get(); }
    GaussPoint *giveSlaveGaussPoint1() { return this->slaveGp1.get(); }
};
} // end namespace oofem
#endif // twofluidmaterial_h
