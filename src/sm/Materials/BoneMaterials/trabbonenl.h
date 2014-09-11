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
#include "Materials/structuralnonlocalmaterialext.h"
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
    double localCumPlastStrainForAverage;

public:
    TrabBoneNLStatus(int n, Domain * d, GaussPoint * g);
    virtual ~TrabBoneNLStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    double giveLocalCumPlastStrainForAverage()     { return localCumPlastStrainForAverage; }
    void setLocalCumPlastStrainForAverage(double ls) { localCumPlastStrainForAverage = ls; }

    virtual const char *giveClassName() const { return "TrabBoneNLStatus"; }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *tStep);

    virtual Interface *giveInterface(InterfaceType);
};


/**
 * Trabecular bone nonlocal material.
 */
class TrabBoneNL : public TrabBoneMaterial, public StructuralNonlocalMaterialExtensionInterface
{
protected:
    double R;
    double mParam;

public:
    TrabBoneNL(int n, Domain * d);
    virtual ~TrabBoneNL();

    virtual const char *giveClassName() const { return "TrabBoneNL"; }
    virtual const char *giveInputRecordName() const { return _IFT_TrabBoneNL_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual Interface *giveInterface(InterfaceType);

    virtual void computeCumPlastStrain(double &alpha, GaussPoint *gp, TimeStep *tStep);

    virtual void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &strainVector, TimeStep *tStep);

    void computeLocalCumPlastStrain(double &alpha, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
    {
        TrabBoneMaterial :: computeCumPlastStrain(alpha, gp, tStep);
    }

    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep);

    virtual double computeWeightFunction(const FloatArray &src, const FloatArray &coord);

    virtual int hasBoundedSupport() { return 1; }

    virtual void giveSupportRadius(double &radius) { radius = this->R; }

protected:
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new TrabBoneNLStatus(1, TrabBoneMaterial :: domain, gp); }
};
} // end namespace oofem
#endif
