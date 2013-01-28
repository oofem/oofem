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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef trabbonenl_h
#define trabbonenl_h

#include "trabbonematerial.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "cltypes.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

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
    TrabBoneNLStatus(int n, Domain *d, GaussPoint *g);
    virtual ~TrabBoneNLStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    double giveLocalCumPlastStrainForAverage()     { return localCumPlastStrainForAverage; }
    void setLocalCumPlastStrainForAverage(double ls) { localCumPlastStrainForAverage = ls; }

    virtual const char *giveClassName() const { return "TrabBoneNLStatus"; }
    virtual classType   giveClassID() const { return TrabBoneMaterialStatusClass; }

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
    TrabBoneNL(int n, Domain *d);
    virtual ~TrabBoneNL();

    virtual const char *giveClassName() const { return "TrabBoneNL"; }
    virtual classType   giveClassID()   const { return TrabBoneNLClass; }
    virtual const char *giveInputRecordName() const { return "trabbonenl"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    virtual Interface *giveInterface(InterfaceType);

    virtual void computeCumPlastStrain(double &alpha, GaussPoint *gp, TimeStep *atTime);

    virtual void giveRealStressVector(FloatArray &answer,  MatResponseForm form, GaussPoint *gp, const FloatArray &strainVector, TimeStep *atTime);

    void computeLocalCumPlastStrain(double &alpha, const StrainVector &strain, GaussPoint *gp, TimeStep *atTime)
    {
        TrabBoneMaterial :: computeCumPlastStrain(alpha, gp, atTime);
    }

    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *atTime);

    virtual double computeWeightFunction(const FloatArray &src, const FloatArray &coord);

    virtual int hasBoundedSupport() { return 1; }

    virtual void giveSupportRadius(double &radius) { radius = this->R; }

protected:
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new TrabBoneNLStatus(1, TrabBoneMaterial :: domain, gp); }
};
} // end namespace oofem
#endif
