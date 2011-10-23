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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef MisesMatGrad_h

#include "misesmat.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "cltypes.h"

#include "sparsemtrx.h"
#include "dynalist.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {

/**
 * Gradient Mises maaterial status.
 */
class MisesMatGradStatus : public MisesMatStatus
{
protected:
    double localCumPlastStrainForAverage;

public:
    MisesMatGradStatus(int n, Domain *d, GaussPoint *g);
    ~MisesMatGradStatus();

    void printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    const char *giveClassName() const { return "MisesMatGradStatus"; }
    classType giveClassID() const { return MisesMatClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);
};


/**
 * Gradient Mises material.
 */
class MisesMatGrad : public MisesMat
{
protected:
    double R;
    double mParam;

public:
    MisesMatGrad(int n, Domain *d);
    ~MisesMatGrad();

    // definition
    const char *giveClassName() const { return "MisesMatGrad"; }
    classType giveClassID() const { return MisesMatClass; }
    const char *giveInputRecordName() const { return "MisesMatGrad"; }

    IRResultType initializeFrom(InputRecord *ir);
    int hasMaterialModeCapability(MaterialMode mode);

    void giveCharacteristicMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual void give1dStressStiffMtrx(FloatMatrix & answer,  MatResponseForm, MatResponseMode, GaussPoint * gp,  TimeStep * tStep);
    virtual void givePlaneStrainStiffMtrx(FloatMatrix & answer,  MatResponseForm, MatResponseMode, GaussPoint * gp,  TimeStep * tStep);
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer,  MatResponseForm, MatResponseMode, GaussPoint * gp,  TimeStep * tStep);

    void give1dKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStrainKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void give3dKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void give1dGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStrainGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void give3dGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void giveInternalLength(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void giveRealStressVector(FloatArray &answer,  MatResponseForm form, GaussPoint *gp, const FloatArray &strainVector, TimeStep *tStep);
    void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep);
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

protected:
    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new MisesMatGradStatus(1, MisesMat :: domain, gp); }
};
} // end namespace oofem
#define MisesMatGrad_h
#endif
