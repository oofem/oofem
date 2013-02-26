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

#ifndef TrabBoneGrad3D_h
#define TrabBoneGrad3D_h

#include "trabbone3d.h"

namespace oofem {

class LinearElasticMaterial;

/**
 * Gradient bone damage-plastic material status.
 */
class TrabBoneGrad3DStatus : public TrabBone3DStatus
{
protected:
    /// Equivalent strain for avaraging
    double nlKappa;
     /// Reference to the basic elastic material
    LinearElasticMaterial *linearElasticMaterial;

public:
    TrabBoneGrad3DStatus(int n, Domain *d, GaussPoint *g);
    virtual ~TrabBoneGrad3DStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual const char *giveClassName() const { return "TrabBoneGrad3DStatus"; }
    virtual classType giveClassID() const { return TrabBoneGrad3DClass; }

    double giveNlKappa(){return nlKappa;}
    void setNlKappa(double kappa){nlKappa = kappa;}
    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *);
};


/**
 * Gradient bone damage-plastic material model.
 */
class TrabBoneGrad3D : public TrabBone3D
{
protected:
    double l;
    double mParam;

public:
    TrabBoneGrad3D(int n, Domain *d);
    virtual ~TrabBoneGrad3D();

    virtual const char *giveClassName() const { return "TrabBoneGrad3D"; }
    virtual classType giveClassID() const { return TrabBoneGrad3DClass; }
    virtual const char *giveInputRecordName() const { return "TrabBoneGrad3D"; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int hasMaterialModeCapability(MaterialMode mode);

    void giveCharacteristicMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime);
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer,  MatResponseForm, MatResponseMode, GaussPoint * gp,  TimeStep * atTime);
    void give3dKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void give3dGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void giveInternalLength(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    virtual void giveRealStressVector(FloatArray &answer,  MatResponseForm form, GaussPoint *gp, const FloatArray &strainVector, TimeStep *atTime);
    virtual void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *atTime);
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);
    //LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

protected:
    // Creates the corresponding material status
    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new TrabBoneGrad3DStatus(1, TrabBone3D :: domain, gp); }
};
} // end namespace oofem

#endif
