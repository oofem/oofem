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

#ifndef MisesMatGrad_h

#include "misesmat.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "graddpmaterialextensioninterface.h"
#include "cltypes.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

///@name Input fields for MisesMatGrad
//@{
#define _IFT_MisesMatGrad_Name "misesmatgrad"
#define _IFT_MisesMatGrad_r "r"
#define _IFT_MisesMatGrad_m "m"
//@}

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
    virtual ~MisesMatGradStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "MisesMatGradStatus"; }
    virtual classType giveClassID() const { return MisesMatClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);
};


/**
 * Gradient Mises material.
 */
class MisesMatGrad : public MisesMat, GradDpMaterialExtensionInterface
{
protected:
    double R;
    double mParam;

public:
    MisesMatGrad(int n, Domain *d);
    virtual ~MisesMatGrad();

    // definition
    virtual const char *giveInputRecordName() const { return _IFT_MisesMatGrad_Name; }
    virtual const char *giveClassName() const { return "MisesMatGrad"; }
    virtual classType giveClassID() const { return MisesMatClass; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int hasMaterialModeCapability(MaterialMode mode);
    
    virtual Interface *giveInterface(InterfaceType t) { if ( t == GradDpMaterialExtensionInterfaceType ) return static_cast< GradDpMaterialExtensionInterface* >(this); else return NULL; }

    virtual void giveCharacteristicMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual void give1dStressStiffMtrx(FloatMatrix & answer, MatResponseMode, GaussPoint * gp, TimeStep * tStep);
    virtual void givePlaneStrainStiffMtrx(FloatMatrix & answer, MatResponseMode, GaussPoint * gp, TimeStep * tStep);
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer, MatResponseMode, GaussPoint * gp, TimeStep * tStep);

    virtual void givePDGradMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void givePDGradMatrix_ku(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void givePDGradMatrix_uk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void givePDGradMatrix_kk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void givePDGradMatrix_LD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void give1dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStrainKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void give3dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void give1dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStrainGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void give3dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void giveInternalLength(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void giveRealStressVector(FloatArray &answer,  GaussPoint *gp, const FloatArray &strainVector, TimeStep *tStep);
    virtual void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep);
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

protected:
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new MisesMatGradStatus(1, MisesMat :: domain, gp); }
};
} // end namespace oofem
#define MisesMatGrad_h
#endif
