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

#ifndef idmgrad1_h
#define idmgrad1_h

#include "idm1.h"

namespace oofem {

/**
 * Gradient-enhancedl Isotropic Damage model for concrete in tension,
 */
class IDGMaterial : public IsotropicDamageMaterial1
{
protected:
    /// Length scale parameter
    double internalLength;
    /// Parameter specifiying which averaging type(0 - classical, 1 - distance based, 2 - stress based)
    int averType;
    /// Parameters specifying how the length scale parameter l is adjusted
    double beta,t;
    

public:
    /// Constructor
    IDGMaterial(int n, Domain *d);
    /// Destructor
    virtual ~IDGMaterial();

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "IDGMaterial"; }
    virtual classType giveClassID() const { return IDGMaterialClass; }
    virtual const char *giveInputRecordName() const { return "idmgrad1"; }
    virtual void giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,const FloatArray &totalStrain,TimeStep *atTime);
    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime);
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int hasMaterialModeCapability(MaterialMode mode);

    void giveCharacteristicMatrix(FloatMatrix &answer,  MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime);
    virtual void give1dStressStiffMtrx(FloatMatrix & answer,  MatResponseForm, MatResponseMode, GaussPoint * gp,  TimeStep * tStep);
    void give1dKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void give1dGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void givePlaneStressKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void givePlaneStressGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void giveInternalLength(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void giveInternalLengthDerivative(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);

    virtual void initDamaged(double kappa, FloatArray &strainVector, GaussPoint *gp);
 
    void computeEta(FloatMatrix &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime);
};


/**
 * Material status for gradient-enhanced isotropic damage model for concrete in tension.
 */
class IDGMaterialStatus : public IsotropicDamageMaterial1Status
{
public:
    IDGMaterialStatus(int n, Domain *d, GaussPoint *g);
    virtual ~IDGMaterialStatus();

    virtual const char *giveClassName() const { return "IDGMaterialStatus"; }
    virtual classType giveClassID() const { return IDGMaterialClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};

} // end namespace oofem
#endif // idmgrad1_h
