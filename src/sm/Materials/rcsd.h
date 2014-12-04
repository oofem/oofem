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

//   ************************************************************************************
//   *** CLASS ROTATING SMEARED CRACK MODEL WITH TRANSITION TO SCALAR DAMAGE ************
//   ************************************************************************************

#ifndef rcsd_h
#define rcsd_h

#include "rcm2.h"

///@name Input fields for RCSDMaterial
//@{
#define _IFT_RCSDMaterial_Name "rcsd"
#define _IFT_RCSDMaterial_sdtransitioncoeff "sdtransitioncoeff"
//@}

namespace oofem {
#define rcsd_Omega 300
#define pscm_SDTransitionCoeff 306
#define RCSD_DAMAGE_EPS 1.e-4

class GaussPoint;

/**
 * This class implements associated Material Status to RCSDMaterial
 */
class RCSDMaterialStatus : public RCM2MaterialStatus
{
public:
    enum rcsdMode { rcMode, sdMode };

protected:
    double maxEquivStrain, tempMaxEquivStrain;
    double damageCoeff, tempDamageCoeff;
    FloatMatrix Ds0;
    double damageStiffCoeff, depsf, depsp;
    rcsdMode mode, tempMode;

public:
    RCSDMaterialStatus(int n, Domain * d, GaussPoint * g);
    virtual ~RCSDMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    double giveTempMaxEquivStrain() { return tempMaxEquivStrain; }
    void setTempMaxEquivStrain(double val) { tempMaxEquivStrain = val; }
    double giveDamageStiffCoeff() { return damageStiffCoeff; }
    void setDamageStiffCoeff(double val) { damageStiffCoeff = val; }
    double giveTempDamageCoeff() { return tempDamageCoeff; }
    void setTempDamageCoeff(double val) { tempDamageCoeff = val; }
    const FloatMatrix *giveDs0Matrix() { return & Ds0; }
    void setDs0Matrix(FloatMatrix &mtrx) { Ds0 = mtrx; }
    double giveDamageEpsfCoeff() { return depsf; }
    void setDamageEpsfCoeff(double val) { depsf = val; }
    double giveDamageEpspCoeff() { return depsp; }
    void setDamageEpspCoeff(double val) { depsp = val; }

    rcsdMode giveTempMode() { return tempMode; }
    void setTempMode(rcsdMode mode) { tempMode = mode; }

    // query for non-tem variables (usefull for postprocessing)
    double giveMaxEquivStrain() { return maxEquivStrain; }
    double giveDamageCoeff() { return damageCoeff; }

    rcsdMode giveMode() { return mode; }

    // definition
    virtual const char *giveClassName() const { return "RCSDMaterialStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    // saves current context(state) into stream
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};


/**
 * This class implements a Rotating Crack Model with transition to scalar damage
 * for fracture in smeared fashion
 * ( only material stiffness modification is required, no changes in
 * mesh topology).
 * Model according to Milan Jirasek RC-SD model.
 */
class RCSDMaterial : public RCM2Material
{
protected:
    double SDTransitionCoeff;

public:
    RCSDMaterial(int n, Domain * d);
    virtual ~RCSDMaterial();

    // identification and auxiliary functions
    virtual const char *giveInputRecordName() const { return _IFT_RCSDMaterial_Name; }
    virtual const char *giveClassName() const { return "RCSDMaterial"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual double give(int aProperty, GaussPoint *gp);

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *,
                                      const FloatArray &, TimeStep *);

#ifdef __OOFEG
#endif

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new RCSDMaterialStatus(1, domain, gp); }

protected:
    double computeCurrEquivStrain(GaussPoint *, const FloatArray &, double, TimeStep *);
    // two functions used to initialize and updating temporary variables in
    // gp's status. These variables are used to control process, when
    // we try to find equilibrium state.

    virtual void giveEffectiveMaterialStiffnessMatrix(FloatMatrix &answer,
                                                      MatResponseMode rMode,
                                                      GaussPoint *gp, TimeStep *tStep);

    double computeDamageCoeff(double, double, double, double);
    virtual double giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp,
                                       double crackStrain, int i);
    //virtual double giveShearRetentionFactor(GaussPoint* gp, double eps_cr, int i);
    virtual double giveNormalCrackingStress(GaussPoint *gp, double eps_cr, int i);
    virtual double giveMinCrackStrainsForFullyOpenCrack(GaussPoint *gp, int i);
    //virtual void updateStatusForNewCrack( GaussPoint*, int, double);
    virtual double computeStrength(GaussPoint *, double);
    virtual int checkSizeLimit(GaussPoint *gp, double);
    ////
};
} // end namespace oofem
#endif // rcsd_h
