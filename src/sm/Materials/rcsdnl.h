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

//   *********************************************************************************************
//   *** CLASS NONLOCAL ROTATING SMEARED CRACK MODEL WITH TRANSITION TO SCALAR DAMAGE ************
//   *********************************************************************************************

#ifndef rcsdnl_h
#define rcsdnl_h

#include "rcsde.h"
#include "Materials/structuralnonlocalmaterialext.h"
#include "../sm/Materials/structuralmaterial.h"
///@name Input fields for RCSDNLMaterial
//@{
#define _IFT_RCSDNLMaterial_Name "rcsdnl"
#define _IFT_RCSDNLMaterial_ft "ft"
#define _IFT_RCSDNLMaterial_sdtransitioncoeff "sdtransitioncoeff"
#define _IFT_RCSDNLMaterial_sdtransitioncoeff2 "sdtransitioncoeff2"
#define _IFT_RCSDNLMaterial_r "r"
#define _IFT_RCSDNLMaterial_ef "ef"
#define _IFT_RCSDNLMaterial_gf "gf"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to RCSDNLMaterial.
 */
class RCSDNLMaterialStatus : public RCSDEMaterialStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:
    FloatArray nonlocalStrainVector, tempNonlocalStrainVector, localStrainVectorForAverage;

public:
    RCSDNLMaterialStatus(int n, Domain * d, GaussPoint * g);
    virtual ~RCSDNLMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    const FloatArray &giveNonlocalStrainVector() { return nonlocalStrainVector; }
    const FloatArray &giveTempNonlocalStrainVector() { return tempNonlocalStrainVector; }
    void setTempNonlocalStrainVector(FloatArray ls) { tempNonlocalStrainVector = std :: move(ls); }

    const FloatArray &giveLocalStrainVectorForAverage() { return localStrainVectorForAverage; }
    void setLocalStrainVectorForAverage(FloatArray ls) { localStrainVectorForAverage = std :: move(ls); }

    // definition
    virtual const char *giveInputRecordName() const { return _IFT_RCSDNLMaterial_Name; }
    virtual const char *giveClassName() const { return "RCSDNLMaterialStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *);

    // saves current context(state) into stream
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual Interface *giveInterface(InterfaceType);
};


/**
 * This class implements a Nonlocal Rotating Crack Model with transition to scalar damage
 * for fracture in smeared fashion
 * Only material stiffness modification is required, no changes in mesh topology.
 * Model according to Milan Jirasek RC-SD model.
 */
class RCSDNLMaterial : public RCSDEMaterial, public StructuralNonlocalMaterialExtensionInterface
{
protected:
    /**
     * Nondimensional parameter controlling the transition between rc and sd model,
     * with respect to shear stifness degradation.
     */
    double SDTransitionCoeff2;
    /**
     * Interaction radius, related to the nonlocal characteristic length of material.
     */
    double R;
    /**
     * Strain at complete failure. For exponential law, ef is the strain at the intersection
     * of the horizontal axis with the tangent to the softening curve at peak stress.
     */
    double ef;

public:
    RCSDNLMaterial(int n, Domain * d);
    virtual ~RCSDNLMaterial();

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "RCSDNLMaterial"; }

    virtual Interface *giveInterface(InterfaceType t);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *,
                                      const FloatArray &, TimeStep *);

    /**
     * Implements the service updating local variables in given integration points,
     * which take part in nonlocal average process. Actually, no update is necessary,
     * because the value used for nonlocal averaging is strain vector used for nonlocal secant stiffness
     * computation. It is therefore necessary only to store local strain in corresponding status.
     * This service is declared at StructuralNonlocalMaterial level.
     * @param strainVector total strain vector in given integration point.
     * @param gp integration point to update.
     * @param tStep solution step indicating time of update.
     */
    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep);

    virtual double computeWeightFunction(const FloatArray &src, const FloatArray &coord);
    virtual int hasBoundedSupport() { return 1; }
    /**
     * Determines the width (radius) of limited support of weighting function
     */
    virtual void giveSupportRadius(double &radius) { radius = this->R; }

    virtual int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip);
    virtual int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip);
    virtual int estimatePackSize(DataStream &buff, GaussPoint *ip);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new RCSDNLMaterialStatus(1, RCSDEMaterial :: domain, gp); }

protected:
    virtual double giveCharacteristicElementLength(GaussPoint *gp, const FloatArray &) { return 1.0; }
    virtual double giveMinCrackStrainsForFullyOpenCrack(GaussPoint *gp, int i);
    virtual double computeStrength(GaussPoint *, double)  { return this->Ft; }
};
} // end namespace oofem
#endif // rcsdnl_h
