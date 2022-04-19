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

#ifndef mazarsmodelnl_h
#define mazarsmodelnl_h

#include "mazarsmodel.h"
#include "sm/Materials/structuralnonlocalmaterialext.h"

///@name Input fields for MazarsNLMaterial
//@{
#define _IFT_MazarsNLMaterial_Name "mazarsmodelnl"
#define _IFT_MazarsNLMaterial_r "r"
//@}

namespace oofem {
class GaussPoint;

/**
 * This class implements associated Material Status to MazarsNLModel.
 */
class MazarsNLMaterialStatus : public MazarsMaterialStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:
    /// Equivalent strain for averaging.
    double localEquivalentStrainForAverage = 0.;

public:
    /// Constructor
    MazarsNLMaterialStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    /// Returns the local equivalent strain to be averaged.
    double giveLocalEquivalentStrainForAverage() { return localEquivalentStrainForAverage; }
    /// Sets the local equivalent strain for average to given value.
    void setLocalEquivalentStrainForAverage(double ls) { localEquivalentStrainForAverage = ls; }

    const char *giveClassName() const override { return "MazarsNLMaterialStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    /**
     * Interface requesting service.
     * In the case of nonlocal constitutive models,
     * the use of multiple inheritance is assumed. Typically, the class representing nonlocal
     * constitutive model status is derived both from class representing local status and from class
     * NonlocalMaterialStatusExtension or from one of its derived classes
     * (which declare services and variables corresponding to specific analysis type).
     * @return In both cases, this function returns pointer to this object, obtained by
     * returning adress of component or using pointer conversion from receiver to base class
     * NonlocalMaterialStatusExtension. If no nonlocal extension exists, NULL pointer is returned.
     */
    Interface *giveInterface(InterfaceType it) override;
};


/**
 * This class implements a Nonlocal Mazars Damage Model for Concrete
 * Model based on nonlocal averaging of equivalent strain.
 */
class MazarsNLMaterial : public MazarsMaterial, public StructuralNonlocalMaterialExtensionInterface
{
protected:
    /// Interaction radius, related to the nonlocal characteristic length of material.
    double R = 0.;

public:
    /// Constructor
    MazarsNLMaterial(int n, Domain * d);

    const char *giveClassName() const override { return "MazarsNLMaterial"; }

    void initializeFrom(InputRecord &ir) override;
    Interface *giveInterface(InterfaceType it) override;

    double computeEquivalentStrain(const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) const override;
    /**
     * Computes the equivalent local strain measure from given strain vector (full form).
     * @param[out] kappa Return parameter, containing the corresponding equivalent strain
     * @param strain Total strain vector in full form
     * @param gp Integration point.
     * @param tStep Time step.
     */
    double computeLocalEquivalentStrain(const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) const
    { return MazarsMaterial :: computeEquivalentStrain(strain, gp, tStep); }

    void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep) const override;
  double computeWeightFunction(const double cl, const FloatArray &src, const FloatArray &coord) const override;
    int hasBoundedSupport() const override { return 1; }
    /**
     * Determines the width (radius) of limited support of weighting function
     */
    virtual double giveSupportRadius() const { return this->R; }

    int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override;
    int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override;
    int estimatePackSize(DataStream &buff, GaussPoint *ip) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new MazarsNLMaterialStatus(gp); }

protected:
    void initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp) const override;
};
} // end namespace oofem
#endif // mazarsmodelnl_h
