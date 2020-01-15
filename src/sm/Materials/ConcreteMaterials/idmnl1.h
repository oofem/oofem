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

#ifndef idmnl1_h
#define idmnl1_h

#include "sm/Materials/ConcreteMaterials/idm1.h"
#include "sm/Materials/structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"

///@name Input fields for IDNLMaterial
//@{
#define _IFT_IDNLMaterial_Name "idmnl1"
#define _IFT_IDNLMaterial_r "r"
//@}

namespace oofem {
class GaussPoint;

/**
 * This class implements associated Material Status to IDNLMaterial (Nonlocal isotropic damage).
 * Stores local equivalent strain for averaging.
 */
class IDNLMaterialStatus : public IsotropicDamageMaterial1Status, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:
    /// Equivalent strain for averaging.
    double localEquivalentStrainForAverage = 0.;

    /* // Variables used to track loading/reloading
     * public:
     * enum LastStateType {LST_elastic, LST_loading, LST_unloading};
     * LastStateType lst;
     */

public:
    /// Constructor.
    IDNLMaterialStatus(GaussPoint *g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    /// Returns the local  equivalent strain to be averaged.
    double giveLocalEquivalentStrainForAverage() { return localEquivalentStrainForAverage; }
    /// Sets the localEquivalentStrainForAverage to given value.
    void setLocalEquivalentStrainForAverage(double ls) { localEquivalentStrainForAverage = ls; }

    const char *giveClassName() const override { return "IDNLMaterialStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    /**
     * Interface requesting service.
     * In the case of nonlocal constitutive models,
     * the use of multiple inheritance is assumed. Typically, the class representing nonlocal
     * constitutive model status is derived both from class representing local status and from class
     * NonlocalMaterialStatusExtensionInterface or from one of its derived classes
     * (which declare services and variables corresponding to specific analysis type).
     * @return In both cases, this function returns pointer to this object, obtained by
     * returning address of component or using pointer conversion from receiver to base class
     * NonlocalMaterialStatusExtensionInterface.
     */
    Interface *giveInterface(InterfaceType) override;
};


/**
 * This class implements a Nonlocal Isotropic Damage Model for Concrete in Tension
 * Model based on nonlocal averaging of equivalent strain.
 */
class IDNLMaterial : public IsotropicDamageMaterial1, public StructuralNonlocalMaterialExtensionInterface,
    public NonlocalMaterialStiffnessInterface
{
public:
    /// Constructor
    IDNLMaterial(int n, Domain *d);

    const char *giveClassName() const override { return "IDNLMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_IDNLMaterial_Name; }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    Interface *giveInterface(InterfaceType it) override;

    double computeEquivalentStrain(const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) const override;
    /**
     * Function used in the Stress based nonlocal variation.In this function the ratio of the first two
     * eigenvalues is and the angle of the first eigenvector with respect to the horizontal axis is calculated
     * @param[out] nx x-component of the first eigenvector of effective stress.
     * @param[out] ny y-component of the first eigenvector of effective stress.
     * @param[out] ratio Value of the ratio of the second over the first eigenvalue of the stress tensor (sigma2/sigma1)
     * @param gp Gauss Point whose nonlocal interactions domain is modified
     * @param flag showing whether stress based averaging is activated (flag=true).For zero strain states the stress-based averaging is deactivated (flag=false)
     */
    void computeAngleAndSigmaRatio(double &nx, double &ny, double &ratio, GaussPoint *gp, bool &flag) const;
    /**
     * Function used to compute the new weight based on stress-based averaging.
     * @param nx x-component of the first eigenvector of effective stress.
     * @param ny y-component of the first eigenvector of effective stress.
     * @param ratio Value of the ratio of the second over the first eigenvalue of the stress tensor (sigma2/sigma1).
     * @param gp Gauss Point whose nonlocal interactions domain is modified.
     * @param jGp Gauss Point which contributes to the nonlocal interactions domain of gp.
     * @param weight Original weight.
     * @return New weight based on stress-based averaging.
     */
    double computeStressBasedWeight(double &nx, double &ny, double &ratio, GaussPoint *gp, GaussPoint *jGp, double weight) const;
    double computeStressBasedWeightForPeriodicCell(double &nx, double &ny, double &ratio, GaussPoint *gp, GaussPoint *jGp) const;

    double computeLocalEquivalentStrain(const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) const
    { return IsotropicDamageMaterial1 :: computeEquivalentStrain(strain, gp, tStep); }

    void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep) const override;

    /// Compute the factor that specifies how the interaction length should be modified (by eikonal nonlocal damage models)
    double giveNonlocalMetricModifierAt(GaussPoint *gp) const override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

#ifdef __OOFEG
    /// Plots the sparse structure of stiffness contribution.
    void NonlocalMaterialStiffnessInterface_showSparseMtrxStructure(GaussPoint *gp, oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

    double computeDamageParam(double kappa, const FloatArray &strain, GaussPoint *gp) const override;

    /**@name Services required by NonlocalMaterialStiffnessInterface and related ones to support Nonlocal Stiffness*/
    //@{
    void NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx &dest, const UnknownNumberingScheme &s,
                                                              GaussPoint *gp, TimeStep *tStep) override;
    /**
     * Returns integration list of receiver. Contains localIntegrationRecord structures, containing
     * references to integration points and their weights that influence to nonlocal average in
     * receiver's associated integration point.
     */
    std :: vector< localIntegrationRecord > *NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp) override;
    /**
     * Computes the "local" part of nonlocal stiffness contribution assembled for given integration point.
     * @param gp Source integration point.
     * @param loc Local code numbers.
     * @param s Determines the equation numbering scheme.
     * @param lcontrib "Local" contribution.
     * @param tStep Time step.
     * @return Nonzero if local point contributes (loading) or zero if not (unloading in elastic range, elastic)
     */
    int giveLocalNonlocalStiffnessContribution(GaussPoint *gp, IntArray &loc, const UnknownNumberingScheme &s,
                                               FloatArray &lcontrib, TimeStep *tStep);
    /**
     * Computes the "remote" part of nonlocal stiffness contribution assembled for given integration point.
     * @param gp Remote integration point.
     * @param rloc Remote element code numbers.
     * @param s Determines the equation numbering scheme.
     * @param rcontrib "Remote" contribution.
     * @param tStep Time step.
     */
    void giveRemoteNonlocalStiffnessContribution(GaussPoint *gp, IntArray &rloc, const UnknownNumberingScheme &s,
                                                 FloatArray &rcontrib, TimeStep *tStep);
    /**
     * Computes elastic stiffness for normal stress components.
     * @param answer Result of size (3,3).
     * @param mode Determines the MatResponseMode.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    void giveNormalElasticStiffnessMatrix(FloatMatrix &answer,
                                          MatResponseMode mode,
                                          GaussPoint *gp, TimeStep *tStep);
    //@}

    int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override;
    int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override;
    int estimatePackSize(DataStream &buff, GaussPoint *ip) override;
    double predictRelativeComputationalCost(GaussPoint *gp) override;
    double predictRelativeRedistributionCost(GaussPoint *gp) override { return 1.0; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new IDNLMaterialStatus(gp); }

protected:
    void initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp) const override { }
};
} // end namespace oofem
#endif // idmnl1_h
