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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef rankinematnl_h
#define rankinematnl_h

#include "rankinemat.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "nonlocalmaterialext.h"
#include "cltypes.h"

#define _IFT_RankineMatNl_Name "rankmatnl"

namespace oofem {
/**
 * Rankine nonlocal material status.
 */
class RankineMatNlStatus : public RankineMatStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:
    /// Equivalent strain for averaging.
    double localCumPlasticStrainForAverage = 0.;

    /// For printing only
    double kappa_nl = 0.;
    double kappa_hat = 0.;

public:
    RankineMatNlStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    double giveLocalCumPlasticStrainForAverage() const { return localCumPlasticStrainForAverage; }
    void setLocalCumPlasticStrainForAverage(double ls) { localCumPlasticStrainForAverage = ls; }

    const char *giveClassName() const override { return "RankineMatNlStatus"; }

    void initTempStatus() override;

    void updateYourself(TimeStep *tStep) override;

    void setKappa_nl(double kap) { kappa_nl = kap; }
    void setKappa_hat(double kap) { kappa_hat = kap; }
    double giveKappa_nl() const { return kappa_nl; }
    double giveKappa_hat() const { return kappa_hat; }

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    Interface *giveInterface(InterfaceType) override;
};


/**
 * Rankine nonlocal material
 */
class RankineMatNl : public RankineMat, public StructuralNonlocalMaterialExtensionInterface,
public NonlocalMaterialStiffnessInterface
{
public:
    RankineMatNl(int n, Domain * d);

    const char *giveClassName() const override { return "RankineMatNl"; }
    const char *giveInputRecordName() const override { return _IFT_RankineMatNl_Name; }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    Interface *giveInterface(InterfaceType) override;

    /**
     * Computes the nonlocal cumulated plastic strain from its local form.
     * @param kappa return param, containing the corresponding cumulated plastic strain.
     * @param gp integration point.
     * @param tStep time step.
     */
    virtual double computeCumPlasticStrain(GaussPoint *gp, TimeStep *tStep) const;
    double computeDamage(GaussPoint *gp, TimeStep *tStep) const;
    //void modifyNonlocalWeightFunctionAround(GaussPoint *gp);
    double computeDistanceModifier(double damage) const;
    double computeLocalCumPlasticStrain(GaussPoint *gp, TimeStep *tStep) const
    {
        return RankineMat :: computeCumPlastStrain(gp, tStep);
    }


    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;
    //void givePlaneStrainStiffMtrx(FloatMatrix& answer, MatResponseMode,GaussPoint * gp,TimeStep * tStep) override;
    //void give3dMaterialStiffnessMatrix(FloatMatrix& answer,  MatResponseMode,GaussPoint* gp, TimeStep* tStep) override;

#ifdef __OOFEG
    // Plots the sparse structure of stiffness contribution.
    //void NonlocalMaterialStiffnessInterface_showSparseMtrxStructure(GaussPoint *gp, oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

    void NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx &dest, const UnknownNumberingScheme &s,
                                                              GaussPoint *gp, TimeStep *tStep) override;

    std :: vector< localIntegrationRecord > *NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp) override;

    /**
     * Computes the "local" part of nonlocal stiffness contribution assembled for given integration point.
     * @param gp Source integration point.
     * @param loc Local code numbers.
     * @param lcontrib "Local" contribution.
     * @param s Numbering scheme that determines assembly.
     * @param tStep Time step.
     * @return Nonzero if local point contributes (loading) or zero if not (unloading in elastic range, elastic).
     */
    int giveLocalNonlocalStiffnessContribution(GaussPoint *gp, IntArray &loc, const UnknownNumberingScheme &s,
                                               FloatArray &lcontrib, TimeStep *tStep);

    /**
     * Computes the "remote" part of nonlocal stiffness contribution assembled for given integration point.
     * @param gp Remote integration point.
     * @param rloc Remote element code numbers.
     * @param rcontrib "Remote" contribution.
     * @param s Numbering scheme that determines assembly.
     * @param tStep Time step.
     */
    void giveRemoteNonlocalStiffnessContribution(GaussPoint *gp, IntArray &rloc, const UnknownNumberingScheme &s,
                                                 FloatArray &rcontrib, TimeStep *tStep);

    // Computes elastic stiffness for normal stress components
    // @param answer result of size (3,3)
    // @param mode determines the MatResponseMode
    // @param gp integration point
    // @param tStep time step
    //  void giveNormalElasticStiffnessMatrix (FloatMatrix& answer, MatResponseMode rMode, GaussPoint*gp, TimeStep* tStep) ;

    FloatArrayF<3> giveRealStressVector_PlaneStress(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const override;

    // Computes 1D stress
    FloatArrayF<1> giveRealStressVector_1d(const FloatArrayF<1> &strain, GaussPoint *gp, TimeStep *tStep) const override;

    void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep) const override;

    /// Compute the factor that specifies how the interaction length should be modified (by eikonal nonlocal damage models)
    double giveNonlocalMetricModifierAt(GaussPoint *gp) const override;

    int hasBoundedSupport() const override { return 1; }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override;
    int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override;
    int estimatePackSize(DataStream &buff, GaussPoint *ip) override;

protected:
    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<RankineMatNlStatus>(gp); }
};
} // end namespace oofem
#endif
