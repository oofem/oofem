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

#ifndef rankinematnl_h
#define rankinematnl_h

#include "rankinemat.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
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
    double localCumPlasticStrainForAverage;

    /// For printing only
    double kappa_nl;
    double kappa_hat;

public:
    RankineMatNlStatus(int n, Domain * d, GaussPoint * g);
    virtual ~RankineMatNlStatus();

    void printOutputAt(FILE *file, TimeStep *tStep) override;

    double giveLocalCumPlasticStrainForAverage() { return localCumPlasticStrainForAverage; }
    const FloatArray *giveLTangentContrib();
    void setLocalCumPlasticStrainForAverage(double ls) { localCumPlasticStrainForAverage = ls; }

    const char *giveClassName() const override { return "RankineMatNlStatus"; }

    void initTempStatus() override;

    void updateYourself(TimeStep *tStep) override;

    void setKappa_nl(double kap) { kappa_nl = kap; }
    void setKappa_hat(double kap) { kappa_hat = kap; }
    double giveKappa_nl() { return kappa_nl; }
    double giveKappa_hat() { return kappa_hat; }

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
    virtual ~RankineMatNl() { }

    const char *giveClassName() const override { return "RankineMatNl"; }
    const char *giveInputRecordName() const override { return _IFT_RankineMatNl_Name; }

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    Interface *giveInterface(InterfaceType) override;

    /**
     * Computes the nonlocal cumulated plastic strain from its local form.
     * @param kappa return param, containing the corresponding cumulated plastic strain.
     * @param gp integration point.
     * @param tStep time step.
     */
    virtual void computeCumPlasticStrain(double &kappa, GaussPoint *gp, TimeStep *tStep);
    double computeDamage(GaussPoint *gp, TimeStep *tStep);
    void modifyNonlocalWeightFunctionAround(GaussPoint *gp);
    double computeDistanceModifier(double damage);
    void computeLocalCumPlasticStrain(double &kappa, GaussPoint *gp, TimeStep *tStep)
    {
        RankineMat :: computeCumPlastStrain(kappa, gp, tStep);
    }


    void givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) override;
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

    void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &strainVector, TimeStep *tStep) override;

    // Computes 1D stress
    void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &strainVector, TimeStep *tStep) override;

    void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep) override;

    int hasBoundedSupport() override { return 1; }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override;
    int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override;
    int estimatePackSize(DataStream &buff, GaussPoint *ip) override;

protected:
    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new RankineMatNlStatus(1, RankineMat :: domain, gp); }
};
} // end namespace oofem
#endif
