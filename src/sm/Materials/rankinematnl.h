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

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    double giveLocalCumPlasticStrainForAverage() { return localCumPlasticStrainForAverage; }
    const FloatArray *giveLTangentContrib();
    void setLocalCumPlasticStrainForAverage(double ls) { localCumPlasticStrainForAverage = ls; }

    virtual const char *giveClassName() const { return "RankineMatNlStatus"; }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *tStep);

    void setKappa_nl(double kap) { kappa_nl = kap; }
    void setKappa_hat(double kap) { kappa_hat = kap; }
    double giveKappa_nl() { return kappa_nl; }
    double giveKappa_hat() { return kappa_hat; }

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual Interface *giveInterface(InterfaceType);
};


/**
 * Rankine nonlocal material
 */
class RankineMatNl : public RankineMat, public StructuralNonlocalMaterialExtensionInterface,
public NonlocalMaterialStiffnessInterface
{
public:
    RankineMatNl(int n, Domain * d);
    virtual ~RankineMatNl() {
        ;
    }

    virtual const char *giveClassName() const { return "RankineMatNl"; }
    virtual const char *giveInputRecordName() const { return _IFT_RankineMatNl_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual Interface *giveInterface(InterfaceType);

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


    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep);
    //virtual void givePlaneStrainStiffMtrx(FloatMatrix& answer, MatResponseMode,GaussPoint * gp,TimeStep * tStep);
    //virtual void give3dMaterialStiffnessMatrix(FloatMatrix& answer,  MatResponseMode,GaussPoint* gp, TimeStep* tStep);

#ifdef __OOFEG
    // Plots the sparse structure of stiffness contribution.
    //virtual void NonlocalMaterialStiffnessInterface_showSparseMtrxStructure(GaussPoint *gp, oofegGraphicContext &gc, TimeStep *tStep);
#endif

    virtual void NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx &dest, const UnknownNumberingScheme &s,
                                                                      GaussPoint *gp, TimeStep *tStep);

    virtual std :: list< localIntegrationRecord > *NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp);

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

    virtual void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &strainVector, TimeStep *tStep);

    // Computes 1D stress
    virtual void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &strainVector, TimeStep *tStep);

    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep);

    virtual int hasBoundedSupport() { return 1; }

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip);
    virtual int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip);
    virtual int estimatePackSize(DataStream &buff, GaussPoint *ip);

protected:
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new RankineMatNlStatus(1, RankineMat :: domain, gp); }
};
} // end namespace oofem
#endif
