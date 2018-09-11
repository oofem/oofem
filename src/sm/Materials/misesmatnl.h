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

#ifndef misesmatnl_h

#include "sm/Materials/misesmat.h"
#include "sm/Materials/structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "cltypes.h"

///@name Input fields for MisesMatNl
//@{
#define _IFT_MisesMatNl_Name "misesmatnl"
#define _IFT_MisesMatNl_averagingtype "averagingtype"
#define _IFT_MisesMatNl_exp "exp"
#define _IFT_MisesMatNl_rf "rf"
//@}

namespace oofem {
/**
 * Mises Nonlocal material status.
 * @author Milan
 */
class MisesMatNlStatus : public MisesMatStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:
    // STATE VARIABLE DECLARATION
    // Equivalent strain for avaraging
    double localCumPlasticStrainForAverage;

public:
    MisesMatNlStatus(int n, Domain * d, GaussPoint * g);
    virtual ~MisesMatNlStatus();

    void printOutputAt(FILE *file, TimeStep *tStep) override;

    // STATE VARIABLE
    // declare state variable access and modification methods
    double giveLocalCumPlasticStrainForAverage() { return localCumPlasticStrainForAverage; }
    const FloatArray *giveLTangentContrib();
    void setLocalCumPlasticStrainForAverage(double ls) { localCumPlasticStrainForAverage = ls; }

    const char *giveClassName() const override { return "MisesMatNlStatus"; }

    void initTempStatus() override;

    void updateYourself(TimeStep *tStep) override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    Interface *giveInterface(InterfaceType) override;
};


/**
 * Mises nonlocal material.
 * @author Milan
 */
class MisesMatNl : public MisesMat, public StructuralNonlocalMaterialExtensionInterface,
public NonlocalMaterialStiffnessInterface
{
protected:
    double Rf;
    double exponent;
    int averType;

public:
    MisesMatNl(int n, Domain * d);
    virtual ~MisesMatNl();

    const char *giveClassName() const override { return "MisesMatNl"; }
    const char *giveInputRecordName() const override { return _IFT_MisesMatNl_Name; }

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    Interface *giveInterface(InterfaceType) override;

    /**
     * Computes the nonlocal cumulated plastic strain from its local form.
     * @param[out] kappa Return param, containing the corresponding cumulated plastic strain.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void computeCumPlasticStrain(double &kappa, GaussPoint *gp, TimeStep *tStep);
    double computeDamage(GaussPoint *gp, TimeStep *tStep);
    void modifyNonlocalWeightFunctionAround(GaussPoint *gp);
    double computeDistanceModifier(double damage);
    void computeLocalCumPlasticStrain(double &kappa, GaussPoint *gp, TimeStep *tStep)
    {
        MisesMat :: computeCumPlastStrain(kappa, gp, tStep);
    }

    void give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) override;
    //void givePlaneStrainStiffMtrx(FloatMatrix& answer, MatResponseMode, GaussPoint *gp,TimeStep *tStep) override;
    //void give3dMaterialStiffnessMatrix(FloatMatrix& answer, MatResponseMode, GaussPoint *gp, TimeStep *tStep) override;

    void NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx &dest, const UnknownNumberingScheme &s,
                                                              GaussPoint *gp, TimeStep *tStep) override;

    std :: vector< localIntegrationRecord > *NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp) override;

    /**
     * Computes the "local" part of nonlocal stiffness contribution assembled for given integration point.
     * @param gp Source integration point
     * @param loc Local code numbers
     * @param lcontrib "Local" contribution
     * @param s Numbering scheme that determines assembly.
     * @param tStep Time step.
     * @return Nonzero if local point contributes (loading) or zero if not (unloading in elastic range, elastic)
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

    void giveRealStressVector_3d(FloatArray &answer,  GaussPoint *gp, const FloatArray &strainVector, TimeStep *tStep) override;
    void giveRealStressVector_1d(FloatArray &answer,  GaussPoint *gp, const FloatArray &strainVector, TimeStep *tStep) override;

    void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep) override;

    int hasBoundedSupport() override { return 1; }

    int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override;
    int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override;
    int estimatePackSize(DataStream &buff, GaussPoint *ip) override;

protected:
    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new MisesMatNlStatus(1, MisesMat :: domain, gp); }
};
} // end namespace oofem
#define misesmatnl_h
#endif
