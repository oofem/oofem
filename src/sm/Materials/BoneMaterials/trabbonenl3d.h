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


#ifndef trabbonenl3d_h
#define trabbonenl3d_h

#include "trabbone3d.h"
#include "sm/Materials/structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "cltypes.h"

///@name Input fields for TrabBoneNL3D
//@{
#define _IFT_TrabBoneNL3D_Name "trabbonenl3d"
#define _IFT_TrabBoneNL3D_r "r"
#define _IFT_TrabBoneNL3D_m "m"
//@}

namespace oofem {
/**
 * Trabecular bone nonlocal material status.
 */
class TrabBoneNL3DStatus : public TrabBone3DStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:
    /// Equivalent strain for averaging.
    double localCumPlastStrainForAverage = 0.;

public:
    TrabBoneNL3DStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    double giveLocalCumPlastStrainForAverage() const { return localCumPlastStrainForAverage; }
    void setLocalCumPlastStrainForAverage(double ls) { localCumPlastStrainForAverage = ls; }

    const char *giveClassName() const override { return "TrabBoneNL3DStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;
    Interface *giveInterface(InterfaceType it) override;
};


/**
 * Trabecular bone nonlocal material model.
 */
class TrabBoneNL3D : public TrabBone3D,
public StructuralNonlocalMaterialExtensionInterface,
public NonlocalMaterialStiffnessInterface
{
protected:
    double R = 0.;
    double mParam = 0.;

public:
    TrabBoneNL3D(int n, Domain * d);

    const char *giveClassName() const override { return "TrabBoneNL3D"; }
    const char *giveInputRecordName() const override { return _IFT_TrabBoneNL3D_Name; }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    Interface *giveInterface(InterfaceType it) override;

    double computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const override;

    /**
     * Computes the local cumulated plastic strain from given strain vector (full form).
     * @param[out] kappa Return parameter, containing the corresponding cumulated plastic strain.
     * @param strain Total strain vector in full form.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    double computeLocalCumPlastStrain(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
    {
        return TrabBone3D :: computeCumPlastStrain(gp, tStep);
    }

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

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
     * @param s  Numbering scheme for unknowns.
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
     * @param s  Numbering scheme for unknowns.
     * @param tStep Time step.
     */
    void giveRemoteNonlocalStiffnessContribution(GaussPoint *gp, IntArray &rloc, const UnknownNumberingScheme &s,
                                                 FloatArray &rcontrib, TimeStep *tStep);

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override;

    void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep) const override;

  double computeWeightFunction(const double cl, const FloatArray &src, const FloatArray &coord) const override;

    int hasBoundedSupport() const override { return 1; }

    /// Determines the width (radius) of limited support of weighting function.
    virtual void giveSupportRadius(double &radius) { radius = this->R; }

    int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override;
    int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override;
    int estimatePackSize(DataStream &buff, GaussPoint *ip) override;

protected:
    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new TrabBoneNL3DStatus(gp); }
};
} // end namespace oofem
#endif
