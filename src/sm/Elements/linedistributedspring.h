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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#ifndef linedistributedspring_H
#define linedistributedspring_H

#include "sm/Elements/structuralelement.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

#define _IFT_LineDistributedSpring_Name "linedistributedspring"
#define _IFT_LineDistributedSpring_Dofs "dofs"
#define _IFT_LineDistributedSpring_Stifnesses "k"

namespace oofem {
class FEI3dLineLin;

/**
 * This class implements two-node subsoil element with linear interpolation.
 * In each node, a user selected DOFs can be interpolated.
 * At present only a linear spring behaviour is implemented directly at element level.
 * Can be generalized to nonlinear one, but then a material model should be used instead.
 *
 */
class LineDistributedSpring : public StructuralElement,
public ZZNodalRecoveryModelInterface,
public SPRNodalRecoveryModelInterface
{
protected:
    static FEI3dLineLin interp_lin;
    IntArray dofs;
    FloatArray springStiffnesses;

public:
    LineDistributedSpring(int n, Domain * d);
    virtual ~LineDistributedSpring() { }

    FEInterpolation *giveInterpolation() const override;
    FEInterpolation *giveInterpolation(DofIDItem id) const override;

    MaterialMode giveMaterialMode() override { return _Unknown; }

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_LineDistributedSpring_Name; }
    const char *giveClassName() const override { return "LineDistributedSpring"; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_line_1;}

    void initializeFrom(InputRecord &ir) override;

    int computeNumberOfDofs() override { return this->dofs.giveSize(); }
    void giveDofManDofIDMask(int inode, IntArray &) const override;

    double computeVolumeAround(GaussPoint *gp) override;

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override
    { computeLumpedMassMatrix(answer, tStep); }
    void giveInternalForcesVector(FloatArray &answer,
                                  TimeStep *tStep, int useUpdatedGpRecord) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    Interface *giveInterface(InterfaceType it) override;
    int checkConsistency() override;
    void printOutputAt(FILE *File, TimeStep *tStep) override;

protected:
    void computeGaussPoints() override;
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode) override;
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;

    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    int SPRNodalRecoveryMI_giveNumberOfIP() override { return this->numberOfGaussPoints; }
    SPRPatchType SPRNodalRecoveryMI_givePatchType() override { return SPRPatchType_2dxy; }
};
} // end namespace oofem
#endif // linedistributedspring_H
