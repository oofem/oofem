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

#include "../sm/Elements/structuralelement.h"
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

    virtual FEInterpolation *giveInterpolation() const;
    virtual FEInterpolation *giveInterpolation(DofIDItem id) const;

    virtual MaterialMode giveMaterialMode()  { return _Unknown; }

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_LineDistributedSpring_Name; }
    virtual const char *giveClassName() const { return "LineDistributedSpring"; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int computeNumberOfDofs() { return this->dofs.giveSize(); }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;

    virtual double computeVolumeAround(GaussPoint *gp);

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }
    void giveInternalForcesVector(FloatArray &answer,
                                  TimeStep *tStep, int useUpdatedGpRecord);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual Interface *giveInterface(InterfaceType it);
    int checkConsistency();
    void printOutputAt(FILE *File, TimeStep *tStep);

protected:
    virtual void computeGaussPoints();
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);

    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP() { return this->numberOfGaussPoints; }
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType() { return SPRPatchType_2dxy; }
};
} // end namespace oofem
#endif // linedistributedspring_H
