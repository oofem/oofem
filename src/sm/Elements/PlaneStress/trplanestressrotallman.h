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

#ifndef trplanestressrotallman_h
#define trplanestressrotallman_h

#include "sm/Elements/PlaneStress/trplanstrss.h"


///@name Input fields for TrPlaneStrRotAllman
//@{
#define _IFT_TrPlanestressRotAllman_Name "trplanestressrotallman"
//@}

namespace oofem {
class FEI2dTrQuad;

/**
 * Class implements an triangular three-node  plane-
 * stress elasticity finite element with independentvertex rotations.
 * Each node has 3 degrees of freedom.
 * For reference, see:
 * Allman, D.J.
 * A compatible triangular element including vertex rotations for plane elasticity analysis.
 * Computers & Structures Vol. 19, No. 1-2, pp. 1-8, 1984.
 */
class TrPlanestressRotAllman : public TrPlaneStress2d
{
protected:
    static FEI2dTrQuad qinterpolation; // quadratic interpolation for constructing shape functons

public:
    TrPlanestressRotAllman(int, Domain *);
    virtual ~TrPlanestressRotAllman() { }

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override;

    double giveArea() override;
    virtual void computeLocalNodalCoordinates(std::vector< FloatArray > &lxy);
    /**
     * Computes the stiffness matrix stabilization of zero energy mode (equal rotations)
     *
     * @param answer Computed stiffness matrix (symmetric).
     * @param rMode Response mode.
     * @param tStep Time step.
     */
    void computeStiffnessMatrixZeroEnergyStabilization(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
public:
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_TrPlanestressRotAllman_Name; }
    const char *giveClassName() const override { return "TrPlanestressRotAllman"; }
    void initializeFrom(InputRecord &ir) override;
    MaterialMode giveMaterialMode() override { return _PlaneStress; }
    integrationDomain giveIntegrationDomain() const override { return _Triangle; }
    /** Computes the stiffness matrix of receiver. Overloaded to add stabilization of zero-energy mode (equal rotations) */
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void computeGaussPoints() override;
    int computeNumberOfDofs() override { return 9; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;

    Interface *giveInterface(InterfaceType interface) override;

    void computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global) override;
    void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const override;
    int SPRNodalRecoveryMI_giveNumberOfIP() override;
};
} // end namespace oofem
#endif //  trplanestressrotallman_h
