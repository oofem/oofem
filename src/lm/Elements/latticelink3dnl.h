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
 *               Copyright (C) 1993 - 2026   Borek Patzak
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef latticelink3dnl_h
#define latticelink3dnl_h

#include "latticelink3d.h"

///@name Input fields for LatticeLink3dNL
//@{
#define _IFT_LatticeLink3dNL_Name "latticelink3dnl"
//@}

namespace oofem {
/**
 * Geometrically nonlinear (corotational) version of LatticeLink3d for large rotations /
 * catenary action. Same input as LatticeLink3d.
 *
 * Kinematics (two rotation routes):
 *   - node-A rigid arm rotates with node A's spin: cA = R(spin_A) * (coordB - coordA);
 *   - bond frame (local axis 1 = bar direction) rotates with node B's spin only, since the
 *     bar is a material direction of the rebar and node B is the rebar node.
 * The jump uses cA (node A) and node B's translation (its arm is zero); the slip is that jump
 * projected onto the current bond frame R(spin_B)*triad0. Rotation is handled internally
 * (computeGtoLRotationMatrix returns false), mirroring Lattice3dNL.
 */
class LatticeLink3dNL : public LatticeLink3d
{
public:
    LatticeLink3dNL(int n, Domain *d) : LatticeLink3d(n, d) { }
    virtual ~LatticeLink3dNL() { }

    const char *giveInputRecordName() const override { return _IFT_LatticeLink3dNL_Name; }
    const char *giveClassName() const override { return "LatticeLink3dNL"; }

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override;

protected:
    bool computeGtoLRotationMatrix(FloatMatrix &answer) override { return false; }

    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;

    /// Geometric B in global coords: the node-A rigid-arm entries use the rotated arm cA.
    void computeNLBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep);

    /// Current bond frame (3x3, rows = current local axes in global): reference triad rotated
    /// by node-B spin only (bar direction is a material direction, tracked by the rebar node).
    void computeCurrentFrame(FloatMatrix &Qcur, const FloatArray &spinB);

    /// Rotated node-A rigid arm cA = R(spinA) * (coordB - coordA), global coords.
    void computeRotatedArm(FloatArray &cA, const FloatArray &spinA);
};
} // end namespace oofem
#endif
