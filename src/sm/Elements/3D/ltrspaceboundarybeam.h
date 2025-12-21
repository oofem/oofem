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

#ifndef ltrspaceboundarybeam_h
#define ltrspaceboundarybeam_h

#include "sm/Elements/3D/ltrspaceboundary.h"

#define _IFT_LTRSpaceBoundaryBeam_Name "ltrspaceboundarybeam"
#define _IFT_LTRSpaceBoundaryBeam_Location "location"

namespace oofem {
class FEI3dTetLin;
class ParamKey;

/**
 * This class implements a linear tetrahedral four-node finite element.
 * Each node has 3 degrees of freedom. This element is used for 3D RVE analyses with Periodic Boundary Conditions.
 * At least one node is located at the image boundary.
 * These nodes are replaced with a periodic mirror nodes and a control node is used to impose the macroscopic (average) strain.
 * MACROSCOPIC INPUT: DEFORMATIONS AND CURVATURE (BEAM, 3 COMPONENTS: Exx Ezx Kxx)
 * Exx = du/dx - axial strain
 * Ezx = dw/dx - slope
 * Kxx = d^2(w)/dx^2 (Euler-Bernoulli) or d(phi)/dx (Timoshenko) - curvature
 *
 * @author: Adam Sciegaj
 */
class LTRSpaceBoundaryBeam : public LTRSpaceBoundary
{
protected:
    void computeTransformationMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    static ParamKey IPK_LTRSpaceBoundaryBeam_location;
public:
    LTRSpaceBoundaryBeam(int n, Domain *d);
    virtual ~LTRSpaceBoundaryBeam() { }

    int computeNumberOfDofs() override { return 15; };
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    // definition & identification
    void initializeFrom(InputRecord &ir, int priority) override;
    const char *giveInputRecordName() const override { return _IFT_LTRSpaceBoundaryBeam_Name; }
    const char *giveClassName() const override { return "LTRSpaceBoundaryBeam"; }
};
} // end namespace oofem
#endif // LTRSpaceBoundaryBeam_h
