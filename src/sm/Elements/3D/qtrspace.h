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

#ifndef qtrspace_h
#define qtrspace_h

#include "sm/Elements/structural3delement.h"
#include "sm/ErrorEstimators/huertaerrorestimator.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

#define _IFT_QTRSpace_Name "qtrspace"

namespace oofem {
class FEI3dTetQuad;

/**
 * This class implements an Quadratic 3d  10 - node
 * elasticity finite element.
 *
 * Each node has 3 degrees of freedom.
 * One single additional attribute is needed for Gauss integration purpose :
 * 'jacobianMatrix'. This 3x3 matrix contains polynomials.
 * TASKS :
 * - calculating its Gauss points ;
 * - calculating its B,D,N matrices and dV.
 */
class QTRSpace : public Structural3DElement, public SPRNodalRecoveryModelInterface, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface
{
protected:
    static FEI3dTetQuad interpolation;

public:
    QTRSpace(int, Domain *);
    virtual ~QTRSpace() { }

    FEInterpolation *giveInterpolation() const override;

    void initializeFrom(InputRecord &ir) override;

    Interface *giveInterface(InterfaceType) override;
    int testElementExtension(ElementExtension ext) override { return ( ( ext == Element_SurfaceLoadSupport ) ? 1 : 0 ); }
    int giveNumberOfIPForMassMtrxIntegration() override { return 5; }

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    int SPRNodalRecoveryMI_giveNumberOfIP() override;
    SPRPatchType SPRNodalRecoveryMI_givePatchType() override;

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep) override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_QTRSpace_Name; }
    const char *giveClassName() const override { return "QTRSpace"; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_tetra_2;}

};
} // end namespace oofem
#endif
