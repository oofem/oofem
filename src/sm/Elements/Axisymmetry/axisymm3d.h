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

#ifndef axisymm3d_h
#define axisymm3d_h

#include "Elements/structural2delement.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "spatiallocalizer.h"

///@name Input fields for Axisymm3d
//@{
#define _IFT_Axisymm3d_Name "axisymm3d"
//#define _IFT_Axisymm3d_nipfish "nipfish"
//@}

namespace oofem {
class FEI2dTrLin;

/**
 * This class implements an triangular three-node finite element for axisymmetric continuum.
 * Each node has 2 degrees of freedom.
 *
 * Tasks:
 * - calculating its B,D,N matrices and dV.
 */
class Axisymm3d : public AxisymElement, public ZZNodalRecoveryModelInterface,
public NodalAveragingRecoveryModelInterface, public SPRNodalRecoveryModelInterface,
public SpatialLocalizerInterface
{
protected:
    static FEI2dTrLin interpolation;

    double area;

public:
    Axisymm3d(int n, Domain * d);
    virtual ~Axisymm3d();

    virtual double giveArea();
    virtual Interface *giveInterface(InterfaceType it);
    virtual FEInterpolation *giveInterpolation() const;

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type);
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
#endif

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                            InternalStateType type, TimeStep *tStep);

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();

    virtual const char *giveClassName() const { return "Axisymm3d"; }
    virtual const char *giveInputRecordName() const { return _IFT_Axisymm3d_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

};
} // end namespace oofem
#endif // axisymm3d_h
