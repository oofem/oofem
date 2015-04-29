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

#ifndef qtrplstr_h
#define qtrplstr_h

#include "Elements/structural2delement.h"
#include "ErrorEstimators/directerrorindicatorrc.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

#define _IFT_QTrPlaneStress2d_Name "qtrplstr"

namespace oofem {
class FEI2dTrQuad;

/**
 * This class implements a quadratic triangular 6-node plane-
 * stress elasticity finite element. Each node has 2 degrees of freedom.
 */
class QTrPlaneStress2d : public PlaneStressElement, public SpatialLocalizerInterface,
public SPRNodalRecoveryModelInterface
{
protected:
    static FEI2dTrQuad interpolation;

public:
    QTrPlaneStress2d(int n, Domain * d);
    virtual ~QTrPlaneStress2d() { }

    virtual Interface *giveInterface(InterfaceType it);
    virtual FEInterpolation *giveInterpolation() const;

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawSpecial(oofegGraphicContext &gc, TimeStep  *tStep);
#endif

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_QTrPlaneStress2d_Name; }
    virtual const char *giveClassName() const { return "QTrPlaneStress2d"; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();

protected:
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 4; }

};
} // end namespace oofem
#endif // qtrplstr_h
