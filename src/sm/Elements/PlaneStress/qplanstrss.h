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

#ifndef qplanstrss_h
#define qplanstrss_h

#include "Elements/structural2delement.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"

#define _IFT_QPlaneStress2d_Name "qplanestress2d"

namespace oofem {
class FEI2dQuadQuad;

/**
 * This class implements an Quadratic isoparametric eight-node quadrilateral plane-
 * stress elasticity finite element. Each node has 2 degrees of freedom.
 */
class QPlaneStress2d : public PlaneStressElement, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface
{
protected:
    static FEI2dQuadQuad interpolation;

public:
    QPlaneStress2d(int n, Domain * d);
    virtual ~QPlaneStress2d() { }

    virtual FEInterpolation *giveInterpolation() const;

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_QPlaneStress2d_Name; }
    virtual const char *giveClassName() const { return "QPlaneStress2d"; }
    
    virtual Interface *giveInterface(InterfaceType it);

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                            InternalStateType type, TimeStep *tStep);

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
#endif

};
} // end namespace oofem
#endif // qplanstrss_h
