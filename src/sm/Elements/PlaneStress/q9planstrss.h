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

#ifndef Q9PLANSTRSS_H_
#define Q9PLANSTRSS_H_

#include "Elements/structural2delement.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"

#define _IFT_Q9PlaneStress2d_Name "q9planestress2d"

namespace oofem {
class FEI2dQuadBiQuad;

/**
 * 9-node plane stress element.
 *
 * @date May 22, 2013
 * @author Erik Svenning
 */
class Q9PlaneStress2d : public PlaneStressElement, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface
{
protected:
    static FEI2dQuadBiQuad interpolation;

public:
    Q9PlaneStress2d(int n, Domain * d);
    virtual ~Q9PlaneStress2d() { }

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_Q9PlaneStress2d_Name; }
    virtual const char *giveClassName() const { return "Q9PlaneStress2d"; }
    virtual FEInterpolation *giveInterpolation() const;

    virtual Interface *giveInterface(InterfaceType it);

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                            InternalStateType type, TimeStep *tStep);

};
} // end namespace oofem
#endif // qplanstrss_h
