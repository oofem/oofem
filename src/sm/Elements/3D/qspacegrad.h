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

#ifndef qspacegrad_h
#define qspacegrad_h

#include "sm/Elements/3D/qspace.h"
#include "sm/Elements/graddpelement.h"

#define _IFT_QSpaceGrad_Name "qspacegrad"

namespace oofem {
class FEI3dHexaLin;

/**
 * Quadratic 3d  20 - node element with quadratic approximation of displacements and linear approximation of gradient
 *
 * @author Ladislav Svoboda
 */
class QSpaceGrad : public QSpace, public GradDpElement
{
protected:
    static FEI3dHexaLin interpolation_lin;

public:
    QSpaceGrad(int n, Domain * d);
    virtual ~QSpaceGrad() { }

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_QSpaceGrad_Name; }
    const char *giveClassName() const override { return "QSpaceGrad"; }
    int computeNumberOfDofs() override { return 68; }
    MaterialMode giveMaterialMode() override { return _3dMat; }

protected:
    void computeGaussPoints() override;
    void computeNkappaMatrixAt(GaussPoint *gp, FloatArray &answer) override;
     void computeBkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    StructuralElement *giveStructuralElement() override { return this; }
    NLStructuralElement *giveNLStructuralElement() override { return this; }
};
}
#endif // end namespace oofem
