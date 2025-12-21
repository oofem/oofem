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

#ifndef qspacegraddamage_h
#define qspacegraddamage_h

#include "../sm/Elements/3D/qspace.h"
#include "../sm/Elements/GradientDamage/graddamageelement.h"

#define _IFT_QSpaceGradDamage_Name "qspacegraddamage"

namespace oofem {
class FEI3dHexaLin;

/**
 * Quadratic 3d  20 - node element with quadratic approximation of displacements and linear approximation of gradient damage driving variable
 *
 * @author Martin Horak
 */
class QSpaceGradDamage : public QSpace, public GradientDamageElement
{
protected:
    static FEI3dHexaLin interpolation_lin;

public:
    QSpaceGradDamage(int n, Domain * d);
    virtual ~QSpaceGradDamage() { }

    void initializeFrom(InputRecord &ir, int prio) override;
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    void giveDofManDofIDMask_u(IntArray &answer) const override;
    void giveDofManDofIDMask_d(IntArray &answer) const override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_QSpaceGradDamage_Name; }
    const char *giveClassName() const override { return "QSpaceGradDamage"; }
    int computeNumberOfDofs() override { return 68; }
    MaterialMode giveMaterialMode() override { return _3dMat; }

protected:
    void computeGaussPoints() override;
    void computeNdMatrixAt(GaussPoint *gp, FloatArray &answer) override;
    void computeBdMatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    StructuralElement *giveStructuralElement() override { return this; }
    NLStructuralElement *giveNLStructuralElement() override { return this; }
    void giveLocationArray_u(IntArray &answer) override { }
    void giveLocationArray_d(IntArray &answer) override { }
};
}
#endif // end namespace oofem
