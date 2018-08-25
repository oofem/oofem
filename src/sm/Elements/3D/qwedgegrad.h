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

#ifndef qwedgegrad_h
#define qwedgegrad_h

#include "sm/Elements/3D/qwedge.h"
#include "sm/Elements/graddpelement.h"

#define _IFT_QWedgeGrad_Name "qwedgegrad"

namespace oofem {
class FEI3dWedgeLin;

/**
 * Quadratic 3D element
 * @author L. Svoboda
 */
class QWedgeGrad : public QWedge, public GradDpElement
{
protected:
    static FEI3dWedgeLin interpolation_lin;

public:
    QWedgeGrad(int, Domain *);
    virtual ~QWedgeGrad() { }

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_QWedgeGrad_Name; }
    const char *giveClassName() const override { return "QWedgeGrad"; }
    int computeNumberOfDofs() override { return 51; }
    MaterialMode giveMaterialMode() override { return _3dMat; }

protected:
    void computeGaussPoints() override;
    void computeNkappaMatrixAt(GaussPoint *gp, FloatArray &answer) override;
     void computeBkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    StructuralElement *giveStructuralElement() override { return this; }
    NLStructuralElement *giveNLStructuralElement() override { return this; }

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override { GradDpElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord); }
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override { GradDpElement :: computeStiffnessMatrix(answer, rMode, tStep); }
};
}
#endif // end namespace oofem
