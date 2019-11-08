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

#ifndef qplanestraingrad_h
#define qplanestraingrad_h

#include "sm/Elements/PlaneStrain/qplanestrain.h"
#include "sm/Elements/graddpelement.h"

#define _IFT_QPlaneStrainGrad_Name "qplanestraingrad"

namespace oofem {
class FEI2dQuadLin;

class QPlaneStrainGrad : public QPlaneStrain, public GradDpElement
{
protected:
    static FEI2dQuadLin interpolation_lin;

public:
    QPlaneStrainGrad(int n, Domain * d);
    virtual ~QPlaneStrainGrad() { }

    void initializeFrom(InputRecord &ir) override;

    const char *giveInputRecordName() const override { return _IFT_QPlaneStrainGrad_Name; }
    const char *giveClassName() const override { return "QPlaneStrainGrad"; }

    MaterialMode giveMaterialMode() override { return _PlaneStrain; }
    int computeNumberOfDofs() override { return 20; }

protected:
     void computeBkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    void computeNkappaMatrixAt(GaussPoint *gp, FloatArray &answer) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override { GradDpElement :: computeStiffnessMatrix(answer, rMode, tStep); }
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override { GradDpElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord); }

    void computeGaussPoints() override;
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    StructuralElement *giveStructuralElement() override { return this; }
    NLStructuralElement *giveNLStructuralElement() override { return this; }
};
} // end namespace oofem
#endif // qplanestraingrad_h
