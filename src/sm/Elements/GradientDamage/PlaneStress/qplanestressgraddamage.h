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

#ifndef qplanestressgraddamage_h
#define qplanestressgraddamage_h

#include "../sm/Elements/PlaneStress/qplanstrss.h"
#include "../sm/Elements/GradientDamage/graddamageelement.h"

#define _IFT_QPlaneStressGradDamage_Name "qplanestressgraddamage"

namespace oofem {
class FEI2dQuadLin;

class QPlaneStressGradDamage : public QPlaneStress2d, public GradientDamageElement
{
protected:
      static IntArray locationArray_u;
      static IntArray locationArray_d;
      static FEI2dQuadLin interpolation_lin;

public:
    QPlaneStressGradDamage(int n, Domain * d);
    virtual ~QPlaneStressGradDamage() { }

    void initializeFrom(InputRecord &ir, int priority) override;

    const char *giveInputRecordName() const override { return _IFT_QPlaneStressGradDamage_Name; }
    const char *giveClassName() const override { return "QPlaneStressGradDamage"; }
    MaterialMode giveMaterialMode() override { return _PlaneStress; }
    int computeNumberOfDofs() override { return 20; }

protected:
    void computeBdMatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    void computeNdMatrixAt(GaussPoint *gp, FloatArray &answer) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override { GradientDamageElement :: computeStiffnessMatrix(answer, rMode, tStep); }
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override { GradientDamageElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord); }

    void computeGaussPoints() override;
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    void giveDofManDofIDMask_u(IntArray &answer) const override;
    void giveDofManDofIDMask_d(IntArray &answer) const override;

    
    StructuralElement *giveStructuralElement() override { return this; }
    NLStructuralElement *giveNLStructuralElement() override { return this; }
    void giveLocationArray_u(IntArray &answer) override;
    void giveLocationArray_d(IntArray &answer) override;
    void postInitialize() override;

};
} // end namespace oofem
#endif // qplanestressgrad_h
