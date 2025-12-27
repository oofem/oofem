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

#ifndef quad1planestraingraddamage_h
#define quad1planestraingraddamage_h

#include "../sm/Elements/PlaneStrain/quad1planestrain.h"
#include "../sm/Elements/GradientDamage/graddamageelement.h"

#define _IFT_Quad1PlaneStrainGradDamage_Name "quad1planestraingraddamage"

namespace oofem {

class Quad1PlaneStrainGradDamage : public Quad1PlaneStrain, public GradientDamageElement
{
protected:
  static IntArray locationArray_u;
  static IntArray locationArray_d;

public:
    Quad1PlaneStrainGradDamage(int n, Domain * d);
    virtual ~Quad1PlaneStrainGradDamage() { }

    //void initializeFrom(InputRecord &ir) override;

    const char *giveInputRecordName() const override { return _IFT_Quad1PlaneStrainGradDamage_Name; }
    const char *giveClassName() const override { return "Quad1PlaneStrainGradDamage"; }

    MaterialMode giveMaterialMode() override { return _PlaneStrain; }
    int computeNumberOfDofs() override { return 12; }

protected:
    void computeBdMatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    void computeNdMatrixAt(GaussPoint *gp, FloatArray &answer) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override { GradientDamageElement :: computeStiffnessMatrix(answer, rMode, tStep); }
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override { GradientDamageElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord); }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    void giveDofManDofIDMask_u(IntArray &answer) const override;
    void giveDofManDofIDMask_d(IntArray &answer) const override;
    
    StructuralElement *giveStructuralElement() override { return this; }
    NLStructuralElement *giveNLStructuralElement() override { return this; }
    void postInitialize() override;
    void giveLocationArray_u(IntArray &answer) override;
    void giveLocationArray_d(IntArray &answer) override;
};
} // end namespace oofem
#endif // quad1planestraingrad_h
