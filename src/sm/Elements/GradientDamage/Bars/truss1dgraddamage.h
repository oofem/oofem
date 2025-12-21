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

#ifndef truss1dgraddamage_h
#define truss1dgraddamage_h


#include "../sm/Elements/Bars/truss1d.h"
#include "../sm/Elements/GradientDamage/graddamageelement.h"

#define _IFT_Truss1dGradDamage_Name "truss1dgraddamage"

namespace oofem {
class FEI1dLin;

/**
 * This class implements a two-node truss bar element for one-dimensional
 * analysis.
 */
 class Truss1dGradDamage : public Truss1d, public GradientDamageElement 
{
protected:
  static IntArray locationArray_u;
  static IntArray locationArray_d;
public:
    Truss1dGradDamage(int n, Domain * d);
    virtual ~Truss1dGradDamage() {;}

    const char *giveInputRecordName() const override { return _IFT_Truss1dGradDamage_Name; }
    const char *giveClassName() const override { return "Truss1dGradDamage"; }

    MaterialMode giveMaterialMode() override { return _1dMat; }
    void initializeFrom(InputRecord &ir, int priority) override;
    int computeNumberOfDofs() override { return 4; }

protected:
    void computeBdMatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    void computeNdMatrixAt(GaussPoint *gp, FloatArray &answer) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;
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
#endif // truss1dgraddamage_h
