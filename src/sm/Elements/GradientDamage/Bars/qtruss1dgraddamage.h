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

#ifndef qtruss1dgraddamage_h
#define qtruss1dgraddamage_h

#include "../sm/Elements/Bars/qtruss1d.h"
#include "../sm/Elements/GradientDamage/graddamageelement.h"

#define _IFT_QTruss1dGradDamage_Name "qtruss1dgraddamage"

namespace oofem {
class FEI1dLin;

/**
 * This class implements a three-node gradient truss bar element for one-dimensional
 * analysis.
 */
class QTruss1dGradDamage : public QTruss1d, public GradientDamageElement
{
protected:
    static FEI1dLin interpolation_lin;

public:
    QTruss1dGradDamage(int n, Domain * d);
    virtual ~QTruss1dGradDamage() { }

    virtual const char *giveInputRecordName() const { return _IFT_QTruss1dGradDamage_Name; }
    virtual const char *giveClassName() const { return "QTruss1dGradDamage"; }

    virtual MaterialMode giveMaterialMode() { return _1dMat; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int computeNumberOfDofs() { return 5; }

protected:
    virtual void computeBdMatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeNdMatrixAt(GaussPoint *gp, FloatArray &answer);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
    virtual void computeGaussPoints();
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    void giveDofManDofIDMask_u(IntArray &answer) const;
    void giveDofManDofIDMask_d(IntArray &answer) const;
    
    virtual StructuralElement *giveStructuralElement() { return this; }
    virtual NLStructuralElement *giveNLStructuralElement() { return this; }
    virtual void giveLocationArray_u(IntArray &answer){;}
    virtual void giveLocationArray_d(IntArray &answer){;}
};
} // end namespace oofem
#endif // truss1d_h
