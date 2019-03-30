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

#ifndef qwedgegraddamage_h
#define qwedgegraddamage_h

#include "../sm/Elements/3D/qwedge.h"
#include "../sm/Elements/GradientDamage/graddamageelement.h"

#define _IFT_QWedgeGradDamage_Name "qwedgegraddamage"

namespace oofem {
class FEI3dWedgeLin;

/**
 * @author M. Horak
 */
class QWedgeGradDamage : public QWedge, public GradientDamageElement
{
protected:
    static FEI3dWedgeLin interpolation_lin;

public:
    QWedgeGradDamage(int, Domain *);
    virtual ~QWedgeGradDamage() { }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    void giveDofManDofIDMask_u(IntArray &answer) const;
    void giveDofManDofIDMask_d(IntArray &answer) const;
    
    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_QWedgeGradDamage_Name; }
    virtual const char *giveClassName() const { return "QWedgeGradDamage"; }
    virtual int computeNumberOfDofs() { return 51; }
    virtual MaterialMode giveMaterialMode() { return _3dMat; }

protected:
    virtual void computeGaussPoints();
    virtual void computeNdMatrixAt(GaussPoint *gp, FloatArray &answer);
    virtual void computeBdMatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual StructuralElement *giveStructuralElement() { return this; }
    virtual NLStructuralElement *giveNLStructuralElement() { return this; }

    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) { GradientDamageElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord); }
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) { GradientDamageElement :: computeStiffnessMatrix(answer, rMode, tStep); }
    virtual void giveLocationArray_u(IntArray &answer){;}
    virtual void giveLocationArray_d(IntArray &answer){;}

};
}
#endif // end namespace oofem
