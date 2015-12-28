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

#ifndef qtrspacegrad_h
#define qtrspacegrad_h

#include "../sm/Elements/3D/qtrspace.h"
#include "../sm/Elements/graddpelement.h"

#define _IFT_QTRSpaceGrad_Name "qtrspacegrad"

namespace oofem {
class FEI3dTetLin;

/**
 * Quadratic 3D element
 * @author L. Svoboda
 */
class QTRSpaceGrad : public QTRSpace, public GradDpElement
{
protected:
    static FEI3dTetLin interpolation_lin;
    
public:
    QTRSpaceGrad(int, Domain *);
    virtual ~QTRSpaceGrad() { }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_QTRSpaceGrad_Name; }
    virtual const char *giveClassName() const { return "QTRSpaceGrad"; }
    virtual int computeNumberOfDofs() { return 34; }

protected:
    ///////////////////////////////////////////////////////////////////////////////
    void computeGaussPoints();
    void computeNkappaMatrixAt(GaussPoint *gp, FloatArray &answer);
    void computeBkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer);
    StructuralElement *giveStructuralElement() { return this; }
    NLStructuralElement *giveNLStructuralElement() { return this; }

    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) { GradDpElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord); }
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) { GradDpElement :: computeStiffnessMatrix(answer, rMode, tStep); }
    virtual void computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode) { GradDpElement :: computeForceLoadVector(answer, tStep, mode); }
};
}
#endif // end namespace oofem
