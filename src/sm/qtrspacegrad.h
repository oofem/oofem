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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef qtrspacegrad_h
#define qtrspacegrad_h

#include "qtrspace.h"
#include "graddpelement.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "eleminterpmapperinterface.h"
#include "huertaerrorestimator.h"
#include "sprnodalrecoverymodel.h"

namespace oofem {

class FEI3dTetLin;

/**
 * Quadratic 3D element
 * @author L. Svoboda
 */
class QTRSpaceGrad : public QTRSpace,public GradDpElement
{
protected:
    int numberOfGaussPoints;
    static FEI3dTetLin interpolation;

public:
    QTRSpaceGrad(int,Domain*);
    virtual ~QTRSpaceGrad() {}

    virtual IRResultType initializeFrom(InputRecord* ir);
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray& answer) const;

    // definition & identification
    virtual const char* giveClassName() const { return "QTRSpaceGrad"; }
    virtual classType giveClassID() const { return QTRSpaceClass; }
    virtual int computeNumberOfDofs(EquationID ut) { return 34; }
  
protected:
    ///////////////////////////////////////////////////////////////////////////////
    void computeGaussPoints();
    void computeNkappaMatrixAt(GaussPoint*, FloatMatrix&);
    void computeBkappaMatrixAt(GaussPoint*, FloatMatrix&);
    StructuralElement* giveStructuralElement() { return this; }
    NLStructuralElement* giveNLStructuralElement() { return this; }

    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) { GradDpElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord); }
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) { GradDpElement :: computeStiffnessMatrix(answer, rMode, tStep); }
    virtual void computeForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode) { GradDpElement :: computeForceLoadVector(answer, stepN, mode); }
    //virtual void computeNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode) { GradDpElement :: computeNonForceLoadVector(answer, stepN, mode); }
    virtual void computeNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode) { ; }

};
 
}
#endif // end namespace oofem

