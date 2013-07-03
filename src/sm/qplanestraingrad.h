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

#ifndef qplanestraingrad_h
#define qplanestraingrad_h
#include "structuralelement.h"
#include "gaussintegrationrule.h"
#include "graddpelement.h"
#include "qplanestrain.h"
#include "fei2dquadlin.h"

#define _IFT_QPlaneStrainGrad_Name "qplanestraingrad"

namespace oofem {
class QPlaneStrainGrad : public QPlaneStrain, public GradDpElement
{
protected:
    int numberOfGaussPoints;
    static FEI2dQuadLin interpolation;

public:
    QPlaneStrainGrad(int n, Domain *d);
    virtual ~QPlaneStrainGrad() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveInputRecordName() const { return _IFT_QPlaneStrainGrad_Name; }
    virtual const char *giveClassName() const { return "QPlaneStrainGrad"; }
    virtual classType giveClassID() const { return QPlaneStrainGradClass; }

    virtual MaterialMode giveMaterialMode() { return _PlaneStrain; }
    virtual int computeNumberOfDofs(EquationID ut) { return 20; }

protected:
    virtual void computeBkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeNkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) { GradDpElement :: computeStiffnessMatrix(answer, rMode, tStep); }
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) { GradDpElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord); }
    virtual void computeForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode) { GradDpElement :: computeForceLoadVector(answer, stepN, mode); }
    virtual void computeGaussPoints();
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    virtual StructuralElement *giveStructuralElement() { return this; }
    virtual NLStructuralElement *giveNLStructuralElement() { return this; }
};
} // end namespace oofem
#endif // qplanestraingrad_h
