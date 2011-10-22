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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef qtrplstrgrad_h
#define qtrplstrgrad_h
#include "structuralelement.h"
#include "gaussintegrationrule.h"
#include "graddpelement.h"
#include "qtrplstr.h"
#include "fei2dquadlin.h"

namespace oofem {
class QTrPlaneStressGrad : public QTrPlaneStress2d, public GradDpElement
{
protected:
    int numberOfGaussPoints;
    static FEI2dQuadLin interpolation;

public:
    QTrPlaneStressGrad(int n, Domain *d);
    ~QTrPlaneStressGrad() { }

    IRResultType initializeFrom(InputRecord *ir);

    const char *giveClassName() const { return "QTrPlaneStressGrad"; }
    classType giveClassID() const { return QTrPlaneStressGradClass; }

    Element_Geometry_Type giveGeometryType() const { return EGT_quad_2; }
    integrationDomain giveIntegrationDomain() { return _Triangle; }
    MaterialMode giveMaterialMode() { return _PlaneStressGrad; }
    virtual int computeNumberOfDofs(EquationID ut) { return 15; }

protected:
    virtual void computeBkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeNkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) { GradDpElement :: computeStiffnessMatrix(answer, rMode, tStep); }
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) { GradDpElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord); }
    virtual void computeForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode) { GradDpElement :: computeForceLoadVector(answer, stepN, mode); }
    virtual void computeNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode) { GradDpElement :: computeNonForceLoadVector(answer, stepN, mode); }
    void computeGaussPoints();
    void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    StructuralElement *giveStructuralElement() { return this; }
};
} // end namespace oofem
#endif // qtrplstrgrad_h
