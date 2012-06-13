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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef qtruss1dgrad_h
#define qtruss1dgrad_h

#include "fei1dlin.h"
#include "structuralelement.h"
#include "gaussintegrationrule.h"
#include "graddpelement.h"
#include "qtruss1d.h"

namespace oofem {

/**
 * This class implements a three-node gradient truss bar element for one-dimensional
 * analysis.
 */
class QTruss1dGrad : public QTruss1d, public GradDpElement
{
protected:
    double length;
    static FEI1dLin interpolation;
    int numberOfGaussPoints;

public:
    QTruss1dGrad(int n, Domain *d);
    virtual ~QTruss1dGrad() { }

    virtual const char *giveClassName() const { return "QTruss1dGrad"; }
    virtual classType giveClassID() const { return QTruss1dGradClass; }

    virtual Element_Geometry_Type giveGeometryType() const { return EGT_line_2; }
    virtual integrationDomain  giveIntegrationDomain() { return _Line; }
    virtual MaterialMode giveMaterialMode() { return _1dMatGrad; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    int getNprimNodes() { return 3; }
    int getNprimVars() { return 1; }
    int getNsecNodes() { return 2; }
    int getNsecVars() { return 1; }

protected:
    virtual void computeBkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeNkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
    virtual void computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    virtual void computeGaussPoints();
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    virtual StructuralElement* giveStructuralElement(){return this;}
};
} // end namespace oofem
#endif // truss1d_h
