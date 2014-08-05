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

#ifndef q4axisymm_h
#define q4axisymm_h

#include "structuralelement.h"
#include "fei2dquadquad.h"
#include "zznodalrecoverymodel.h"

///@name Input fields for Q4Axisymm
//@{
#define _IFT_Q4Axisymm_Name "q4axisymm"
#define _IFT_Q4Axisymm_nipfish "nipfish"
//@}

namespace oofem {
/**
 * This class implements an Quadratic isoparametric eight-node quadrilateral -
 * elasticity finite element for axisymmetric 3d continuum.
 * Each node has 2 degrees of freedom.
 * @todo Use FEI classes.
 */
class Q4Axisymm : public StructuralElement, public ZZNodalRecoveryModelInterface
{
protected:
    static FEI2dQuadQuad interp;
    int numberOfGaussPoints, numberOfFiAndShGaussPoints;

public:
    Q4Axisymm(int n, Domain * d);
    virtual ~Q4Axisymm();

    virtual FEInterpolation *giveInterpolation() const { return & interp; }

    virtual int computeNumberOfDofs() { return 16; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;
    virtual double giveCharacteristicLength(const FloatArray &normalToCrackPlane);
    virtual double computeVolumeAround(GaussPoint *gp);
    virtual void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    // definition & identification
    virtual Interface *giveInterface(InterfaceType);
    virtual const char *giveInputRecordName() const { return _IFT_Q4Axisymm_Name; }
    virtual const char *giveClassName() const { return "Q4axisymm"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialMode giveMaterialMode() { return _3dMat; }

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer);

    virtual void computeGaussPoints();

    FloatArray *GiveDerivativeKsi(double, double);
    FloatArray *GiveDerivativeEta(double, double);
    void computeJacobianMatrixAt(FloatMatrix &answer, GaussPoint *gp);

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc);
    void drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type);
#endif
};
} // end namespace oofem
#endif // q4axisymm_h
