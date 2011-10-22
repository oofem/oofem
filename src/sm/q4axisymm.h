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

#ifndef q4axisymm_h
#define q4axisymm_h

#include "structuralelement.h"

namespace oofem {

/**
 * This class implements an Quadratic isoparametric eight-node quadrilateral -
 * elasticity finite element for axisimmetric 3d continuum.
 * Each node has 2 degrees of freedom.
 * @todo Use FEI classes.
 */
class Q4Axisymm : public StructuralElement
{
protected:
    int numberOfGaussPoints, numberOfFiAndShGaussPoints;

public:
    Q4Axisymm(int n, Domain *d);
    ~Q4Axisymm();

    virtual int computeNumberOfDofs(EquationID ut) { return 16; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    double computeVolumeAround(GaussPoint *gp);
    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    // definition & identification
    Interface *giveInterface(InterfaceType) { return NULL; }
    const char *giveClassName() const { return "Q4axisymm"; }
    classType giveClassID() const { return Q4AxisymmClass; }
    IRResultType initializeFrom(InputRecord *ir);
    integrationDomain giveIntegrationDomain() { return _Square; }
    MaterialMode giveMaterialMode() { return _3dMat; }

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    void computeGaussPoints();

    FloatArray *GiveDerivativeKsi(double, double);
    FloatArray *GiveDerivativeEta(double, double);
    void computeJacobianMatrixAt(FloatMatrix &answer, GaussPoint *gp);
};
} // end namespace oofem
#endif // q4axisymm_h

