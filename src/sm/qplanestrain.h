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

#ifndef qplanstrain_h
#define qplanstrain_h

#include "nlstructuralelement.h"
#include "fei2dquadquad.h"
#include "zznodalrecoverymodel.h"
#include "mathfem.h"

#define _IFT_QPlaneStrain_Name "qplanestrain"

namespace oofem {
/**
 * This class implements an Quadratic isoparametric eight-node quadrilateral plane-
 * stress elasticity finite element. Each node has 2 degrees of freedom.
 */
class QPlaneStrain : public NLStructuralElement, public ZZNodalRecoveryModelInterface
{
protected:
    static FEI2dQuadQuad interpolation;

public:
    QPlaneStrain(int N, Domain *d);
    virtual ~QPlaneStrain() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual FEInterpolation *giveInterpolation() const { return & interpolation; }

    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_QPlaneStrain_Name; }
    virtual const char *giveClassName() const { return "QPlaneStrain"; }
    virtual classType giveClassID() const { return QPlaneStrainClass; }
    virtual int computeNumberOfDofs() { return 16; }
    virtual MaterialMode giveMaterialMode() { return _PlaneStrain; }

    virtual int testElementExtension(ElementExtension ext) { return 0; }

    virtual Interface *giveInterface(InterfaceType it);

    virtual double computeVolumeAround(GaussPoint *gp);

    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane);

    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &);
    virtual void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    virtual void drawScalar(oofegGraphicContext &context);
#endif

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);
    void computeGaussPoints();
};
} // end namespace oofem
#endif // qplanstrain_h
