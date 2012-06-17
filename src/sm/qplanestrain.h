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

#ifndef qplanstrain_h
#define qplanstrain_h

#include "structuralelement.h"
#include "fei2dquadquad.h"
#include "zznodalrecoverymodel.h"
#include "mathfem.h"

namespace oofem {
/**
 * This class implements an Quadratic isoparametric eight-node quadrilateral plane-
 * stress elasticity finite element. Each node has 2 degrees of freedom.
 */
class QPlaneStrain : public StructuralElement, public ZZNodalRecoveryModelInterface
{
protected:
    int numberOfGaussPoints;
    static FEI2dQuadQuad interpolation;

public:
    QPlaneStrain(int N, Domain *d);
    virtual ~QPlaneStrain() { }

    virtual FEInterpolation *giveInterpolation() { return &interpolation; }

    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    // definition & identification
    virtual const char *giveClassName() const { return "QPlaneStrain"; }
    virtual classType giveClassID() const { return QPlaneStrainClass; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_quad_2; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int computeNumberOfDofs(EquationID ut) { return 16; }

    virtual int testElementExtension(ElementExtension ext) { return 0; }

    virtual Interface *giveInterface(InterfaceType it);

    virtual double computeVolumeAround(GaussPoint *gp);

    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane) {
        return this->giveLenghtInDir(normalToCrackPlane) / sqrt( ( double ) this->numberOfGaussPoints );
    }

    virtual int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    void drawScalar(oofegGraphicContext &context);
    //void drawInternalState(DrawMode mode);
#endif

    virtual integrationDomain giveIntegrationDomain() { return _Square; }
    virtual MaterialMode giveMaterialMode() { return _PlaneStrain; }

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    void computeGaussPoints();
};
} // end namespace oofem
#endif // qplanstrain_h
