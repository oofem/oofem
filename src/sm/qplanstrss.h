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

#ifndef qplanstrss_h
#define qplanstrss_h

#include "structuralelement.h"
#include "fei2dquadquad.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "mathfem.h"

namespace oofem {
/**
 * This class implements an Quadratic isoparametric eight-node quadrilateral plane-
 * stress elasticity finite element. Each node has 2 degrees of freedom.
 */
class QPlaneStress2d : public StructuralElement, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface
{
protected:
    int numberOfGaussPoints;
    static FEI2dQuadQuad interpolation;

public:
    QPlaneStress2d(int n, Domain *d);
    ~QPlaneStress2d() { }

    virtual int computeNumberOfDofs(EquationID ut) { return 16; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    // definition & identification
    const char *giveClassName() const { return "QPlaneStress2d"; }
    classType giveClassID() const { return QPlaneStress2dClass; }
    Element_Geometry_Type giveGeometryType() const { return EGT_quad_2; }
    FEInterpolation *giveInterpolation() { return & interpolation; }
    IRResultType initializeFrom(InputRecord *ir);

    virtual int testElementExtension(ElementExtension ext) { return 0; }

    Interface *giveInterface(InterfaceType it);

    double computeVolumeAround(GaussPoint *gp);

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane);

    Element *ZZNodalRecoveryMI_giveElement() { return this; }
    int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint, InternalStateType type);
    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);
    void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    void drawScalar(oofegGraphicContext &context);
    //void drawInternalState(DrawMode mode);
#endif
    integrationDomain giveIntegrationDomain() { return _Square; }
    MaterialMode giveMaterialMode() { return _PlaneStress; }

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    void computeGaussPoints();
};
} // end namespace oofem
#endif // qplanstrss_h
