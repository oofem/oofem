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

#ifndef interfaceelement1d_h
#define interfaceelement1d_h

#include "structuralelement.h"
#include "gaussintegrationrule.h"

namespace oofem {
/**
 * This class implements a one-dimensional interface element connecting two nodes (with the same position)
 * In order to compute normal and tangential direction of the slip plane, a reference node or specific direction is needed.
 */
class InterfaceElem1d : public StructuralElement
{

protected:
    enum cmode { ie1d_1d, ie1d_2d, ie1d_3d } mode;
    int referenceNode;
    FloatArray normal;

public:
    InterfaceElem1d(int n, Domain *d);
    ~InterfaceElem1d() { }

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)  { computeLumpedMassMatrix(answer, tStep); }

    /**
     * Computes the global coordinates from given element's local coordinates.
     * Required by nonlocal material models.
     * @returns nonzero if successful
     *
     */
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    /**
     * Computes the element local coordinates from given global coordinates.
     * @returns nonzero if successful (if point is inside element); zero otherwise
     */
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);

    virtual int computeNumberOfDofs(EquationID ut);
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;

    double computeVolumeAround(GaussPoint *gp);


    virtual int testElementExtension(ElementExtension ext) { return 0; }

    /** Interface requesting service */
    Interface *giveInterface(InterfaceType it) { return NULL; }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    void drawScalar(oofegGraphicContext &context);
#endif

    // definition & identification
    const char *giveClassName() const { return "InterfaceElem1d"; }
    classType giveClassID() const { return InterfaceElem1dClass; }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_point; }

    integrationDomain giveIntegrationDomain() { return _Point; }
    MaterialMode giveMaterialMode();

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer) { }
    void computeGaussPoints();
    void          evaluateLocalCoordinateSystem(FloatMatrix &);

    int giveApproxOrder() { return 1; }

    void computeLocalSlipDir(FloatArray &normal);
    cmode giveCoordMode() const { return this->mode; }
    void setCoordMode();
};
} // end namespace oofem
#endif // interfaceelement1d_h
