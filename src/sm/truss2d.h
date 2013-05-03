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

#ifndef truss2d_h
#define truss2d_h

#include "nlstructuralelement.h"

///@name Input fields for 2D truss element
//@{
#define _IFT_Truss2d_Name "truss2d"
#define _IFT_Truss2d_cs "cs"
//@}

namespace oofem {

/**
 * This class implements a two-node truss bar element for two-dimensional
 * analysis.
 *
 * A truss bar element is characterized by its 'length' and its 'pitch'. The
 * pitch is the angle in radians between the X-axis and the axis of the
 * element (oriented node1 to node2).
 *
 * Can be used in xy, xz or yz planes.
 *
 * @author Peter Grassl
 */
class Truss2d : public NLStructuralElement
{
protected:
    double length;
    double pitch;
    int cs_mode;
    ///@todo Use interpolator class for lines in 2D

public:
    Truss2d(int n, Domain *d);
    virtual ~Truss2d() { }

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }
    virtual int giveLocalCoordinateSystem(FloatMatrix &answer);

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    virtual int computeNumberOfDofs(EquationID ut) { return 4; }
    virtual void giveDofManDofIDMask(int inode, EquationID eid, IntArray &answer) const;

    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &)
    { return this->giveLength(); }

    virtual double computeVolumeAround(GaussPoint *gp);

    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
#endif

    // definition & identification
    virtual const char *giveClassName() const { return "Truss2d"; }
    virtual classType giveClassID() const { return Truss2dClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }
    virtual integrationDomain giveIntegrationDomain() { return _Line; }
    virtual MaterialMode giveMaterialMode() { return _1dMat; }

protected:
    // edge load support
    void resolveCoordIndices(int &c1, int &c2);
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    virtual void giveEdgeDofMapping(IntArray &answer, int) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &, int, GaussPoint *gp);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &, int = 1, int = ALL_STRAINS);
    virtual void computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *gp, int);
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &);
    virtual void computeGaussPoints();

    double giveLength();
    double givePitch();
    virtual int giveApproxOrder() { return 1; }
};
} // end namespace oofem
#endif // truss2d_h
