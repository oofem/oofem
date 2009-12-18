/* $Header: /home/cvs/bp/oofem/tm/src/tr1_ht.h,v 1.2 2003/04/23 14:22:15 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   *******************************************************************************
//   *** CLASS Tr1_ht: Triangle(2d), linear approximation, Heat Transfer element ***
//   *******************************************************************************

#ifndef tr1_ht_h
#define tr1_ht_h

#include "transportelement.h"
#include "spatiallocalizer.h"

namespace oofem {

class Tr1_ht : public TransportElement, public SpatialLocalizerInterface
{
protected:
    int numberOfGaussPoints;
    double area;

public:

    // constructor
    Tr1_ht(int, Domain *, ElementMode em = HeatTransferEM);
    ~Tr1_ht();                       // destructor

    /** Computes the contribution to balance equation(s) due to internal sources */
    virtual void          computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *, ValueModeType mode);
    double                computeVolumeAround(GaussPoint *);
    /**
     * Computes the global coordinates from given element's local coordinates.
     * Required by nonlocal material models. Child classes should overload this function only
     * if they can be used together with nonlocal materil (where nonlocal averaging over
     * surronding volume is used).
     * @returns nonzero if successful; zero otherwise
     */
    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    /**
     * Computes the element local coordinates from given global coordinates.
     * @returns nonzero if successful (if point is inside element); zero otherwise
     */
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);
    // definition
    const char *giveClassName() const { return "Tr1_htElement"; }
    classType                giveClassID() const { return Tr1_htClass; }

    virtual int            computeNumberOfDofs(EquationID ut) { return 3; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    IRResultType           initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }

    /** Interface requesting service */
    Interface *giveInterface(InterfaceType);

    /**
     * @name The element interface required by SpatialLocalizerInterface
     */
    //@{
    /// Returns reference to corresponding element
    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    /// Returns nonzero if given element contains given point
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    /// Returns distance of given point from element parametric center
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);
    //@}

#ifdef __OOFEG
    //
    // Graphics output
    //
    //void          drawYourself (oofegGraphicContext&);
    //virtual void  drawRawGeometry (oofegGraphicContext&) {}
    //virtual void  drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

protected:
    void                  computeGaussPoints();

    virtual void  computeGradientMatrixAt(FloatMatrix &answer, GaussPoint *);
    virtual void  computeNmatrixAt(FloatMatrix &n, FloatArray *);
    /* computes the submatrix of interpolation matrix cooresponding to single unknown.
     */
    virtual void  computeNSubMatrixAt(FloatMatrix &n, FloatArray *);

    double  giveArea();

    void computeEgdeNMatrixAt(FloatMatrix &n, GaussPoint *gp);
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    void giveEdgeDofMapping(IntArray &mask, int iEdge);
    void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);

    int giveApproxOrder(int unknownIndx) { return 1; }
};

} // end namespace oofem
#endif // tr1_ht_h
