/* $Header: /home/cvs/bp/oofem/sm/src/truss2d.h,v 1.6 2003/04/06 14:08:32 bp Exp $ */
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

//   *********************
//   *** CLASS TRUSS2D ***
//   *********************

/*
 * Contributions:
 * Peter Grassl - has added cs_mode, allowing the put this element in xy, xz, or yz planes
 */

#ifndef truss2d_h
#define truss2d_h


#include "nlstructuralelement.h"
#include "gaussintegrationrule.h"
#include "gausspnt.h"

namespace oofem {
/**
 * This class implements a two-node truss bar element for two-dimensional
 * analysis.
 * DESCRIPTION :
 * A truss bar element is characterized by its 'length' and its 'pitch'. The
 * pitch is the angle in radians between the X-axis anf the axis of the
 * element (oriented node1 to node2).
 */
class Truss2d : public NLStructuralElement
{
protected:
    double length;
    double pitch;
    // FloatMatrix*  rotationMatrix ;
    int cs_mode;
public:

    Truss2d(int, Domain *);                     // constructor
    ~Truss2d()   { }                            // destructor

    // FloatArray*   ComputeBodyLoadVectorAt (TimeStep*) ;
    void          computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    void          computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }
    // FloatArray*   ComputeResultingBodyForceAt (TimeStep*) ;
    // FloatMatrix*  computeStiffnessMatrix () ;
    // FloatArray*   ComputeStrainVector (GaussPoint*,TimeStep*) ;
    // FloatArray*   ComputeInitialStrainVector (TimeStep* );
    int              giveLocalCoordinateSystem(FloatMatrix &answer);
    /**
     * Computes the global coordinates from given element's local coordinates.
     * Required by nonlocal material models.
     * @returns nonzero if successful
     */
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    virtual int            computeNumberOfDofs(EquationID ut) { return 4; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    // characteristic length in gp (for some material models)
    double        giveCharacteristicLenght(GaussPoint *, const FloatArray &)
    { return this->giveLength(); }

    double        computeVolumeAround(GaussPoint *);

    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }
    //int    hasEdgeLoadSupport () {return 1;}

#ifdef __OOFEG
    void          drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
#endif
    //
    // definition & identification
    //
    const char *giveClassName() const { return "Truss2d"; }
    classType            giveClassID() const { return Truss2dClass; }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }
    integrationDomain  giveIntegrationDomain() { return _Line; }
    MaterialMode          giveMaterialMode()  { return _1dMat; }

protected:
    // edge load support
    void resolveCoordIndices(int &c1, int &c2);
    void  computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *);
    void  giveEdgeDofMapping(IntArray &answer, int) const;
    double        computeEdgeVolumeAround(GaussPoint *, int);
    void          computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
    { computeGlobalCoordinates( answer, * ( gp->giveCoordinates() ) ); }
    int   computeLoadLEToLRotationMatrix(FloatMatrix &, int, GaussPoint *);
    void          computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    void          computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *, int);
    void          computeNmatrixAt(GaussPoint *, FloatMatrix &);
    //  int           computeGtoNRotationMatrix (FloatMatrix&);
    void          computeGaussPoints();

    double        giveLength();
    double        givePitch();
    int           giveApproxOrder() { return 1; }
};
} // end namespace oofem
#endif // truss2d_h
