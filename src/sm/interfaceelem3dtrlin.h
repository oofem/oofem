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

#ifndef interfaceelem3dtrlin_h
#define interfaceelem3dtrlin_h


#include "structuralelement.h"
#include "gaussintegrationrule.h"
#include "fei2dtrlin.h"

/*
 * This class implements 3d triangular surface interface element with linear interpolation.
 */
class InterfaceElement3dTrLin : public StructuralElement
{
protected:
    static FEI2dTrLin interpolation;

public:
    InterfaceElement3dTrLin(int, Domain *);                     // constructor
    ~InterfaceElement3dTrLin()   { }                            // destructor

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

    virtual int            computeNumberOfDofs(EquationID ut) { return 18; }
    virtual void           giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    double        computeVolumeAround(GaussPoint *);


    virtual int testElementExtension(ElementExtension ext) { return 0; }

    /** Interface requesting service */
    Interface *giveInterface(InterfaceType) { return NULL; }

#ifdef __OOFEG
    void          drawRawGeometry(oofegGraphicContext &);
    void          drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    void          drawScalar(oofegGraphicContext &context);
#endif
    //
    // definition & identification
    //
    const char *giveClassName() const { return "InterfaceElement3dTrLin"; }
    classType            giveClassID() const { return InterfaceElement3dTrLinClass; }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }

protected:
    void          computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    void          computeNmatrixAt(GaussPoint *, FloatMatrix &) { }
    void          computeGaussPoints();
    integrationDomain  giveIntegrationDomain() { return _Triangle; }

    int           giveApproxOrder() { return 1; }
    /*
     * internal
     * Note: this is not overloaded computeGtoLRotationMatrix
     * the gp parameter is needed, since element geometry can be (in general) curved
     *
     * //void          computeGtoLRotationMatrix(FloatMatrix& answer, GaussPoint* gp);
     */
    int          computeGtoLRotationMatrix(FloatMatrix &answer);
    void         computeLCS(FloatMatrix &answer);
};

#endif // interfaceelem3dtrlin_h








