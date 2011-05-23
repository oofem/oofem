/* $Header: /home/cvs/bp/oofem/sm/src/qplanstrss.h,v 1.4.4.1 2004/04/05 15:19:47 bp Exp $ */
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

//   ************************************
//   *** CLASS QUADRATIC PLANE STRAIN ***
//   ************************************

#ifndef qplanstrss_h
#define qplanstrss_h

#include "structuralelement.h"
#include "fei2dquadquad.h"
#include "zznodalrecoverymodel.h"
#include "mathfem.h"

namespace oofem {
class QPlaneStress2d : public StructuralElement, public ZZNodalRecoveryModelInterface
{
    /*
     * This class implements an Quadratic isoparametric eight-node quadrilateral plane-
     * stress elasticity finite element. Each node has 2 degrees of freedom.
     *
     * DESCRIPTION :
     *
     * One single additional attribute is needed for Gauss integration purpose :
     * 'jacobianMatrix'. This 2x2 matrix contains polynomials.
     *
     * TASKS :
     *
     * - calculating its Gauss points ;
     * - calculating its B,D,N matrices and dV.
     */

protected:
    int numberOfGaussPoints;
    static FEI2dQuadQuad interpolation;
public:
    QPlaneStress2d(int, Domain *);                       // constructor
    ~QPlaneStress2d()  { }                               // destructor

    // characteristic length in gp (for some material models)
    // double        giveCharacteristicLenght (GaussPoint*, const FloatArray&) {return 0.;}
    virtual int            computeNumberOfDofs(EquationID ut) { return 16; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    //
    // definition & identification
    //
    const char *giveClassName() const { return "QPlaneStress2d"; }
    classType        giveClassID()   const { return QPlaneStress2dClass; }
    Element_Geometry_Type giveGeometryType() const { return EGT_quad_2; }
    FEInterpolation *giveInterpolation() { return & interpolation; } 
    IRResultType initializeFrom(InputRecord *ir);

    virtual int testElementExtension(ElementExtension ext) { return 0; }
    /** Interface requesting service */
    Interface *giveInterface(InterfaceType);
    //int    hasEdgeLoadSupport () {return 0;}
    double                computeVolumeAround(GaussPoint *);
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    /**
     * Returns characteristic length of element in given integration point and in
     * given direction. Required by material models relying on crack-band approach to achieve
     * objectivity with respect to mesh size
     */
    double        giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane) {
        return this->giveLenghtInDir(normalToCrackPlane) / sqrt( ( double ) this->numberOfGaussPoints );
    }

    /**
     * @name The element interface required by ZZNodalRecoveryModel
     */
    //@{
    Element *ZZNodalRecoveryMI_giveElement() { return this; }
    //@}

#ifdef __OOFEG
    void          drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    void drawScalar(oofegGraphicContext &context);
    //      void          drawInternalState (DrawMode mode);
#endif
    integrationDomain  giveIntegrationDomain() { return _Square; }
    MaterialMode          giveMaterialMode()  { return _PlaneStress; }

protected:
    void             computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    void             computeNmatrixAt(GaussPoint *, FloatMatrix &);
    void             computeGaussPoints();
};
} // end namespace oofem
#endif // qplanstrss_h
