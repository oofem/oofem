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

#ifndef interfaceelem3dtrlin_h
#define interfaceelem3dtrlin_h

#include "structuralelement.h"
#include "fei2dtrlin.h"

#define _IFT_InterfaceElement3dTrLin_Name "interface3dtrlin"

namespace oofem {
/**
 * This class implements 3d triangular surface interface element with linear interpolation.
 */
class InterfaceElement3dTrLin : public StructuralElement
{
protected:
    ///@todo Implement FEI3dTrLin, then remove giveIntegrationDomain and giveElementGeometry
    static FEI2dTrLin interpolation;

public:
    InterfaceElement3dTrLin(int n, Domain *d);
    virtual ~InterfaceElement3dTrLin() { }

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    virtual bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);

    virtual int computeNumberOfDofs() { return 18; }
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;

    virtual double computeVolumeAround(GaussPoint *gp);

    virtual int testElementExtension(ElementExtension ext) { return 0; }

    virtual Interface *giveInterface(InterfaceType) { return NULL; }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    void drawScalar(oofegGraphicContext &context);
#endif

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_InterfaceElement3dTrLin_Name; }
    virtual const char *giveClassName() const { return "InterfaceElement3dTrLin"; }
    virtual classType giveClassID() const { return InterfaceElement3dTrLinClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_wedge_1; }
    virtual integrationDomain giveIntegrationDomain() const { return _Triangle; }
    virtual MaterialMode giveMaterialMode() { return _3dInterface; }

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) { }
    virtual void computeGaussPoints();

    virtual int giveApproxOrder() { return 1; }

    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
    void computeLCS(FloatMatrix &answer);
};
} // end namespace oofem
#endif // interfaceelem3dtrlin_h
