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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#ifndef interfaceelem2dlin_h
#define interfaceelem2dlin_h

#include "../sm/Elements/structuralelement.h"

#define _IFT_InterfaceElem2dLin_Name "interface2dlin"
#define _IFT_InterfaceElem2dLin_axisymmode "axisymmode"

namespace oofem {
class FEI2dLineLin;

/**
 * This class implements a two dimensional interface element.
 * The approximatuion of geometry and unknowns is linear.
 */
class InterfaceElem2dLin : public StructuralElement
{
protected:
    static FEI2dLineLin interp;
    /// Flag controlling axisymmetric mode (integration over unit circumferential angle)
    bool axisymmode;

public:
    InterfaceElem2dLin(int n, Domain * d);
    virtual ~InterfaceElem2dLin() { }

    virtual FEInterpolation *giveInterpolation() const;

    virtual int computeNumberOfDofs() { return 8; }
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;

    virtual double computeVolumeAround(GaussPoint *gp);


    virtual int testElementExtension(ElementExtension ext) { return 0; }

    virtual Interface *giveInterface(InterfaceType) { return NULL; }

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
#endif

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_InterfaceElem2dLin_Name; }
    virtual const char *giveClassName() const { return "InterfaceElem2dLin"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual MaterialMode giveMaterialMode() { return _2dInterface; }

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) { }
    virtual void computeGaussPoints();

    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
};
} // end namespace oofem
#endif // interfaceelem2dlin_h
