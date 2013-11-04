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

#ifndef intelline2_h
#define intelline2_h

#include "structuralinterfaceelement.h"

#define _IFT_IntElLine2_Name "intelline2"

namespace oofem {

class FEI2dLineQuad;

/**
 * This class implements a two dimensional interface element.
 * Even if geometry approx is quadratic, the element is assumed straight
 * If not straight, the rotation matrix depends on actual integration point
 * and stiffness and strain computations should be modified.
 */
class IntElLine2 : public StructuralInterfaceElement
{
protected:
    static FEI2dLineQuad interp;

public:
    IntElLine2(int n, Domain *d);
    virtual ~IntElLine2() { }

    virtual FEInterpolation *giveInterpolation() const;

    virtual int computeNumberOfDofs() { return 12; }
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;

    virtual double computeAreaAround(GaussPoint *gp);
    virtual void computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeCovarBaseVectorAt(GaussPoint *gp, FloatArray &G);

    virtual int testElementExtension(ElementExtension ext) { return 0; }

    virtual Interface *giveInterface(InterfaceType) { return NULL; }



    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_IntElLine2_Name; }
    virtual const char *giveClassName() const { return "IntElLine2"; }
    virtual classType giveClassID() const { return IntElLine2Class; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    //virtual MaterialMode giveMaterialMode() { return _2dInterface; }

protected:
    //virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeGaussPoints();

    virtual int giveApproxOrder() { return 1; }
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
};
} // end namespace oofem
#endif 
