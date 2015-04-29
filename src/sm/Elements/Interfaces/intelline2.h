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

#include "../sm/Elements/Interfaces/intelline1.h"

#define _IFT_IntElLine2_Name "intelline2"
#define _IFT_IntElLine2_LinearTraction "linear"

namespace oofem {
class FEI2dLineQuad;

/**
 * This class implements a two dimensional interface element and is simply an extension 
 * of IntElLine1 to a quadratic approximation.
 * @author Jim Brouzoulis
 * @author Borek Patzak
 */
class IntElLine2 : public IntElLine1
{
protected:
    static FEI2dLineQuad interp;
    static FEI2dLineLin interpLin;

public:
    IntElLine2(int n, Domain * d);
    virtual ~IntElLine2() { }
    virtual FEInterpolation *giveInterpolation() const;
    virtual int computeNumberOfDofs() { return 12; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    
    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_IntElLine2_Name; }
    virtual const char *giveClassName() const { return "IntElLine2"; }

    #ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
    #endif    
    
protected:
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeGaussPoints();

    Element_Geometry_Type giveGeometryType() const { return EGT_quad_21_interface; }

    /// If linear interpolation should be used.
	bool linear;

};
} // end namespace oofem
#endif
