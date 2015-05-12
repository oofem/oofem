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

#ifndef intelline1_h
#define intelline1_h

#include "../sm/Elements/Interfaces/structuralinterfaceelement.h"

#define _IFT_IntElLine1_Name "intelline1"
#define _IFT_IntElLine1_axisymmode "axisymmode"

namespace oofem {
class FEI2dLineLin;

/**
 * This class implements a two dimensional interface element.
 * Even if geometry approx is quadratic, the element is assumed straight
 * If not straight, the rotation matrix depends on actual integration point
 * and stiffness and strain computations should be modified.
 * 
 * @author Jim Brouzoulis
 * @author Borek Patzak
 */
class IntElLine1 : public StructuralInterfaceElement
{
protected:
    static FEI2dLineLin interp;
    /// Flag controlling axisymmetric mode (integration over unit circumferential angle)
    bool axisymmode;
    
public:
    IntElLine1(int n, Domain * d);
    virtual ~IntElLine1() { }

    virtual FEInterpolation *giveInterpolation() const;

    virtual int computeNumberOfDofs() { return 8; }
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;

    virtual double computeAreaAround(GaussPoint *gp);
    virtual void computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeCovarBaseVectorAt(GaussPoint *gp, FloatArray &G);

    virtual int testElementExtension(ElementExtension ext) { return 0; }

    //virtual Interface *giveInterface(InterfaceType) { return NULL; }

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_IntElLine1_Name; }
    virtual const char *giveClassName() const { return "IntElLine1"; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveEngTraction(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
    {
        this->giveInterfaceCrossSection()->giveEngTraction_2d(answer, gp, jump, tStep);
    }

    virtual void giveStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, IntegrationPoint *ip, TimeStep *tStep)
    {
        this->giveInterfaceCrossSection()->give2dStiffnessMatrix_Eng(answer, rMode, ip, tStep);
    }

    #ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
    #endif
    
protected:
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeGaussPoints();

    Element_Geometry_Type giveGeometryType() const { return EGT_quad_1_interface; }
};
} // end namespace oofem
#endif
