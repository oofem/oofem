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

#ifndef intelpoint_h
#define intelpoint_h

#include "../sm/Elements/Interfaces/structuralinterfaceelement.h"
#include "gaussintegrationrule.h"

///@name Input fields for Material
//@{
#define _IFT_IntElPoint_Name "intelpoint"
#define _IFT_IntElPoint_refnode "refnode"
#define _IFT_IntElPoint_normal "normal"
#define _IFT_IntElPoint_area "area"
//@}

namespace oofem {
/**
 * This class implements an interface element that connects two nodes.
 * In order to compute the normal and tangential direction of the slip plane, 
 * a reference node or specific direction must be specified in the input.
 * @author Jim Brouzoulis
 * @author Borek Patzak
 */
class IntElPoint : public StructuralInterfaceElement
{
protected:
    enum cmode { ie1d_1d, ie1d_2d, ie1d_3d } mode;
    int referenceNode;
    FloatArray normal;
    double area;
public:
    IntElPoint(int n, Domain *d);
    virtual ~IntElPoint() { }

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    virtual int computeNumberOfDofs();
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;

    virtual double computeAreaAround(GaussPoint *gp);



    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_IntElPoint_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }

    virtual MaterialMode giveMaterialMode();


    virtual void computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer);

    virtual void giveEngTraction(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
    {
        if ( this->giveDomain()->giveNumberOfSpatialDimensions() == 3 ) {
            this->giveInterfaceCrossSection()->giveEngTraction_3d(answer, gp, jump, tStep);
        } else if ( this->giveDomain()->giveNumberOfSpatialDimensions() == 2 ) {
            this->giveInterfaceCrossSection()->giveEngTraction_2d(answer, gp, jump, tStep);
        } else if ( this->giveDomain()->giveNumberOfSpatialDimensions() == 1 ) {
            this->giveInterfaceCrossSection()->giveEngTraction_1d(answer, gp, jump, tStep);
        }
    }

    virtual void giveStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, IntegrationPoint *ip, TimeStep *tStep)
    {
        if ( this->giveDomain()->giveNumberOfSpatialDimensions() == 3 ) {
            this->giveInterfaceCrossSection()->give3dStiffnessMatrix_Eng(answer, rMode, ip, tStep);
        } else if ( this->giveDomain()->giveNumberOfSpatialDimensions() == 2 ) {
            this->giveInterfaceCrossSection()->give2dStiffnessMatrix_Eng(answer, rMode, ip, tStep);
        } else if ( this->giveDomain()->giveNumberOfSpatialDimensions() == 1 ) {
            this->giveInterfaceCrossSection()->give1dStiffnessMatrix_Eng(answer, rMode, ip, tStep);
        }
        
    }

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
#endif    
    
protected:
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeGaussPoints();

    void computeLocalSlipDir(FloatArray &normal);
    cmode giveCoordMode() const { return this->mode; }
    void setCoordMode();
};
} // end namespace oofem
#endif // interfaceelement1d_h
