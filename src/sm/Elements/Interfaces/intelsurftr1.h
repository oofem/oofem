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

#ifndef intelsurftr1_h
#define intelsurftr1_h

#include "sm/Elements/Interfaces/structuralinterfaceelement.h"
#include "fei3dtrlin.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"

#define _IFT_IntElSurfTr1_Name "intelsurftr1"

namespace oofem {
/**
 * This class implements 3d triangular surface interface element with linear interpolation.
 */
class IntElSurfTr1 : public StructuralInterfaceElement
{
protected:
    static FEI3dTrLin interpolation;

public:
    IntElSurfTr1(int n, Domain *d);

    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;
    bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords) override;
    virtual void computeCovarBaseVectorsAt(IntegrationPoint *ip, FloatArray &G1, FloatArray &G2);
    //bool computeGtoLRotationMatrix(FloatMatrix &answer) override;
    void computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer) override;

    int computeNumberOfDofs() override { return 18; }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    double computeAreaAround(IntegrationPoint *ip) override;

    void giveEngTraction(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep) override
    {
        answer = this->giveInterfaceCrossSection()->giveEngTraction_3d(jump, gp, tStep);
    }

    void giveStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, IntegrationPoint *ip, TimeStep *tStep) override
    {
        answer = this->giveInterfaceCrossSection()->give3dStiffnessMatrix_Eng(rMode, ip, tStep);
    }

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_IntElSurfTr1_Name; }
    void initializeFrom(InputRecord &ir) override;
    Element_Geometry_Type giveGeometryType() const override { return EGT_wedge_1; }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

protected:
    void computeNmatrixAt(GaussPoint *ip, FloatMatrix &answer) override;
    void computeGaussPoints() override;

};
} // end namespace oofem
#endif // intelsurftr1_h
