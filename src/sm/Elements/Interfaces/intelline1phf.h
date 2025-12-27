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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef intelline1phf_h
#define intelline1phf_h

#include "sm/Elements/Interfaces/structuralinterfaceelementphf.h"

#define _IFT_IntElLine1PhF_Name "intelline1phf"
#define _IFT_IntElLine1PhF_axisymmode "axisymmode"

namespace oofem {
class FEI2dLineLin;
class ParamKey;

/**
 * 
 * @author Jim Brouzoulis
 */
class IntElLine1PhF : public StructuralInterfaceElementPhF
{
protected:
    static FEI2dLineLin interp;
    /// Flag controlling axisymmetric mode (integration over unit circumferential angle)
    bool axisymmode = false;

    static ParamKey IPK_IntElLine1PhF_axisymmode;

public:
    IntElLine1PhF(int n, Domain * d);

    FEInterpolation *giveInterpolation() const override;

    //int computeNumberOfDofs() override { return 8; }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    double computeAreaAround(GaussPoint *gp) override;
    void computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    FloatArrayF<2> computeCovarBaseVectorAt(GaussPoint *gp) const;

    int testElementExtension(ElementExtension ext) override { return 0; }

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_IntElLine1PhF_Name; }
    const char *giveClassName() const override { return "IntElLine1PhF"; }
    void initializeFrom(InputRecord &ir, int priority) override;

    //void giveEngTraction(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep) override;
    void giveEngTraction(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, const double damage, TimeStep *tStep) override;

    void giveStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, IntegrationPoint *ip, TimeStep *tStep) override
    {
        answer = this->giveInterfaceCrossSection()->give2dStiffnessMatrix_Eng(rMode, ip, tStep);
    }

    void giveDofManDofIDMask_u(IntArray &answer) override;
    void giveDofManDofIDMask_d(IntArray &answer) override;
    void getLocationArray_u( IntArray &answer) override;
    void getLocationArray_d( IntArray &answer) override;
    void computeCovarBaseVectorsAt(GaussPoint *gp, FloatMatrix &G) override;

protected:
    void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    void computeGaussPoints() override;

    Element_Geometry_Type giveGeometryType() const override { return EGT_quad_1_interface; };
};
} // end namespace oofem
#endif
