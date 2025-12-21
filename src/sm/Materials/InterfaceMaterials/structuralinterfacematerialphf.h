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
#ifndef structuralinterfacematerialphf_h
#define structuralinterfacematerialphf_h

#include "sm/Materials/InterfaceMaterials/structuralinterfacematerial.h"
#include "floatmatrixf.h"

namespace oofem {
class GaussPoint;

/**
 *
 * @author Jim Brouzoulis
 */
class StructuralInterfaceMaterialPhF : public StructuralInterfaceMaterial
{
    struct Tangents {
        FloatMatrixF<3,3> jj;
        FloatArrayF<3> jd;
        FloatArrayF<3> dj;
        double dd; 
    };

public:
    StructuralInterfaceMaterialPhF(int n, Domain * d);

    //virtual FloatArrayF<1> giveEngTraction_1d(const FloatArrayF<1> &jump, double damage, GaussPoint *gp, TimeStep *tStep);
    virtual FloatArrayF<2> giveEngTraction_2d(const FloatArrayF<2> &jump, double damage, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatArrayF<3> giveEngTraction_3d(const FloatArrayF<3> &jump, double damage, GaussPoint *gp, TimeStep *tStep) const
    { OOFEM_ERROR("not implemented "); };

    virtual double giveDrivingForce(GaussPoint *gp) const { OOFEM_ERROR("not implemented "); return -1; }
    virtual double giveDrivingForcePrime(GaussPoint *gp) const { OOFEM_ERROR("not implemented "); return -1; }

    // Df/Dd
    //virtual void giveDrivingForceTangent(FloatMatrix &answer,  MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveTangents(FloatMatrix &jj, FloatMatrix &jd, FloatMatrix &dj, FloatMatrix &dd, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
    { OOFEM_ERROR("not implemented "); };

    const char *giveClassName() const override { return "StructuralInterfaceMaterialPhF"; }
};
} // end namespace oofem
#endif // StructuralInterfaceMaterialPhF_h
