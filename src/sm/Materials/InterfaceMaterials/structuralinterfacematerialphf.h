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
#ifndef structuralinterfacematerialphf_h
#define structuralinterfacematerialphf_h

#include "Materials/InterfaceMaterials/structuralinterfacematerial.h"

namespace oofem {
class GaussPoint;

/**
 *
 * @author Jim Brouzoulis
 */
class StructuralInterfaceMaterialPhF : public StructuralInterfaceMaterial
{
public:

    StructuralInterfaceMaterialPhF(int n, Domain * d);
    /// Destructor.
    virtual ~StructuralInterfaceMaterialPhF() { }


    //virtual void giveEngTraction_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, const double damage, TimeStep *tStep);
    virtual void giveEngTraction_2d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, const double damage, TimeStep *tStep);
    virtual void giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, const double damage, TimeStep *tStep)
    { OOFEM_ERROR("not implemented "); };
    

    virtual double giveDrivingForce(GaussPoint *gp) { OOFEM_ERROR("not implemented "); return -1;};
    virtual double giveDrivingForcePrime(GaussPoint *gp) { OOFEM_ERROR("not implemented "); return -1; };
    
    // Df/Dd
    //virtual void giveDrivingForceTangent(FloatMatrix &answer,  MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveTangents(FloatMatrix &jj, FloatMatrix &jd, FloatMatrix &dj, FloatMatrix &dd, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
    { OOFEM_ERROR("not implemented "); };
    
    
    
    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "StructuralInterfaceMaterialPhF"; }
    
};
} // end namespace oofem
#endif // StructuralInterfaceMaterialPhF_h
