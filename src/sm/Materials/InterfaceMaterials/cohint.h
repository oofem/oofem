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

#ifndef cohint_h
#define cohint_h

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

namespace oofem {
///@name Input fields for cohesive interface material
//@{
#define _IFT_CohesiveInterfaceMaterial_Name "cohint"
#define _IFT_CohesiveInterfaceMaterial_kn "kn"
#define _IFT_CohesiveInterfaceMaterial_ks "ks"
//@}

/**
 * Class representing cohesive interface material model.
 * @author Milan Jirasek
 */
class CohesiveInterfaceMaterial : public StructuralInterfaceMaterial
{
protected:
    /// Elastic properties (normal and shear moduli).
    double kn, ks;

public:
    /// Constructor
    CohesiveInterfaceMaterial(int n, Domain * d);
    /// Destructor
    virtual ~CohesiveInterfaceMaterial() { }

    virtual int hasNonLinearBehaviour() { return 0; }
    virtual bool hasAnalyticalTangentStiffness() const { return true; }

    virtual const char *giveClassName() const { return "CohesiveInterfaceMaterial"; }
    virtual const char *giveInputRecordName() const { return _IFT_CohesiveInterfaceMaterial_Name; }

    virtual void giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep);
    virtual void give3dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new StructuralInterfaceMaterialStatus(1, FEMComponent :: domain, gp); }
};
} // namespace oofem
#endif
