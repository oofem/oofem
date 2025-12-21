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

#ifndef cohint_h
#define cohint_h

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

namespace oofem {
///@name Input fields for cohesive interface material
//@{
#define _IFT_CohesiveInterfaceMaterial_Name "cohint"
#define _IFT_CohesiveInterfaceMaterial_kn "kn"
#define _IFT_CohesiveInterfaceMaterial_transitionopening "transitionopening"
#define _IFT_CohesiveInterfaceMaterial_ks "ks"
#define _IFT_CohesiveInterfaceMaterial_stiffCoeffKn "stiffcoeffkn"
#define _IFT_CohesiveInterfaceMaterial_smoothMag "smoothmag"
//@}

/**
 * Class representing cohesive interface material model.
 * @author Milan Jirasek
 */
class CohesiveInterfaceMaterial : public StructuralInterfaceMaterial
{
protected:
    /// Elastic properties (normal and shear moduli).
    double kn = 0., ks = 0.;

    /// Reduction of normal stiffness when in tension
    double stiffCoeffKn = 0.;

    /// Opening when material stiffness changes from kn to kn*stiffCoeffKn
    double transitionOpening = 0.;

    /// Smoothing region between tension and compression stiffness.
    /// Uses atan(smoothMag*x) function, where (smoothMag*x) needs to be 6
    /// to be 10% off the asymptotic value. SmoothMag is by default 1.e+4.
    /// Higher values mean sharper transition.
    double smoothMag = 0.;

public:
    /// Constructor
    CohesiveInterfaceMaterial(int n, Domain * d);

    bool hasAnalyticalTangentStiffness() const override { return true; }

    const char *giveClassName() const override { return "CohesiveInterfaceMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_CohesiveInterfaceMaterial_Name; }

    FloatArrayF<3> giveEngTraction_3d(const FloatArrayF<3> &jump, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> give3dStiffnessMatrix_Eng(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<StructuralInterfaceMaterialStatus>(gp); }
};
} // namespace oofem
#endif
