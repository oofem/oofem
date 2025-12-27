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

#ifndef ortholinearelasticmaterial_h
#define ortholinearelasticmaterial_h

#include "sm/Materials/linearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matconst.h"
#include "element.h"

#include <memory>

///@name Input fields for OrthotropicLinearElasticMaterial
//@{
#define _IFT_OrthotropicLinearElasticMaterial_Name "orthole"
#define _IFT_OrthotropicLinearElasticMaterial_ex "ex"
#define _IFT_OrthotropicLinearElasticMaterial_ey "ey"
#define _IFT_OrthotropicLinearElasticMaterial_ez "ez"
#define _IFT_OrthotropicLinearElasticMaterial_nyyz "nyyz"
#define _IFT_OrthotropicLinearElasticMaterial_nyxz "nyxz"
#define _IFT_OrthotropicLinearElasticMaterial_nyxy "nyxy"
#define _IFT_OrthotropicLinearElasticMaterial_gyz "gyz"
#define _IFT_OrthotropicLinearElasticMaterial_gxz "gxz"
#define _IFT_OrthotropicLinearElasticMaterial_gxy "gxy"
#define _IFT_OrthotropicLinearElasticMaterial_talphax "talphax"
#define _IFT_OrthotropicLinearElasticMaterial_talphay "talphay"
#define _IFT_OrthotropicLinearElasticMaterial_talphaz "talphaz"
#define _IFT_OrthotropicLinearElasticMaterial_lcs "lcs"
#define _IFT_OrthotropicLinearElasticMaterial_scs "scs"
//@}

namespace oofem {
class GaussPoint;

/// Coordinate system type.
enum CS_type {
    unknownCS, ///< Unknown coordinate system.
    localCS, ///< Coordinate system of principal axes is specified in global coordinate system (general).
    /**
     * coordinate system of principal axes is specified in shell  coordinate system
     * this is defined as follows: principal z-axis is perpendicular to mid-section *
     * x-axis is perpendicular to z-axis and normal to user specified vector n.
     * (so x-axis is parallel to plane, with n being normal to this plane).
     * y-axis is then perpendicular both to x and z axes.
     * @note This definition of cs is valid only for plates and shells
     * when vector n is parallel to z-axis an error occurs and program is terminated.
     */
    shellCS,
};

/**
 * This class implements a orthotropic linear elastic material in a finite
 * element problem.
 * For large deformations it becomes the orthotropic St. Venant-Kirchoff hyperelasticity model.
 *
 * Tasks:
 * - Returning standard material stiffness marix for 3d-case.
 *   according to current state determined by using data stored
 *   in Gausspoint, and local coordinate system defined in gp.
 * - Methods Give2dPlaneStressMtrx, GivePlaneStrainMtrx, Give1dStressMtrx are
 *   overloaded since form of this matrices is well known, and for
 *   faster response mainly in linear elastic problems.
 * - Returning a material property (method 'give'). Only for non-standard elements.
 * - Returning real stress state vector(tensor) at gauss point for 3d - case.
 * - Transforming stiffness from principal orhotrophy axes to
 *   system used in given GaussPoint.
 */
class OrthotropicLinearElasticMaterial : public LinearElasticMaterial
{
protected:
    CS_type cs_type;
    FloatMatrixF<3,3> localCoordinateSystem;
    FloatArrayF<3> helpPlaneNormal;
    // in localCoordinateSystem the unity vectors are stored
    // COLUMWISE (this is exception, but allows faster numerical
    // implementation)

public:

    OrthotropicLinearElasticMaterial(int n, Domain * d) : LinearElasticMaterial(n, d),
        cs_type(unknownCS)
    { }

    FloatArrayF<6> giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const override;

    const char *giveInputRecordName() const override { return _IFT_OrthotropicLinearElasticMaterial_Name; }
    const char *giveClassName() const override { return "OrthotropicLinearElasticMaterial"; }
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    double give(int aProperty, GaussPoint *gp) const override;

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

protected:
    FloatMatrixF<3,3> giveTensorRotationMatrix(GaussPoint *gp) const;
    FloatMatrixF<6,6> giveRotationMatrix(GaussPoint *gp) const;

    friend class CrossSection;
};
} // end namespace oofem
#endif // ortholinearelasticmaterial_h
