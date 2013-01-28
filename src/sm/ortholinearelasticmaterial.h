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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef ortholinearelasticmaterial_h
#define ortholinearelasticmaterial_h

#include "linearelasticmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

#include "matconst.h"
#include "element.h"

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
    FloatMatrix *localCoordinateSystem;
    FloatArray *helpPlaneNormal;
    // in localCoordinateSystem the unity vectors are stored
    // COLUMWISE (this is exception, but allows faster numerical
    // implementation)

public:

    OrthotropicLinearElasticMaterial(int n, Domain *d) : LinearElasticMaterial(n, d)
    {
        localCoordinateSystem = NULL;
        helpPlaneNormal = NULL;
        cs_type = unknownCS;
    }
    virtual ~OrthotropicLinearElasticMaterial()
    {
        if ( localCoordinateSystem ) { delete localCoordinateSystem; }

        if ( helpPlaneNormal ) { delete helpPlaneNormal; }
    }

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "OrthotropicLinearElasticMaterial"; }
    virtual classType giveClassID() const { return OrthotropicLinearElasticMaterialClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual double give(int aProperty, GaussPoint *gp);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                       MatResponseForm form, MatResponseMode mode, GaussPoint * gp,
                                       TimeStep *tStep);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    /// Computes local 3d stiffness matrix of the receiver.
    virtual void give3dLocalMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseForm form, MatResponseMode mode, GaussPoint *gp,
                                            TimeStep *tStep);

protected:
    void giveTensorRotationMatrix(FloatMatrix &answer, GaussPoint *gp);
    void giveRotationMatrix(FloatMatrix &answer, GaussPoint *gp);

    friend class CrossSection;
};
} // end namespace oofem
#endif // ortholinearelasticmaterial_h
