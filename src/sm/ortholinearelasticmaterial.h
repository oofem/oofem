/* $Header: /home/cvs/bp/oofem/sm/src/ortholinearelasticmaterial.h,v 1.5 2003/04/06 14:08:31 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   **************************************************
//   *** CLASS ORTHOTROPIC LINEAR ELACSTIC MATERIAL ***
//   **************************************************

#ifndef ortholinearelasticmaterial_h
#define ortholinearelasticmaterial_h

#include "linearelasticmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

#include "matconst.h"
#include "element.h"

class GaussPoint;


enum CS_type { unknownCS, localCS, shellCS };
//
// localCS  - coordinate system of principal axes is specified in global coordinate system (general)
// shellCS  - coordinate system of principal axes is specified in shell  coordinate system
//            this is defined as follows: principal z-axis is perpendicular to mid-section
//            x-axis is perpendicular to z-axis and normal to user specified vector n.
//            (so x-axis is parallel to plane, with n beeing normal to this plane).
//            y-axis is then perpendicular both to x and z axes.
//            WARNING: this definition of cs is valid only for plates and shells
//            when vector n is paralel to z-axis an error occurs and program is terminated.
//

class OrthotropicLinearElasticMaterial : public LinearElasticMaterial
{
    /*
     * This class implements a orthotropic linear elastic material in a finite
     * element problem. A material
     * is an attribute of a domain. It is usually also attribute of many elements.
     *
     * DESCRIPTION
     * ORTHOTROPIC Linear Elastic Material
     *
     * TASK
     * - Returning standard material stiffness marix for 3d-case.
     * according to current state determined by using data stored
     * in Gausspoint, and local coordinate system defined in gp.
     * - methods Give2dPlaneStressMtrx, GivePlaneStrainMtrx, Give1dStressMtrx are
     * overloaded since form of this matrices is well known, and for
     * faster response mainly in linear elastic problems.
     * - Returning a material property (method 'give'). Only for non-standard elements.
     * - Returning real stress state vector(tensor) at gauss point for 3d - case.
     * - Transforming stiffness from principal orhotrophy axes to
     * system used in given GaussPoint.
     */

protected:
    CS_type cs_type;
    FloatMatrix *localCoordinateSystem;
    FloatArray *helpPlaneNormal;
    // in localCoordinateSystem the unity vectors are stored
    // COLUMWISE (this is exception, but allows faster numerical
    // implementation)

public:

    OrthotropicLinearElasticMaterial(int n, Domain *d) : LinearElasticMaterial(n, d)
    { localCoordinateSystem = NULL;
      helpPlaneNormal = NULL;
      cs_type = unknownCS; }
    ~OrthotropicLinearElasticMaterial()
    { if ( localCoordinateSystem ) { delete localCoordinateSystem; }

      if ( helpPlaneNormal ) { delete helpPlaneNormal; } }

    void giveThermalDilatationVector(FloatArray &answer, GaussPoint *, TimeStep *);

    // identification and auxiliary functions
    const char *giveClassName() const { return "OrthotropicLinearElasticMaterial"; }
    classType giveClassID()         const { return OrthotropicLinearElasticMaterialClass; }
    IRResultType initializeFrom(InputRecord *ir);

    // non-standard - returns time independent material constant
    double   give(int, GaussPoint*);

    void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
                                       MatResponseForm, MatResponseMode, GaussPoint * gp,
                                       TimeStep * atTime);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    /// Computes local 3d stiffness matirx of the receiver
    void give3dLocalMaterialStiffnessMatrix(FloatMatrix & answer,
                                            MatResponseForm, MatResponseMode, GaussPoint * gp,
                                            TimeStep * atTime);


protected:
    /*
     * FloatMatrix* GivePlaneStressStiffMtrx (MatResponseForm,MatResponseMode,GaussPoint* gp,
     *                                      FloatArray* strainIncrement, TimeStep* atTime);
     * FloatMatrix* GivePlaneStrainStiffMtrx (MatResponseForm,MatResponseMode, GaussPoint* gp,
     *            FloatArray* strainIncrement, TimeStep* atTime);
     * FloatMatrix* Give1dStressStiffMtrx (MatResponseForm,MatResponseMode,GaussPoint* gp,
     *           FloatArray* strainIncrement, TimeStep* atTime);
     */
    FloatMatrix *GiveTensorRotationMatrix(GaussPoint *);
    FloatMatrix *GiveRotationMatrix(GaussPoint *);
    /* Diferent response for every element */

    friend class CrossSection;
};


#endif // ortholinearelasticmaterial_h

