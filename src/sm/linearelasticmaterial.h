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

#ifndef linearelasticmaterial_h
#define linearelasticmaterial_h

#include "structuralmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
class GaussPoint;

/**
 * This class is a abstract base class for all linear elastic material models
 * in a finite element problem. A material
 * is an attribute of a domain. It is usually also attribute of many elements.
 * Efficient implementation of services for obtaining characteristic matrices
 * for several material modes is provided depending on other abstract services.
 * These services include following methods: give2dPlateStiffMtrx, give2dPlateStiffMtrx,
 * give3dShellStiffMtrx and  give2dPlaneStressRotStiffMtrx. These methods can be
 * easily and generally implemented using basic services like givePlaneStressStiffMtrx
 * (must be implemented by parent) from which the response can be easily derived.
 * Also general implementation of giveRealStressVector service is provided,
 * computing the stress increment vector from strain increment multiplied by
 * stiffness.
 *
 * Tasks:
 * - Returning standard material stiffness and flexibility marices for 3d-case.
 *   according to current state determined by using data stored in Gausspoint.
 * - Returning a material property (method 'give'). Only for non-standard elements.
 * - Returning real stress state vector(tensor) at gauss point for 3d - case.
 */
class LinearElasticMaterial : public StructuralMaterial
{
public:
    /// Constructor.
    LinearElasticMaterial(int n, Domain *d) : StructuralMaterial(n, d) { }
    /// Destructor.
    virtual ~LinearElasticMaterial() { }

    virtual void giveCharacteristicMatrix(FloatMatrix &answer,
                                  MatResponseForm form,
                                  MatResponseMode mode,
                                  GaussPoint *gp,
                                  TimeStep *atTime);

    /**
     * Method for computing 2d plate stiffness matrix of receiver using generic givePlaneStressStiffMtrx.
     * @param answer Stiffness matrix.
     * @param form Material response form.
     * @param mode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when atTime is current time step).
     */
    void give2dPlateStiffMtrx(FloatMatrix &answer,
                              MatResponseForm form,
                              MatResponseMode mode,
                              GaussPoint *gp,
                              TimeStep *tStep);
    /**
     * Method for computing 3d shell stiffness matrix of receiver using generic givePlaneStressStiffMtrx.
     * @param answer Stiffness matrix.
     * @param form Material response form.
     * @param mode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when atTime is current time step).
     */
    void give3dShellStiffMtrx(FloatMatrix &answer,
                              MatResponseForm form,
                              MatResponseMode mode,
                              GaussPoint *gp,
                              TimeStep *tStep);
    /**
     * Method for computing 2d plane stress (with rotation field)
     * stiffness matrix of receiver using generic givePlaneStressStiffMtrx.
     * @param answer Stiffness matrix.
     * @param form Material response form.
     * @param mode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when atTime is current time step).
     */
    void give2dPlaneStressRotStiffMtrx(FloatMatrix &answer,
                                       MatResponseForm form,
                                       MatResponseMode mode,
                                       GaussPoint *gp,
                                       TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                              const FloatArray &reducedStrain,
                              TimeStep *tStep);

    virtual int giveStressStrainComponentIndOf(MatResponseForm form, MaterialMode mmode, int ind);
    virtual void giveStressStrainMask(IntArray &, MatResponseForm form, MaterialMode mmode) const;
    virtual int hasNonLinearBehaviour() { return 0; }
    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual const char *giveClassName() const { return "LinearElasticMaterial"; }
    virtual classType giveClassID() const { return LinearElasticMaterialClass; }
};
} // end namespace oofem
#endif // linearelasticmaterial_h
