/* $Header: /home/cvs/bp/oofem/oofemlib/src/linearelasticmaterial.h,v 1.9 2003/04/06 14:08:25 bp Exp $ */
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


//   **************************************
//   *** CLASS LINEAR ELACSTIC MATERIAL ***
//   **************************************

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
 * These sevices include following methods: give2dPlateStiffMtrx, give2dPlateStiffMtrx,
 * give3dShellStiffMtrx and  give2dPlaneStressRotStiffMtrx. These methods can be
 * easily and generally implemented using basic services like givePlaneStressStiffMtrx
 * (must be implemented by parent) from which the response can be easily derived.
 * Also general implementation of giveRealStressVector service is provided,
 * computing the stress increment vector from strain increment multiplied by
 * stiffness.
 */
class LinearElasticMaterial : public StructuralMaterial
{
    /*
     * This class implements a linear elastic material in a finite element problem. A material
     * is an attribute of a domain. It is usually also attribute of many elements.
     *
     * DESCRIPTION
     *
     * TASK
     * - Returning standard material stiffness and flexibility marices for 3d-case.
     * according to current state determined by using data stored
     * in Gausspoint.
     * - Returning a material property (method 'give'). Only for non-standard elements.
     * - Returning real stress state vector(tensor) at gauss point for 3d - case.
     */
public:
    /// Constructor
    LinearElasticMaterial(int n, Domain *d) : StructuralMaterial(n, d) { }
    /// Destructor
    ~LinearElasticMaterial() { }

    /**
     * Computes the stiffness matrix of receiver in given integration point, respecting its history.
     * The algorithm should use temporary or equlibrium  history variables stored in integration point status
     * to compute and return required result. Calls cooresponding material mode related services.
     * @param answer contains result
     * @param form material response form
     * @param mode  material response mode
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    void giveCharacteristicMatrix(FloatMatrix &answer,
                                  MatResponseForm form,
                                  MatResponseMode mode,
                                  GaussPoint *gp,
                                  TimeStep *atTime);

    /**
     * Method for computing 2d plate stiffness matrix of receiver using generic givePlaneStressStiffMtrx.
     * @param answer stiffness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    void give2dPlateStiffMtrx(FloatMatrix &answer,
                              MatResponseForm form,
                              MatResponseMode rMode,
                              GaussPoint *gp,
                              TimeStep *tStep);
    /**
     * Method for computing 3d shell stiffness matrix of receiver using generic givePlaneStressStiffMtrx.
     * @param answer stiffness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */

    void give3dShellStiffMtrx(FloatMatrix &answer,
                              MatResponseForm form,
                              MatResponseMode rMode,
                              GaussPoint *gp,
                              TimeStep *tStep);
    /**
     * Method for computing 2d plane stress (with rotation field)
     * stiffness matrix of receiver using generic givePlaneStressStiffMtrx.
     * @param answer stiffness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    void give2dPlaneStressRotStiffMtrx(FloatMatrix &answer,
                                       MatResponseForm form,
                                       MatResponseMode rMode,
                                       GaussPoint *gp,
                                       TimeStep *tStep);
    /**
     * Computes the real stress vector for given strain increment and integration point.
     * It computes the stress increment vector from strain increment multiplied by
     * stiffness.
     * @param answer contains result
     * @param form material response form
     * @param gp integration point
     * @param reducedStrain strain vector in reduced form
     * @param tStep current time step (most models are able to respond only when atTime is current time step)
     */
    void giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                              const FloatArray &reducedStrain,
                              TimeStep *atTime);

    /**
     * This method returns index of reduced (if form == ReducedForm) or
     * full (if form = FullForm) stres/strain component in Full or Reduced
     * stress/strain vector according to stress/strain mode in given integration point.
     * Supported modes are _PlaneStressRot and _PlaneStressRot, otherwise parent method is called.
     * @param form material response form
     * @gp integration point
     * @ind index of component
     * @return component index or error is generated.
     */
    int  giveStressStrainComponentIndOf(MatResponseForm form, MaterialMode mmode, int ind);
    /**
     * This method returns mask of reduced (if form == ReducedForm)
     * or Full (if form==FullForm) stress/strain vector in full or
     * reduced StressStrainVector acording to stressStrain mode of given gp.
     * Mask has size of reduced or full stress/strain Vector and  i-th component
     * is index to full or reduced stress/strainVector where corresponding
     * component is mapped.
     * Reduced form is sub-vector (of stress or strain components),
     * where components corresponding to imposed zero stress (plane stress,...)
     * are not included. On the other hand, if zero strain component is imposed
     * (Plane strain, ..) this condition must be taken into account in geometrical
     * relations, and corresponding component has to be included in reduced vector.
     *
     * Supported modes are _PlaneStressRot and _PlaneStressRot, otherwise parent method is called.
     *
     * @param answer returned mask
     * @param form material response form
     * @param gp integration point
     * @return for unknown mode error is generated
     */
    void giveStressStrainMask(IntArray &, MatResponseForm form, MaterialMode mmode) const;
    // identification and auxiliary functions
    int hasNonLinearBehaviour()   { return 0; }
    /**
     * Tests for particular material mode capability.
     * Supported are _3dMat, _PlaneStress, _PlaneStrain, _1dMat, _2dPlateLayer, _2dBeamLayer,
     * _3dShellLayer, _2dPlate, _3dShell and _PlaneStressRot porvided that all
     * abstract services are correctly implemented by derived class.
     * @param mode material mode requested
     * @return nonzero if available
     */
    int hasMaterialModeCapability(MaterialMode mode);
    /// Returns "LinearElasticMaterial" - class name of the receiver.
    const char *giveClassName() const { return "LinearElasticMaterial"; }
    /// Returns LinearElasticMaterialClass - classType id of receiver.
    classType giveClassID()         const { return LinearElasticMaterialClass; }

protected:
};

} // end namespace oofem
#endif // linearelasticmaterial_h
