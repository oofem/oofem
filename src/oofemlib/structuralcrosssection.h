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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef structuralcrosssection_h
#define structuralcrosssection_h

#include "crosssection.h"
#include "structuralmaterial.h"
#include "gausspnt.h"
#include "element.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
/**
 * Abstract base class for all structural cross section models. It declares commons services provided by all
 * structural cross section models. The implementation of this services is left on derived classes,
 * which will implement cross section model dependent part. Howewer, some general services are
 * implemented here.
 * For information, how to introduce integration points in cross section volume for
 * macro integration point, see @ref CrossSection reference manual.
 *
 * At structural level of cross section or constitutive models are introduced several stress/strain modes.
 * Full and reduced formats of stress/strain vectors are also introduced for convinience.
 * The full format includes all components, even if they are zero due to stress/strain mode nature,
 * but in the reduced format, only generally nonzero components are stored.
 * (full format must used only if absolutely necessary, to avoid vasting of space. It is used
 * by output routines to print results in general form). Methods for converting vectors between
 * full and reduced format are provided.
 * General full strain vector has one of the following forms:
 * -# strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
 * -# For integrated cross section models (2d and 3d beams, plates and general shells)
 *    strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
 */
class StructuralCrossSection : public CrossSection
{
    /*
     * This class implements a structural cross section in a finite element problem. A cross
     * section  is an attribute of a domain. It is usually also attribute of many
     * elements. This class is base class for SimpleCrossSection, LayeredCrossSection,
     * FibredCrossSection, ....
     *
     * DESCRIPTION
     * The attribute 'propertyDictionary' contains all the properties of a
     * cross section, like its area or thickness.
     *
     * TASK
     * - Returning standard material stiffness marices (like 3dstress-strain, 2d plane ,
     *   plate, 3dbeam, 2d beam ..) according to current state determined by parametr
     *   StressMode by calling gp->material->GiveMaterialStiffnessMatrix (....) and by
     *   possible modifiing returned matrix. (for example in layerde mode aproach
     *   each layer  is asked for 3dMatrialStiffnes and this is integrated for example
     *   over thickness for plate bending broblems)
     * - Returning RealStress state in gauss point and for given Stress mode.
     * - Returning a properties of cross section like thickness or area.
     */
public:
    /**
     * Constructor. Creates cross section with given number, belonging to given domain.
     * @param n Cross section number.
     * @param d Domain to which new cross section will belong.
     */
    StructuralCrossSection(int n, Domain *d) : CrossSection(n, d) { }
    /// Destructor.
    ~StructuralCrossSection() { }

    /**
     * Computes the real stress vector for given strain and integration point.
     * The service should use previously reached equilibrium history variables. Also
     * it should update temporary history variables in status according to newly reached state.
     * The temporary history variables are moved into equilibrium ones after global structure
     * equlibrium has been reached by iteration process.
     * Elements should always pass their requests to their cross section model, which
     * performs necessary integration over its volume and invokes necessary material
     * services for corresponding material model defined for given integration point.
     * @param answer Contains result.
     * @param form Material response form.
     * @param gp Integration point.
     * @param reducedStrainIncrement Strain increment vector in reduced form.
     * @param tStep Current time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveRealStresses(FloatArray & answer, MatResponseForm form,
                                  GaussPoint *gp, const FloatArray &reducedStrainIncrement, TimeStep *tStep);

    /**
     * Computes the stiffness matrix of receiver in given integration point, respecting its history.
     * The algorithm should use temporary or equlibrium  history variables stored in integration point status
     * to compute and return required result.
     * Elements should always pass their requests to their cross section model, which
     * performs necessary integration over its volume and invokes necessary material
     * services for corresponding material model defined for given integration point.
     * @param answer Contains result.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveCharMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp,
                                                 TimeStep *tStep);

    // next function is intended to be used if some object would like to obtain
    // char matrix form different material which is not associated with gp and its element.
    // (mainly for obtaining linear elastic matrix)
    // stress-strain mode is taken from gp.
    // NORMALLY - PLEASE USE GiveCharMaterialStiffnessMatrix function
    /**
     * Computes the stiffness matrix of receiver in given integration point, respecting its history.
     * The algorithm should use temporary or equlibrium history variables stored in integration point status
     * to compute and return required result.
     * Passed material mode is always used instead of mode of given integration point.
     * (Therefore, this function should be used if some object would like to obtain
     * char matrix with different material form than that is associated with gp and its element)
     * But please use GiveCharMaterialStiffnessMatrix, this service will not be supported in future
     * releases.
     *
     * Elements should always pass their requests to their cross section model, which
     * performs necessary integration over its volume and invokes necessary material
     * services for corresponding material model defined for given integration point.
     * @param answer Contains result.
     * @param form Material response form.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param mat Material to evaluate for.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveCharMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                                   MatResponseForm form, MatResponseMode mode,
                                                   GaussPoint *gp, StructuralMaterial *mat,
                                                   TimeStep *tStep);
    /**
     * Computes the compliance matrix of receiver in given integration point, respecting its history.
     * The algorithm should use temporary or equlibrium  history variables stored in integration point status
     * to compute and return required result.
     * Elements should always pass their requests to their cross section model, which
     * performs necessary integration over its volume and invokes necessary material
     * services for corresponding material model defined for given integration point.
     * @param answer Contains result.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveCharMaterialComplianceMatrix(FloatMatrix &answer,
                                                  MatResponseMode mode, GaussPoint *gp,
                                                  TimeStep *tStep);

    // next function is intendet to be used if we would like to obtain
    // char matrix form different material which is not associated with gp and its element.
    // (mainly for obtaining linear elastic matrix)
    // stress-strain mode is taken from gp.
    // NORMALLY - PLEASE USE GiveCharMaterialStiffnessMatrix function
    /**
     * Computes the compliance matrix of receiver in given integration point, respecting its history.
     * The algorithm should use temporary or equlibrium  history variables stored in integration point status
     * to compute and return required result.
     * Passed material mode is always used instead of mode of given integration point.
     * (Therefore, this function should be used if some object would like to obtain
     * char matrix with different material form than that is associated with gp and its element)
     * But please use giveCharMaterialComplianceMatrix, this service will not be supported in future
     * releases.
     *
     * Elements should always pass their requests to their cross section model, which
     * performs necessary integration over its volume and invokes necessary material
     * services for corresponding material model defined for given integration point.
     * @param answer Contains result.
     * @param form Material response form.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param mat Material to evalute for.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveCharMaterialComplianceMatrixOf(FloatMatrix &answer,
                                                    MatResponseForm form,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, StructuralMaterial *mat,
                                                    TimeStep *tStep);

    /**
     * Computes reduced strain vector not dependent on sresses in given integration point. 
     * Returned vector is generated by temperature or shrinkage effects, for example.
     * The load mode (Incremental or Total Load form) passed as parameter is taken into account.
     * Depend on load form, tre resulting strain is total strain or its increment from previous
     * step.
     * @param answer Stress independent strain vector.
     * @param gp Integration point.
     * @param mode Determines load mode.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     * @param mode Determines the response mode.
     */
    virtual void computeStressIndependentStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    /**
     * Computes reduced stress/strain vector from full stress/strain vector.
     * The stress/strain mode is determined form given integration point.
     * @param answer charVector3d reduced.
     * @param gp Integration point.
     * @param charVector3d Full 3d stress/strain vector.
     */
    virtual void giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                                 const FloatArray &charVector3d);
    /**
     * Computes full form of stress/strain from its reduced form, based on stress/strain mode
     * stored in given integration point.
     * @param answer Full form of stress/strain vector.
     * @param gp Integration point.
     * @param strainVector Reduced vector.
     */
    virtual void giveFullCharacteristicVector(FloatArray &answer,  GaussPoint *gp,
                                              const FloatArray &strainVector);
    /**
     * Returns modified gradient of stress vector, which is used to
     * bring stresses back to yield surface.
     * Method imposes zeros on places, where zero stress occurs. if energetically connected
     * strain is zero, we do not impose zero there, because stress exist and
     * must be taken into account when computing yeld function. In such case
     * a problem is assumed to be full 3d with some explicit strain equal to 0.
     * On the other hand, if some stress is imposed to be zero, we understand
     * such case as subspace of 3d case (like a classical plane stess problem, with no
     * tracing of ez, sigma_z)
     * @param gp Integration point.
     * @param gradientStressVector3d General 3d stress gradient.
     */
    virtual FloatArray *imposeStressConstrainsOnGradient(GaussPoint *gp, FloatArray *gradientStressVector3d);
    /**
     * Returns modified gradient of strain vector, which is used to compute plastic strain increment.
     * Imposes zeros on places, where zero strain occurs or energetically connected stress
     * is prescribed to be zero.
     * @see imposeStressConstrainsOnGradient
     * @param gp Integration point.
     * @param gradientStressVector3d General 3d stress gradient.
     */
    virtual FloatArray *imposeStrainConstrainsOnGradient(GaussPoint *gp, FloatArray *gradientStressVector3d);

    /**
     * This method returns mask of reduced(if form == ReducedForm)
     * or Full(if form==FullForm) stressStrain vector in full or
     * reduced StressStrainVector acording to stressStrain mode of given gp.
     * Mask has size of reduced or full StressStrain Vector and  i-th component
     * is index to full or reduced StressStrainVector where corresponding
     * stressStrain resides.
     * @param answer Assembled mask.
     * @param form Material response form.
     * @param mode Determines load mode.
     * @param mat Material to evalute for.
     */
    virtual void giveStressStrainMask(IntArray &answer, MatResponseForm form, MaterialMode mode, StructuralMaterial *mat) const;

    // identification and auxiliary functions
    const char *giveClassName() const { return "StructuralCrossSection"; }
    classType giveClassID() const { return StructuralCrossSectionClass; }

    /**
     * Requests cross section model mode capability.
     * @param ext Cross section to test.
     * @return Nonzero if available.
     */
    int testCrossSectionExtension(CrossSectExtension ext) { return ( ( ext == CS_StructuralCapability ) ? 1 : 0 ); }
    //virtual int hasStructuralCapability () {return 1;}

protected:
    /**
     * For internal usage by cross section model.
     * It is direct interface to material model service giveCharacteristicMatrix.
     * @see Material::giveCharacteristicMatrix
     */
    void giveMaterialStiffnessMatrix(FloatMatrix &answer,
                                     MatResponseForm form,
                                     MatResponseMode mode,
                                     GaussPoint *gp,
                                     TimeStep *tStep);

    /**
     * For internal usage by cross section model.
     * It is direct interface to material model service giveCharacteristicMatrix.
     * Material model passed as parameter is used instead of material, to which given integration point
     * belongs to. It should always be the same material as integration point belongs to, because in
     * integration point are stored load history variables related only to its associated material model.
     * Different model can be used only if it does not depend on any internal history variables, like
     * linear elastic material. Accessing load history variables, which are not in integration point status
     * can lead to segmentation fault error.
     * @see Material::giveCharacteristicMatrix
     * @param answer Contains result.
     * @param form Material response form.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param mat Pointer to material model.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                               MatResponseForm form,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               StructuralMaterial *mat,
                                               TimeStep *tStep);
    friend class StructuralMaterial;
};
} // end namespace oofem
#endif // structuralcrosssection_h

