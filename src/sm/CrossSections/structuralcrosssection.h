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

#ifndef structuralcrosssection_h
#define structuralcrosssection_h

#include "crosssection.h"
#include "../sm/Materials/structuralmaterial.h"

///@name Input fields for CrossSection
//@{
#define _StructuralCrossSection_MaterialNumber "material"
#define _StructuralCrossSection_czMaterialNumber "czmaterial"
//@}


namespace oofem {
class GaussPoint;
typedef GaussPoint IntegrationPoint;
class Element;
class FloatArray;
class FloatMatrix;

/**
 * Abstract base class for all structural cross section models. It declares commons services provided by all
 * structural cross section models. The implementation of this services is left on derived classes,
 * which will implement cross section model dependent part. However, some general services are
 * implemented here.
 * For information, how to introduce integration points in cross section volume for
 * macro integration point, see @ref CrossSection reference manual.
 *
 * At structural level of cross section or constitutive models are introduced several stress/strain modes.
 * Full and reduced formats of stress/strain vectors are also introduced for convenience.
 * The full format includes all components, even if they are zero due to stress/strain mode nature,
 * but in the reduced format, only generally nonzero components are stored.
 * (full format must used only if absolutely necessary, to avoid wasting of space. It is used
 * by output routines to print results in general form). Methods for converting vectors between
 * full and reduced format are provided.
 * General full strain vector has one of the following forms:
 * -# strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
 * -# For integrated cross section models (2d and 3d beams, plates and general shells)
 *    strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
 */
class OOFEM_EXPORT StructuralCrossSection : public CrossSection
{
public:
    /**
     * Constructor. Creates cross section with given number, belonging to given domain.
     * @param n Cross section number.
     * @param d Domain to which new cross section will belong.
     */
    StructuralCrossSection(int n, Domain *d) : CrossSection(n, d)  { }
    /// Destructor.
    virtual ~StructuralCrossSection() { }

    /**
     * Computes the real stress vector for given strain and integration point.
     * The service should use previously reached equilibrium history variables. Also
     * it should update temporary history variables in status according to newly reached state.
     * The temporary history variables are moved into equilibrium ones after global structure
     * equilibrium has been reached by iteration process.
     * Elements should always pass their requests to their cross section model, which
     * performs necessary integration over its volume and invokes necessary material
     * services for corresponding material model defined for given integration point.
     * @param answer Contains result.
     * @param gp Integration point.
     * @param reducedStrain Strain vector in reduced form.
     * @param tStep Current time step (most models are able to respond only when tStep is current time step).
     */
    //@{
    void giveRealStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    virtual void giveRealStress_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) = 0;
    virtual void giveRealStress_3dDegeneratedShell(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) { OOFEM_ERROR("not implemented"); };
    virtual void giveRealStress_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) = 0;
    virtual void giveRealStress_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) = 0;
    virtual void giveRealStress_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) = 0;
    virtual void giveRealStress_Warping(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) = 0;
    //@}

    /**
     * Method for computing the stiffness matrix.
     * @param answer Stiffness matrix.
     * @param mode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    //@{
    virtual void giveStiffnessMatrix_3d(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    virtual void giveStiffnessMatrix_PlaneStress(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    virtual void giveStiffnessMatrix_PlaneStrain(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    virtual void giveStiffnessMatrix_1d(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    //@}

    /**
     * Computes the generalized stress vector for given strain and integration point.
     * @param answer Contains result.
     * @param gp Integration point.
     * @param generalizedStrain Strain vector in reduced generalized form.
     * @param tStep Current time step (most models are able to respond only when tStep is current time step).
     */
    //@{
    virtual void giveGeneralizedStress_Beam2d(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep) = 0;
    virtual void giveGeneralizedStress_Beam3d(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep) = 0;
    virtual void giveGeneralizedStress_Plate(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep) = 0;
    virtual void giveGeneralizedStress_Shell(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep) = 0;
    virtual void giveGeneralizedStress_MembraneRot(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep) = 0;
    virtual void giveGeneralizedStress_PlateSubSoil(FloatArray &answer, GaussPoint *gp, const FloatArray &generalizedStrain, TimeStep *tStep) = 0;
    //@}

    /**
     * Computes the First Piola-Kirchoff stress vector for a given deformation gradient and integration point.
     * The service should use previously reached equilibrium history variables. Also
     * it should update temporary history variables in status according to newly reached state.
     * The temporary history variables are moved into equilibrium ones after global structure
     * equilibrium has been reached by iteration process.
     * Elements should always pass their requests to their cross section model, which
     * performs necessary integration over its volume and invokes necessary material
     * services for corresponding material model defined for given integration point.
     * @param answer Contains the First Piola-Kirchoff stresses.
     * @param gp Integration point.
     * @param reducedFIncrement Increment of the deformation gradient vector in reduced form. @todo should this then be in a multiplicative way? /JB
     * @param tStep Current time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveFirstPKStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedFIncrement, TimeStep *tStep) = 0;


    /**
     * Computes the Cauchy stress vector for a given increment of deformation gradient and given integration point.
     * The service should use previously reached equilibrium history variables. Also
     * it should update temporary history variables in status according to newly reached state.
     * The temporary history variables are moved into equilibrium ones after global structure
     * equilibrium has been reached by iteration process.
     * Elements should always pass their requests to their cross section model, which
     * performs necessary integration over its volume and invokes necessary material
     * services for corresponding material model defined for given integration point.
     * @param answer Contains the Cauchy stress.
     * @param gp Integration point.
     * @param reducedFIncrement Increment of the deformation gradient vector in reduced form. @todo should this then be in a multiplicative way? /JB
     * @param tStep Current time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveCauchyStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedFIncrement, TimeStep *tStep) = 0;

    /**
     * Computes the Eshelby stress vector. Does not update history variables, this is a postprocesing computation.
     * @param answer Contains the Eshelby stress.
     * @param gp Integration point.
     * @param tStep Current time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveEshelbyStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedvF, TimeStep *tStep) {};

    /**
     * Computes the material stiffness matrix dPdF of receiver in a given integration point, respecting its history.
     * The algorithm should use temporary or equilibrium  history variables stored in integration point status
     * to compute and return required result.
     * Elements should always pass their requests to their cross section model, which
     * performs necessary integration over its volume and invokes necessary material
     * services for corresponding material model defined for given integration point.
     * @param answer Contains result.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;

    /**
     * Computes the material stiffness matrix dCde of receiver in a given integration point, respecting its history.
     * The algorithm should use temporary or equilibrium  history variables stored in integration point status
     * to compute and return required result.
     * Elements should always pass their requests to their cross section model, which
     * performs necessary integration over its volume and invokes necessary material
     * services for corresponding material model defined for given integration point.
     * @param answer Contains result.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveStiffnessMatrix_dCde(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;


    /**
     * Computes the stiffness matrix of receiver in given integration point, respecting its history.
     * The algorithm should use temporary or equilibrium  history variables stored in integration point status
     * to compute and return required result.
     * Elements should always pass their requests to their cross section model, which
     * performs necessary integration over its volume and invokes necessary material
     * services for corresponding material model defined for given integration point.
     * @param answer Contains result.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveCharMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;

    /**
     * Computes the stiffness matrix for 2d beams.
     * @param answer The requested matrix.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void give2dBeamStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    /**
     * Computes the stiffness matrix for 2d beams.
     * @param answer The requested matrix.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void give3dBeamStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;

    /**
     * Method for computing 2d plate stiffness matrix.
     * @param answer Stiffness matrix.
     * @param mode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void give2dPlateStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    /**
     * Method for computing 3d shell stiffness matrix.
     * @param answer Stiffness matrix.
     * @param mode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void give3dShellStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    /**
     * Method for computing membrane stiffness matrix with added drilling stiffness.
     * @param answer Stiffness matrix.
     * @param mode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveMembraneRotStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    /**
     * Method for computing subsoil stiffness matrix for plates.
     * @param answer Stiffness matrix.
     * @param mode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void give2dPlateSubSoilStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    /**
     * Returns modified gradient of stress vector, which is used to
     * bring stresses back to yield surface.
     * Method imposes zeros on places, where zero stress occurs. if energetically connected
     * strain is zero, we do not impose zero there, because stress exist and
     * must be taken into account when computing yield function. In such case
     * a problem is assumed to be full 3d with some explicit strain equal to 0.
     * On the other hand, if some stress is imposed to be zero, we understand
     * such case as subspace of 3d case (like a classical plane stress problem, with no
     * tracing of e_z, sigma_z)
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

    virtual int testCrossSectionExtension(CrossSectExtension ext) { return ( ( ext == CS_StructuralCapability ) ? 1 : 0 ); }

    virtual Material *giveMaterial(IntegrationPoint *ip) { OOFEM_ERROR("Missing implementation"); return NULL; }

    virtual void createMaterialStatus(GaussPoint &iGP) = 0;

    virtual int checkConsistency() = 0;
    virtual Interface *giveMaterialInterface(InterfaceType t, IntegrationPoint *ip) { return NULL; }
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode mode) = 0;
};
} // end namespace oofem
#endif // structuralcrosssection_h
