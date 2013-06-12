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

#ifndef structuralmaterial_h
#define structuralmaterial_h

#include "material.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matconst.h"
#include "matstatus.h"
#include "stressstrainprincmode.h"

///@name Input fields for StructuralMaterial
//@{
#define _IFT_StructuralMaterial_referencetemperature "referencetemperature"
//@}

namespace oofem {
#define STRAIN_STEPS 10.0

class GaussPoint;

/**
 * Abstract base class for all "structural" constitutive models. It declares common  services provided
 * by all structural material models. The implementation of these services is partly left on derived classes,
 * which will implement constitutive model dependent part.
 * Some general purpose services are implemented on this level. For details, how to store
 * material model related history variables in integration points, see base class @ref Material documentation.
 *
 * The constitutive model can in general support several material modes (plane stress, plane strain ,... modes).
 * Its capabilities can be examined using hasMaterialModeCapability  service.
 * It is generally assumed, that results obtained from constitutive model services are according to
 * valid material mode. This mode is determined from integration point, which is compulsory parameter of all material
 * services.
 * Structural material introduces several stress/strain modes.
 * Full and reduced formats of stress/strain vectors are introduced for convenience.
 * The full format includes all components, even if they are zero due to stress/strain mode nature,
 * but in the reduced format, only generally nonzero components are stored.
 * (full format must used only if absolutely necessary, to avoid wasting of space. It is used
 * by output routines to print results in general form). Methods for converting vectors between
 * full and reduced format are provided.
 *
 * If in particular mode particular stress component is zero, the corresponding strain is not computed
 * and not stored in reduced vector, and in full vector there is zero value on corresponding position.
 * On the other hand, if some zero strain is imposed,
 * On the other hand, if zero strain component is imposed this condition must be taken into account in geometrical
 * relations (at element level), and corresponding component are included stress/strain reduced vectors.
 *
 * Structural material introduces following basic stress/strain modes
 * - 3d state - all components of general stress/strain vector are generally nonzero.
 *   General 3d strain vector has following components {sig_xx, sig_yy, sig_zz, tau_yz, tau_xz, tau_xy}
 * - plane stress - sig_zz = tau_yz =  tau_xz = 0.
 * - plane strain - eps_z = gamma_xz = gamma_yz = 0.
 *   Note: as already described, if zero strain component is imposed
 *   (Plane strain, ..) this condition must be taken into account in geometrical
 *   relations, and corresponding component has to be included in reduced vector.
 * - 1d uniaxial state - sigma_y = sigma_z = tau_yz = tau_zx = tau_xy  = 0.
 * - 2d beam layer - sigma_y=sigma_z=tau_zy=tau_xy = 0.
 * - 3d shell layer, 2d plate layer - sigma_z = 0.
 *
 * Derived classes can of course extend those modes.
 * Generally speaking, there are following major tasks, covered by declared services.
 * - Computing real stress vector (tensor) at integration point for given strain increment and updating its
 *   state (still temporary state, after overall equilibrium is reached).
 * - Updating its state (final state), when equilibrium has been reached.
 * - Returning its material stiffness (and/or flexibility) matrices for given material mode.
 * - Storing/restoring its context to stream.
 * - Returning its material properties.
 *
 * Structural material services should not be called directly by elements. Instead, they always should
 * pass their requests to corresponding cross section model. Cross section performs all necessary integration over
 * its volume and invokes material model services.
 */
class StructuralMaterial : public Material
{
protected:
    /// Reference temperature (temperature, when material has been built into structure).
    double referenceTemperature;

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    StructuralMaterial(int n, Domain *d) : Material(n, d) { }
    /// Destructor.
    virtual ~StructuralMaterial() { }

    /**
     * Computes the stiffness matrix of receiver in given integration point, respecting its history.
     * The algorithm should use temporary or equilibrium  history variables stored in integration point status
     * to compute and return required result.
     * @param answer Contains result.
     * @param form Material response form.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when atTime is current time step).
     */
    virtual void  giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep);

    /**
     * Computes compliance matrix of receiver in given integration point.
     * @param answer Contains result
     * @param form Material response form.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when atTime is current time step).
     */
    virtual void  giveCharacteristicComplianceMatrix(FloatMatrix &answer,
                                                     MatResponseForm form,
                                                     MatResponseMode mode,
                                                     GaussPoint *gp,
                                                     TimeStep *tStep);

    /**
     * Computes the real stress vector for given total strain and integration point.
     * The total strain is defined as strain computed directly from displacement field at given time.
     * The stress independent parts (temperature, eigenstrains) are subtracted in constitutive
     * driver.
     * The service should use previously reached equilibrium history variables. Also
     * it should update temporary history variables in status according to newly reached state.
     * The temporary history variables are moved into equilibrium ones after global structure
     * equilibrium has been reached by iteration process.
     * @param answer Contains result.
     * @param form Material response form.
     * @param gp Integration point.
     * @param reducedStrain Strain vector in reduced form.
     * @param tStep Current time step (most models are able to respond only when atTime is current time step).
     */
    virtual void giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep) = 0;

    /**
     * Returns a vector of coefficients of thermal dilatation in direction of each material principal (local) axis.
     * @param answer Vector of thermal dilatation coefficients.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when atTime is current time step).
     */
    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
    {
        answer.resize(0);
    }
    /**
     * Returns the reference temperature of receiver.
     */
    double giveReferenceTemperature() { return referenceTemperature; }

    /**
     * Computes reduced strain vector in given integration point, generated by internal processes in
     * material, which are independent on loading in particular integration point.
     * Default implementation takes into account temperature induced strains and eigenstrains.
     * @param answer Returned strain vector.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when atTime is current time step).
     * @param mode Determines response mode (Total or incremental).
     */
    virtual void computeStressIndependentStrainVector(FloatArray &answer,
                                                      GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    // identification and auxiliary functions

    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual const char *giveClassName() const { return "StructuralMaterial"; }
    virtual classType giveClassID() const { return StructuralMaterialClass; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    /**
     * Auxiliary member function that computes principal values of stress/strain vector.
     * @param answer Computed principal values.
     * @param s Stress/strain vector which eigenvalues are computed.
     * @param mode Stress strain principal mode..
     */
    void computePrincipalValues(FloatArray &answer, const FloatArray &s, stressStrainPrincMode mode);
    /**
     * Computes principal values and directions of stress or strain vector.
     * @param answer Computed principal values.
     * @param dir Principal directions (stored column wise).
     * @param s Stress/strain vector.
     * @param mode Stress strain principal mode.
     */
    void computePrincipalValDir(FloatArray &answer, FloatMatrix &dir, const FloatArray &s,
                                stressStrainPrincMode mode);

    /**
     * Computes full 3d material stiffness matrix at given integration point, time, respecting load history
     * in integration point.
     * @param answer Computed results.
     * @param form Material response form.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when atTime is current time step).
     */
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep)
    { _error("give3dMaterialStiffnessMatrix: not implemented "); }

    /**
     * This method returns index of reduced (if form == ReducedForm) or
     * full (if form = FullForm) stress/strain component in Full or Reduced
     * stress/strain vector according to stress/strain mode in given integration point.
     * @param form Material response form.
     * @param mmode Material response mode.
     * @param ind Index of component.
     * @return Component index or 0 or error is generated for unknown Material Mode.
     */
    virtual int giveStressStrainComponentIndOf(MatResponseForm form, MaterialMode mmode, int ind);
    /**
     * This method returns mask of reduced (if form == ReducedForm)
     * or Full (if form==FullForm) stress/strain vector in full or
     * reduced StressStrainVector according to stressStrain mode of given gp.
     * Mask has size of reduced or full stress/strain Vector and  i-th component
     * is index to full or reduced stress/strainVector where corresponding
     * component is mapped.
     * Reduced form is sub-vector (of stress or strain components),
     * where components corresponding to imposed zero stress (plane stress,...)
     * are not included. On the other hand, if zero strain component is imposed
     * (Plane strain, ..) this condition must be taken into account in geometrical
     * relations, and corresponding component has to be included in reduced vector.
     * @param answer Returned mask.
     * @param form Material response form.
     * @param mmode Material response mode.
     * @return For unknown mode error is generated.
     */
    virtual void giveStressStrainMask(IntArray &answer, MatResponseForm form, MaterialMode mmode) const;
    virtual void givePrincipalStressStrainMask(IntArray &answer, MatResponseForm form, MaterialMode mmode) const;
    /**
     * Returns the size of reduced stress/strain vector according to given mode.
     * @param mmode Material response mode.
     */
    virtual int giveSizeOfReducedStressStrainVector(MaterialMode mmode);
    /**
     * Returns the size of reduced principal stress/strain vector according to given mode.
     * @param mmode Material response mode.
     */
    virtual int giveSizeOfReducedPrincipalStressStrainVector(MaterialMode mmode);
    /**
     * Method for subtracting from reduced space strain vector its stress-independent parts
     * (caused by temperature, shrinkage, creep and possibly by other phenomena).
     * Calls StructuralElement::computeStressIndependentStrainVector to obtain stress
     * independent part of strain.
     * @param answer Computed strain vector.
     * @param gp Integration point.
     * @param reducedStrainVector Reduced strain vector.
     * @param tStep Time step.
     * @param mode Determines value mode.
     */
    void giveStressDependentPartOfStrainVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrainVector,
                                               TimeStep *tStep, ValueModeType mode);

    virtual int setIPValue(const FloatArray &value, GaussPoint *gp, InternalStateType type);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *gp);

    /**
     * Computes reduced stress/strain vector from full stress/strain vector.
     * The stress/strain mode is determined form given integration point.
     * @param answer Reduced version of charVector3d.
     * @param gp Integration point.
     * @param charVector3d Full 3d stress/strain vector.
     */
    void giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                         const FloatArray &charVector3d);
    /**
     * Computes full form of stress/strain from its reduced form, based on stress/strain mode
     * stored in given integration point.
     * @param answer Full form of stress/strain vector.
     * @param gp Integration point.
     * @param strainVector Reduced vector.
     */
    void giveFullCharacteristicVector(FloatArray &answer,  GaussPoint *gp,
                                      const FloatArray &strainVector);

protected:
    /**
     * Computes characteristic stiffness matrix corresponding to given material mode
     * (obtained form integration point) by reduction of 3d stiffness matrix.
     * This is general method, how to obtain stiffness matrix corresponding to specific mode
     * from general 3d stiffness. Therefore, it is only necessary to implement algorithm for
     * computing general 3d stiffness. However, this reduction is quite time consuming
     * and if it is possible, it is recommended to provide direct methods for computing
     * particular stiffnesses for supported material modes.
     * @param answer Reduced stiffness.
     * @param form Material response form.
     * @param gp Integration point.
     * @param stiffMtrx3d 3d stiffness matrix (full form) in given integration point.
     */
    void reduceStiffMtrx3d(FloatMatrix &answer, MatResponseForm form, GaussPoint *gp,
                           FloatMatrix &stiffMtrx3d) const;
    /**
     * Computes characteristic compliance matrix corresponding to given material mode
     * (obtained form integration point) by reduction of full 3d compliance matrix.
     * This is general method, how to obtain compliance matrix corresponding to specific mode
     * from general 3d compliance. Therefore, it is only necessary to implement algorithm for
     * computing general 3d compliance. However, this reduction is quite time consuming
     * and if it is possible, it is recommended to provide direct methods for computing
     * particular compliances for supported material modes.
     * @param answer Reduced compliance matrix.
     * @param form Material response form.
     * @param gp Integration point.
     * @param complMtrx3d 3d compliance matrix (full form) in given integration point.
     */
    void reduceComplMtrx3d(FloatMatrix &answer, MatResponseForm form, GaussPoint *gp,
                           FloatMatrix &complMtrx3d) const;
    /**
     * Reduces full 3d stiffness matrix to 2d plane stress matrix.
     * The 3d stiffness should be computed for integration point passed as parameter.
     * @param answer Computed reduced stiffness matrix.
     * @param form Material response form.
     * @param gp Integration point.
     * @param stiffMtrx3d 3d stiffness matrix (in full form).
     */
    void reduceToPlaneStressStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                      GaussPoint *gp, FloatMatrix &stiffMtrx3d) const;
    /**
     * Reduces full 3d stiffness matrix to 2d plane strain matrix.
     * The 3d stiffness should be computed for integration point passed as parameter.
     * Note: as already described, if zero strain component is imposed
     * (Plane strain, ..) this condition must be taken into account in geometrical
     * relations, and corresponding component has to be included in reduced vector.
     * (So plane strain conditions are @f$ \epsilon_z = \gamma_{xz} = \gamma_{yz} = 0 @f$, but relations
     * for @f$\epsilon_z @f$ and @f$ \sigma_z @f$ are included).
     * @param answer Computed reduced stiffness matrix.
     * @param form Material response form.
     * @param gp Integration point.
     * @param stiffMtrx3d 3d stiffness matrix (in full form).
     */
    void reduceToPlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                      GaussPoint *gp, FloatMatrix &stiffMtrx3d) const;
    /**
     * Reduces full 3d stiffness matrix to 1d matrix
     * 1d case: @f$ \sigma_y = \sigma_z = \tau_{yz} = \tau_{zx} = \tau_{xy} = 0 @f$.
     * The 3d stiffness should be computed for integration point passed as parameter.
     * @param answer Computed reduced stiffness matrix.
     * @param form Material response form.
     * @param gp Integration point.
     * @param stiffMtrx3d 3d stiffness matrix (in full form).
     */
    void reduceTo1dStressStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                   GaussPoint *gp, FloatMatrix &stiffMtrx3d) const;
    /**
     * Reduces full 3d stiffness matrix to 2d beam layer matrix.
     * 2dbeamLayer: @f$ \sigma_y = \sigma_z = \tau_{zy} = \tau_{xy} = 0 @f$.
     * The 3d stiffness should be computed for integration point passed as parameter.
     * @param answer Computed reduced stiffness matrix.
     * @param form Material response form.
     * @param gp Integration point.
     * @param stiffMtrx3d 3d stiffness matrix (in full form).
     */
    void reduceTo2dBeamLayerStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                      GaussPoint *gp, FloatMatrix &stiffMtrx3d) const;
    /**
     * Reduces full 3d stiffness matrix to 2d plate layer stiffness matrix.
     * 2dplatelayermode: @f$ \sigma_z = 0 @f$.
     * @param answer Computed reduced stiffness matrix.
     * @param form Material response form.
     * @param gp Integration point.
     * @param stiffMtrx3d 3d stiffness matrix (in full form).
     */
    void reduceTo2dPlateLayerStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                       GaussPoint *gp, FloatMatrix &stiffMtrx3d) const;
    /**
     * Reduces full 3d stiffness matrix to 1d fiber stiffness matrix.
     * 2dplatelayermode: @f$ \sigma_y = \sigma_z = \tau_{yz} = 0 @f$.
     * @param answer Computed reduced stiffness matrix.
     * @param form Material response form.
     * @param gp Integration point.
     * @param stiffMtrx3d 3d stiffness matrix (in full form).
     */
    void reduceTo1dFiberStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                  GaussPoint *gp, FloatMatrix &stiffMtrx3d) const;
    /**
     * Reduces full 3d stiffness matrix to 3d shell layer stiffness matrix.
     * @see StructuralMaterial::reduceTo2dPlateLayerStiffMtrx
     * @param answer Computed reduced stiffness matrix.
     * @param form Material response form.
     * @param gp Integration point.
     * @param stiffMtrx3d 3d stiffness matrix (in full form).
     */
    void reduceTo3dShellLayerStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                       GaussPoint *gp, FloatMatrix &stiffMtrx3d) const;


    /**
     * Reduces full 3d compliance matrix to 2d plane stress matrix.
     * The 3d compliance should be computed for integration point passed as parameter.
     * @param answer Computed reduced stiffness matrix.
     * @param form Material response form.
     * @param gp Integration point.
     * @param complMtrx3d 3d compliance matrix (in full form).
     */
    void reduceToPlaneStressComplMtrx(FloatMatrix &answer, MatResponseForm form,
                                      GaussPoint *gp, FloatMatrix &complMtrx3d) const;
    /**
     * Reduces full 3d compliance matrix to 2d plane strain matrix.
     * The 3d compliance should be computed for integration point passed as parameter.
     * Note: as already described, if zero strain component is imposed
     * (Plane strain, ..) this condition must be taken into account in geometrical
     * relations, and corresponding component has to be included in reduced vector.
     * (So plane strain conditions are @f$ \epsilon_z = \gamma_{xz} = \gamma_{yz} = 0 @f$, but relations
     * for eps_z and sigma_z are included).
     * @param answer Computed reduced compliance matrix.
     * @param form Material response form.
     * @param gp Integration point.
     * @param complMtrx3d 3d compliance matrix (in full form).
     */
    void reduceToPlaneStrainComplMtrx(FloatMatrix &answer, MatResponseForm form,
                                      GaussPoint *gp, FloatMatrix &complMtrx3d) const;
    /**
     * Reduces full 3d compliance matrix to 1d matrix
     * 1d case: @f$ \sigma_y = \sigma_z = \tau_{yz} = \tau_{zx} = \tau_{xy} = 0 @f$.
     * The 3d compliance should be computed for integration point passed as parameter.
     * @param answer Computed reduced compliance matrix.
     * @param form Material response form.
     * @param gp Integration point.
     * @param complMtrx3d 3d compliance matrix (in full form).
     */
    void reduceTo1dStressComplMtrx(FloatMatrix &answer, MatResponseForm form,
                                   GaussPoint *gp, FloatMatrix &complMtrx3d) const;
    /**
     * Reduces full 3d compliance matrix to 2d beam layer matrix.
     * 2dbeamLayer: @f$ \sigma_y = \sigma_z = \tau_{zy} = \tau_{xy} = 0 @f$.
     * @param answer Computed reduced compliance matrix.
     * @param form Material response form.
     * @param gp Integration point.
     * @param complMtrx3d 3d compliance matrix (in full form).
     */
    void reduceTo2dBeamLayerComplMtrx(FloatMatrix &answer, MatResponseForm form,
                                      GaussPoint *gp, FloatMatrix &complMtrx3d) const;
    /**
     * Reduces full 3d compliance matrix to 2d plate layer compliance matrix.
     * 2dplatelayermode: @f$ \sigma_z = 0 @f$.
     * @param answer Computed reduced compliance matrix.
     * @param form Material response form.
     * @param gp Integration point.
     * @param complMtrx3d 3d compliance matrix (in full form).
     */
    void reduceTo2dPlateLayerComplMtrx(FloatMatrix &answer, MatResponseForm form,
                                       GaussPoint *gp, FloatMatrix &complMtrx3d) const;
    /**
     * Reduces full 3d compliance matrix to 3d shell layer compliance matrix.
     * @param answer Computed reduced compliance matrix.
     * @param form Material response form.
     * @param gp Integration point.
     * @param complMtrx3d 3d compliance matrix (in full form).
     */
    void reduceTo3dShellLayerComplMtrx(FloatMatrix &answer, MatResponseForm form,
                                       GaussPoint *gp, FloatMatrix &complMtrx3d) const;
    /**
     * Reduces full 3d compliance matrix to 1d fiber layer compliance matrix.
     * 1dfiber: @f$ \sigma_y = \sigma_z = \tau_{yz} = 0 @f$.
     * @param answer Computed reduced compliance matrix.
     * @param form Material response form.
     * @param gp Integration point.
     * @param complMtrx3d 3d compliance matrix (in full form).
     */
    void reduceTo1dFiberComplMtrx(FloatMatrix &answer, MatResponseForm form,
                                  GaussPoint *gp, FloatMatrix &complMtrx3d) const;


    /**
     * Method for computing plane stress stiffness matrix of receiver.
     * Default implementation computes 3d stiffness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to plane stress stiffness using reduce method described above.
     * However, this reduction is quite time consuming and if it is possible,
     * it is recommended to overload this method and provide direct method for computing
     * particular stiffness matrix.
     * @param answer Stiffness matrix.
     * @param form Material response form.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when atTime is current time step).
     */
    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer,
                                          MatResponseForm form, MatResponseMode mmode, GaussPoint *gp,
                                          TimeStep *tStep);
    /**
     * Method for computing plane strain stiffness matrix of receiver.
     * Default implementation computes 3d stiffness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to plane strain stiffness using reduce method described above.
     * However, this reduction is quite time consuming and if it is possible,
     * it is recommended to overload this method and provide direct method for computing
     * particular stiffness matrix.
     * Note: as already described, if zero strain component is imposed
     * (Plane strain, ..) this condition must be taken into account in geometrical
     * relations, and corresponding component has to be included in reduced vector.
     * (So plane strain conditions are @f$ \epsilon_z = \gamma_{xz} = \gamma_{yz} = 0 @f$, but relations
     * for @f$ \epsilon_z@f$ and @f$\sigma_z@f$ are included).
     * @param answer Stiffness matrix.
     * @param form Material response form.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when atTime is current time step).
     */
    virtual void givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                          MatResponseForm form, MatResponseMode mmode, GaussPoint *gp,
                                          TimeStep *tStep);
    /**
     * Method for computing 1d stiffness matrix of receiver.
     * Default implementation computes 3d stiffness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to 1d stiffness using reduce method described above.
     * However, this reduction is quite time consuming and if it is possible,
     * it is recommended to overload this method and provide direct method for computing
     * particular stiffness matrix.
     * @param answer Stiffness matrix.
     * @param form Material response form.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when atTime is current time step).
     */
    virtual void give1dStressStiffMtrx(FloatMatrix &answer,
                                       MatResponseForm form, MatResponseMode mmode, GaussPoint *gp,
                                       TimeStep *tStep);
    /**
     * Method for computing 2d beam layer stiffness matrix of receiver.
     * Default implementation computes 3d stiffness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to 2d beam layer stiffness using reduce method described above.
     * However, this reduction is quite time consuming and if it is possible,
     * it is recommended to overload this method and provide direct method for computing
     * particular stiffness matrix.
     * @param answer Stiffness matrix.
     * @param form Material response form.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when atTime is current time step).
     */
    virtual void give2dBeamLayerStiffMtrx(FloatMatrix &answer,
                                          MatResponseForm form, MatResponseMode mmode, GaussPoint *gp,
                                          TimeStep *tStep);
    /**
     * Method for computing 2d plate layer stiffness matrix of receiver.
     * Default implementation computes 3d stiffness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to 2d plate layer stiffness using reduce method described above.
     * However, this reduction is quite time consuming and if it is possible,
     * it is recommended to overload this method and provide direct method for computing
     * particular stiffness matrix.
     * @param answer Stiffness matrix.
     * @param form Material response form.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when atTime is current time step).
     */
    virtual void give2dPlateLayerStiffMtrx(FloatMatrix &answer,
                                           MatResponseForm form, MatResponseMode mmode, GaussPoint *gp,
                                           TimeStep *tStep);
    /**
     * Method for computing 3d shell layer stiffness matrix of receiver.
     * Default implementation computes 3d stiffness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to 3d shell layer stiffness using reduce method described above.
     * However, this reduction is quite time consuming and if it is possible,
     * it is recommended to overload this method and provide direct method for computing
     * particular stiffness matrix.
     * @param answer Stiffness matrix.
     * @param form Material response form.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when atTime is current time step).
     */
    virtual void give3dShellLayerStiffMtrx(FloatMatrix & answer,
                                           MatResponseForm form, MatResponseMode mmode, GaussPoint * gp,
                                           TimeStep * tStep);

    /**
     * Method for computing 1d fiber stiffness matrix of receiver.
     * Default implementation computes 3d stiffness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to 1d fiber stiffness using reduce method described above.
     * However, this reduction is quite time consuming and if it is possible,
     * it is recommended to overload this method and provide direct method for computing
     * particular stiffness matrix.
     * @param answer Stiffness matrix.
     * @param form Material response form.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when atTime is current time step).
     */
    virtual void give1dFiberStiffMtrx(FloatMatrix &answer,
                                      MatResponseForm form, MatResponseMode mmode, GaussPoint *gp,
                                      TimeStep *tStep);

    /**
     * Transforms 3d strain vector into another coordinate system.
     * @param answer Transformed strain vector
     * @param base Transformation matrix. There are on each column stored unit vectors of
     * coordinate system (so called base vectors) to which we do transformation. These vectors must
     * be expressed in the same coordinate system as source strainVector.
     * @param strainVector 3d strain.
     * @param transpose Determines if we transpose matrix before transforming.
     */
    void transformStrainVectorTo(FloatArray &answer, const FloatMatrix &base,
                                 const FloatArray &strainVector, bool transpose = false) const;
    /**
     * Transforms 3d stress vector into another coordinate system.
     * @param answer Transformed stress vector.
     * @param base Transformation matrix. There are on each column stored unit vectors of
     * coordinate system (so called base vectors) to which we do transformation. These vectors must
     * be expressed in the same coordinate system as source stressVector.
     * @param stressVector Transformed 3d strain.
     * @param transpose Determines if we transpose matrix before transforming.
     */
    void transformStressVectorTo(FloatArray &answer, const FloatMatrix &base,
                                 const FloatArray &stressVector, bool transpose = false) const;

    /**
     * Computes equivalent of von Mises stress. Returns 0 if six stress components do not exist on the material.
     * @param currentStress Stress vector given by 6 components.
     */
    double computeVonMisesStress(const FloatArray *currentStress);

    /**
     * Computes 3d strain vector transformation matrix from standard vector transformation matrix.
     * @param answer Transformation matrix for strain vector.
     * @param base A (3,3) matrix, where on each column are stored unit direction vectors of
     * local coordinate axes to which we do transformation.
     * @param transpose Determines if we transpose matrix before transforming.
     */
    void giveStrainVectorTranformationMtrx(FloatMatrix &answer, const FloatMatrix &base,
                                           bool transpose = false) const;
    /**
     * Computes 3d stress vector transformation matrix from standard vector transformation matrix.
     * @param answer Transformation matrix for stress vector.
     * @param base A (3,3) matrix, where on each column are stored unit direction vectors of
     * local coordinate axes to which we do transformation.
     * @param transpose Determines if we transpose matrix before transforming.
     */
    void giveStressVectorTranformationMtrx(FloatMatrix &answer, const FloatMatrix &base,
                                           bool transpose = false) const;
    /**
     * Computes 2d stress vector transformation matrix from standard vector transformation matrix.
     * @param answer Transformation matrix for stress vector.
     * @param base A (2,2) matrix, where on each column are stored unit direction vectors of
     * local coordinate axes to which we do transformation.
     * @param transpose Determines if we transpose matrix before transforming.
     */
    void givePlaneStressVectorTranformationMtrx(FloatMatrix &answer, const FloatMatrix &base,
                                                bool transpose = false) const;
    /**
     * Method for sorting newly computed principal values (pVal) and
     * corresponding principal directions (pDir) to be closed
     * to some (often previous) principal directions (toPDir).
     * pDir and toPDir should have eigenvectors stored in columns and normalized.
     * @param pVal New eigenvalues.
     * @param pDir New eigenvectors.
     * @param toPDir Old eigenvector.
     */
    void sortPrincDirAndValCloseTo(FloatArray *pVal, FloatMatrix *pDir, FloatMatrix *toPDir);

    friend class CrossSection;
    friend class StructuralCrossSection;
    friend class SimpleCrossSection;
    friend class LayeredCrossSection;
    friend class RheoChainMaterial;
};
} // end namespace oofem
#endif // structuralmaterial_h
