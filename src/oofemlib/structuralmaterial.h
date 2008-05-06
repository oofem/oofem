/* $Header: /home/cvs/bp/oofem/oofemlib/src/structuralmaterial.h,v 1.11.4.1 2004/04/05 15:19:44 bp Exp $ */
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

//   *********************************
//   *** CLASS STRUCTURAL MATERIAL ***
//   *********************************

#ifndef structuralmaterial_h
#define structuralmaterial_h

#include "material.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

#include "matconst.h"
#include "structuralelement.h"
#include "matstatus.h"
#include "stressstrainprincmode.h"

#define STRAIN_STEPS 10.0

class GaussPoint;


/**
 * Abstract base class for all "structural" constitutive models. It declares common  services provided
 * by all structural material models. The implementation of these services is partly left on derived classes,
 * which will implement constitutive model dependent part.
 * Some general purpose services are implemented on this level. For details, how to store
 * material model related history variables in integration points, see base class \ref Material documentation.
 *
 * The constitutive model can in general support several material modes (plane stress, plane strain ,... modes).
 * Its capabilities can be examined using hasMaterialModeCapability  service.
 * It is generally assumed, that results obtained from constitutive model services are according to
 * valid material mode. This mode is determined from integration point, which is compulsory parameter of all material
 * services.
 * Structural material introduces several stress/strain modes.
 * Full and reduced formats of stress/strain vectors are introduced for convinience.
 * The full format includes all components, even if they are zero due to stress/strain mode nature,
 * but in the reduced format, only generally nonzero components are stored.
 * (full format must used only if absolutely necessary, to avoid vasting of space. It is used
 * by output routines to print results in general form). Methods for converting vectors between
 * full and reduced format are provided.
 *
 * If in particular mode particulat stress component is zero, the corresponding strain is not computed
 * and not stored in reduced vector, and in full vector there is zero value on correspondin position.
 * On the other hand, if some zero strain is imposed,
 * On the other hand, if zero strain component is imposed this condition must be taken into account in geometrical
 * relations (at element level), and corresponding component are included stress/strain reduced vectors.
 *
 * Structural material introduces following basic stress/strain modes
 * <UL>
 * <LI>
 * 3d state - all components of general stress/strain vector are generally nonzero.
 * General 3d strain vector has following components {sig_xx, sig_yy, sig_zz, tau_yz, tau_xz, tau_xy}</LI>
 * <LI>
 * plane stress - sig_zz = tau_yz =  tau_xz = 0.</LI>
 * <LI>
 * plane strain - eps_z = gamma_xz = gamma_yz = 0.
 * Note: as already described, if zero strain component is imposed
 * (Plane strain, ..) this condition must be taken into account in geometrical
 * relations, and corresponding component has to be included in reduced vector.</LI>
 * <LI>
 * 1d uniaxial state - sigma_y = sigma_z = tau_yz = tau_zx = tau_xy  = 0.</LI>
 * <LI>
 * 2d beam layer - sigma_y=sigma_z=tau_zy=tau_xy = 0.</LI>
 * <LI>
 * 3d shell layer, 2d plate layer - sigma_z = 0.</LI>
 * </UL>
 *
 * Derived classes can of course extend those modes.
 * Generally speaking, there are following major tasks, covered by declared services.
 * <UL>
 * <LI>
 * Computing real stress vector (tensor) at integration point for given strain increment and updating its
 * state (still temporary state, after overall equilibrium is reached).</LI>
 * <LI>
 * Updating its state (final state), when equlibrium has been reached.</LI>
 * <LI>
 * Returning its material stiffness (and/or flexibility) matrices for given material mode.</LI>
 * <LI>
 * Storing/restoring its context to stream.</LI>
 * <LI>
 * Returning its material properties.</LI>
 * </UL>
 *
 * Structura material services should not be called directly by elements. Instead, they always should
 * pass their requests to corresponding cross section model. Cross section performs all necessary integration over
 * its volume and invokes material model services.
 */
class StructuralMaterial : public Material
{
    /*
     * This class implements a material in a finite element problem. A material
     * is an attribute of a domain. It is usually also attribute of many elements.
     * DESCRIPTION
     * The attribute 'propertyDictionary' contains all the properties of a mate-
     * rial, like its Young modulus, its mass density or poisson ratio.
     * TASK
     * - Returning standard material stiffness and flexibility marices for 3d-case.
     * according to current state determined by using data stored
     * in Gausspoint.
     * - Returning standard material stiffness for other stress states for point
     * in 3d continua - (2dPlanaStress, 2dPlaneStrain, 1dStress);
     * - Returning a material property (method 'give'). Only for non-standard elements.
     * - Returning real stress state vector(tensor) at gauss point for 3d - case.
     * - Imposing constrains according to stressStrain mode in gp to 3dstiffmessMatrix
     * (function reduceTo).
     * - storing / restoring context
     */
protected:
    /// Raference temperature (temperature, when material has been built into structure).
    double referenceTemperature;
public:

    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n material number
     * @param d domain to which new material will belong
     */
    StructuralMaterial(int n, Domain *d) : Material(n, d) { }
    /// Destructor.
    ~StructuralMaterial()                { }

    // standart matrial stiffness matrices
    /**
     * Computes the stiffness matrix of receiver in given integration point, respecting its history.
     * The algorithm should use temporary or equlibrium  history variables stored in integration point status
     * to compute and return required result.
     * @param answer contains result
     * @param form material response form
     * @param mode  material response mode
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    virtual void  giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);

    // standart matrial compliance matrices
    /**
     * Computes conpliance matrix of receiver in given integration point.
     * @param answer contains result
     * @param form material response form
     * @param mode  material response mode
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    virtual void  giveCharacteristicComplianceMatrix(FloatMatrix &answer,
                                                     MatResponseForm form,
                                                     MatResponseMode mode,
                                                     GaussPoint *gp,
                                                     TimeStep *atTime);


    /**
     * Computes the real stress vector for given total strain and integration point.
     * The total strain is defined as strain computed directly from displacement field at given time.
     * The stress independent parts (temperature, eigen strains) are substracted in constitutive
     * driver.
     * The service should use previously reached equilibrium history variables. Also
     * it should update temporary history variables in status according to newly reached state.
     * The temporary history variables are moved into equilibrium ones after global structure
     * equlibrium has been reached by iteration process.
     * @param answer contains result
     * @param form material response form
     * @param gp integration point
     * @param reducedStrain strain vector in reduced form
     * @param tStep current time step (most models are able to respond only when atTime is current time step)
     */
    virtual void giveRealStressVector(FloatArray & answer, MatResponseForm, GaussPoint *,
                                      const FloatArray &, TimeStep *) = 0;

    // returns a FloatArray(3) of coefficients of thermal dillatation in direction
    // of each (local) axisgiven by principal axis of material
    //
    /**
     * Returns a vector of coefficients of thermal dilatation in direction
     * of each material principal (local) axis.
     * @param answer vector of thermal dilatation coefficients
     * @param gp integration point
     * @param tStep time step (most models are able to respond only when atTime is current time step)
     */
    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *, TimeStep *)
    { answer.resize(0);
      return; }
    /**
     * Returns the reference temperature of receiver.
     */
    double giveReferenceTemperature() { return referenceTemperature; }

    /**
     * Computes reduced strain vector in given integration point, generated by internal processes in
     * material, which are independent on loading in particular integration point.
     * Default implementation takes only into account temperature induced strains.
     * @param answer returned strain vector
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     * @param determines response mode (Total or incremental)
     */
    virtual void computeStressIndependentStrainVector(FloatArray &answer,
                                                      GaussPoint *gp, TimeStep *stepN, ValueModeType mode);

    // identification and auxiliary functions
    /**
     * Requests material mode capability.
     * @param mode material mode requested
     * @return nonzero if available
     */
    virtual int hasMaterialModeCapability(MaterialMode);
    /**
     * Request material extension.
     * @param ext material extension tested
     * @return nonzero if implemented
     */
    virtual int testMaterialExtension(MaterialExtension ext) { return ( ( ext == Material_StructuralCapability ) ? 1 : 0 ); }
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "StructuralMaterial"; }
    /// Returns classType id of receiver.
    classType giveClassID()         const { return StructuralMaterialClass; }
    /**
     * Initializes receiver acording to object description stored in input record.
     * The density of material is read into property dictionary (keyword 'd')
     */
    IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     * @param str string to be filled by input record
     * @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    /**
     * Auxiliary member function that computes principal values of stress/strain vector.
     * @param answer computed principal values
     * @param s stress/strain vector which eigenvalues are computed
     * @param mode stress strain principal mode.
     */
    void computePrincipalValues(FloatArray &answer, const FloatArray &s, stressStrainPrincMode mode);
    /**
     * Computes principal values and directions of stress or strain vector.
     * @param answer computed principal values
     * @param dir principal directions (stored columwise)
     * @param s stress/strain vector
     * @param mode stress strain principal mode.
     */
    void computePrincipalValDir(FloatArray &answer, FloatMatrix &dir, const FloatArray &s,
                                stressStrainPrincMode mode);

    /**
     * Computes full 3d material stiffness matrix at given integration point, time, respecting load history
     * in integration point.
     * @param answer computed results
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseForm, MatResponseMode,
                                               GaussPoint *gp,
                                               TimeStep *atTime)
    { _error("give3dMaterialStiffnessMatrix: not implemented "); }

    /**
     * This method returns index of reduced (if form == ReducedForm) or
     * full (if form = FullForm) stres/strain component in Full or Reduced
     * stress/strain vector according to stress/strain mode in given integration point.
     * @param form material response form
     * @gp integration point
     * @ind index of component
     * @return component index or 0 or error is generated for unknown Material Mode.
     */
    virtual int giveStressStrainComponentIndOf(MatResponseForm, MaterialMode mmode, int);
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
     * @param answer returned mask
     * @param form material response form
     * @param mmode material mode
     * @return for unknown mode error is generated
     */
    virtual void giveStressStrainMask(IntArray & answer, MatResponseForm, MaterialMode mmode) const;
    /**
     * Returns the size of reduced stress/strain vector according to given mode.
     * @param mode material response mode.
     */
    virtual int giveSizeOfReducedStressStrainVector(MaterialMode);
    // FloatArray*  ReduceStressStrainVector (MatResponseForm,GaussPoint* gp,
    //              FloatArray* vec);
    /**
     * Method for sustracting from reduced space strain vector its stress independent parts
     * (caused by temperature, shrinkage, creep and possibly by other phenomenons).
     * Calls StructuralElement::computeStressIndependentStrainVector to obtain stress
     * independent part of strain.
     * @param answer computed strain vector
     * @param gp integration point
     * @param reducedStrainVector reduced strain vector
     * @param stepN time step
     * @param mode determines response mode
     */
    void giveStressDependentPartOfStrainVector(FloatArray &answer, GaussPoint *, const FloatArray &,
                                               TimeStep *, ValueModeType mode);


    /**
     * Returns the integration point corresponding value in Reduced form.
     * @param answer contain corresponding ip value, zero sized if not available
     * @param aGaussPoint integration point
     * @param type determines the type of internal variable
     * @param type determines the type of internal variable
     * @returns nonzero if ok, zero if var not supported
     */
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    /**
     * Returns the mask of reduced indexes of Internal Variable component .
     * @param answer mask of Full VectorSize, with components beeing the indexes to reduced form vectors.
     * @param type determines the internal variable requested (physical meaning)
     * @returns nonzero if ok or error is generated for unknown mat mode.
     */
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);


    /**
     * Returns the type of internal variable (scalar, vector, tensor,...).
     * @param type determines the type of internal variable
     * @returns type of internal variable
     */
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    /**
     * Returns the corresponding integration point  value size in Reduced form.
     * @param type determines the type of internal variable
     * @returns var size, zero if var not supported
     */
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);

    /**
     * Computes reduced stress/strain vector from full stress/strain vector.
     * The stress/strain mode is determined form given integration point.
     * @param answer charVector3d reduced
     * @param gp integration point
     * @param charVector3d full 3d stress/strain  vector
     */
    void giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *,
                                         const FloatArray &charVector3d);
    /**
     * Computes full form of stress/strain from its reduced form, based on stress/strainn mode
     * stored in given integration point.
     * @param answer full form of stress/strain vector
     * @param gp integration point
     * @param strainVector reduced vector
     */
    void giveFullCharacteristicVector(FloatArray &answer,  GaussPoint *,
                                      const FloatArray &);


#ifdef __OOFEG
#endif

protected:
    // virtual MaterialStatus* CreateStatus (GaussPoint* gp);
    /**
     * Computes characteristic stiffness matrix corresponding to given material mode
     * (obtained form integration point) by reduction of 3d stiffness matrix.
     * This is general method, how to obtain stifness matrix corresponding to specific mode
     * from general 3d stiffness. Therefore, it is only necessary to implement algorithm for
     * computing general 3d stiffness. Howewer, this reduction is quite time consuming
     * and if it is possible, it is recomended to provide direct methods for computing
     * particular stifnesses for supported material modes.
     * @param answer reduced stifness
     * @param form material response form
     * @param gp integration point
     * @param stiffMtrx3d 3d stifness matrix (full form) in given integration point
     */
    void reduceStiffMtrx3d(FloatMatrix & answer, MatResponseForm, GaussPoint * gp,
                           FloatMatrix & stiffMtrx3d) const;
    /**
     * Computes characteristic compliance matrix corresponding to given material mode
     * (obtained form integration point) by reduction of full 3d compliance matrix.
     * This is general method, how to obtain compliance matrix corresponding to specific mode
     * from general 3d compliance. Therefore, it is only necessary to implement algorithm for
     * computing general 3d compliance. Howewer, this reduction is quite time consuming
     * and if it is possible, it is recomended to provide direct methods for computing
     * particular compilances for supported material modes.
     * @param answer reduced compliance matrix
     * @param form material response form
     * @param gp integration point
     * @param complMtrx3d 3d compilance matrix (full form) in given integration point
     */
    void reduceComplMtrx3d(FloatMatrix & answer, MatResponseForm, GaussPoint * gp,
                           FloatMatrix & complMtrx3d) const;
    /**
     * Reduces full 3d stifness matrix to 2d plane stress matrix.
     * The 3d stifness should be computed for integration point passed as parameter.
     * @param answer computed 2d plane stress stiffness matrix
     * @param form material response form
     * @param gp integration point
     * @param stiffMtrx3d 3d stifness matrix (in full form)
     */
    void reduceToPlaneStressStiffMtrx(FloatMatrix & answer, MatResponseForm,
                                      GaussPoint * gp, FloatMatrix & stiffMtrx3d) const;
    /**
     * Reduces full 3d stifness matrix to 2d plane strain matrix.
     * The 3d stifness should be computed for integration point passed as parameter.
     * Note: as already described, if zero strain component is imposed
     * (Plane strain, ..) this condition must be taken into account in geometrical
     * relations, and corresponding component has to be included in reduced vector.
     * (So plane strain conditions are eps_z = gamma_xz = gamma_yz = 0, but relations
     * for eps_z and sigma_z are included).
     * @param answer computed 2d plane strain stiffness matrix
     * @param form material response form
     * @param gp integration point
     * @param stiffMtrx3d 3d stifness matrix (in full form)
     */
    void reduceToPlaneStrainStiffMtrx(FloatMatrix & answer, MatResponseForm,
                                      GaussPoint * gp, FloatMatrix & stiffMtrx3d) const;
    /**
     * Reduces full 3d stifness matrix to 1d matrix
     * (1d case ==> sigma_y = sigma_z = tau_yz = tau_zx = tau_xy  = 0.).
     * The 3d stifness should be computed for integration point passed as parameter.
     * @param answer computed 2d plane strain stiffness matrix
     * @param form material response form
     * @param gp integration point
     * @param stiffMtrx3d 3d stifness matrix (in full form)
     */
    void reduceTo1dStressStiffMtrx(FloatMatrix & answer, MatResponseForm,
                                   GaussPoint * gp, FloatMatrix & stiffMtrx3d) const;
    /**
     * Reduces full 3d stifness matrix to 2d beam layer matrix.
     * (2dbeamLayer sigma_y=sigma_z=tau_zy=tau_xy = 0)
     * The 3d stifness should be computed for integration point passed as parameter.
     * @param answer computed 2d plane strain stiffness matrix
     * @param form material response form
     * @param gp integration point
     * @param stiffMtrx3d 3d stifness matrix (in full form)
     */
    void reduceTo2dBeamLayerStiffMtrx(FloatMatrix & answer, MatResponseForm,
                                      GaussPoint * gp, FloatMatrix & stiffMtrx3d) const;
    /**
     * Reduces full 3d stifness matrix to 2d plate layer stifness matrix.
     * (2dplatelayermode assumption sigma_z = 0.)
     * @param answer computed 2d plane strain stiffness matrix
     * @param form material response form
     * @param gp integration point
     * @param stiffMtrx3d 3d stifness matrix (in full form)
     */
    void reduceTo2dPlateLayerStiffMtrx(FloatMatrix & answer, MatResponseForm,
                                       GaussPoint * gp, FloatMatrix & stiffMtrx3d) const;
    /**
     * Reduces full 3d stifness matrix to 1d fiber stifness matrix.
     * (2dplatelayermode assumption sigma_y = sigma_z = tau_yz = 0.)
     * @param answer computed 2d plane strain stiffness matrix
     * @param form material response form
     * @param gp integration point
     * @param stiffMtrx3d 3d stifness matrix (in full form)
     */
    void reduceTo1dFiberStiffMtrx(FloatMatrix & answer, MatResponseForm,
                                  GaussPoint * gp, FloatMatrix & stiffMtrx3d) const;
    /**
     * Reduces full 3d stifness matrix to 3d shell layer stifness matrix.
     * @see StructuralMaterial::reduceTo2dPlateLayerStiffMtrx for reference.
     */
    void reduceTo3dShellLayerStiffMtrx(FloatMatrix & answer, MatResponseForm,
                                       GaussPoint * gp, FloatMatrix & stiffMtrx3d) const;
    /**
     * Reduces full 3d compliance matrix to 2d plane stress matrix.
     * The 3d comliance should be computed for integration point passed as parameter.
     * @param answer computed 2d plane stress compliance matrix
     * @param form material response form
     * @param gp integration point
     * @param complMtrx3d 3d compliance matrix (in full form)
     */
    void reduceToPlaneStressComplMtrx(FloatMatrix & answer, MatResponseForm,
                                      GaussPoint * gp, FloatMatrix & complMtrx3d) const;
    /**
     * Reduces full 3d compliance matrix to 2d plane strain matrix.
     * The 3d compliance should be computed for integration point passed as parameter.
     * Note: as already described, if zero strain component is imposed
     * (Plane strain, ..) this condition must be taken into account in geometrical
     * relations, and corresponding component has to be included in reduced vector.
     * (So plane strain conditions are eps_z = gamma_xz = gamma_yz = 0, but relations
     * for eps_z and sigma_z are included).
     * @param answer computed 2d plane strain compliance matrix
     * @param form material response form
     * @param gp integration point
     * @param complMtrx3d 3d compliance matrix (in full form)
     */
    void reduceToPlaneStrainComplMtrx(FloatMatrix & answer, MatResponseForm,
                                      GaussPoint * gp, FloatMatrix & complMtrx3d) const;
    /**
     * Reduces full 3d compliance matrix to 1d matrix
     * (1d case ==> sigma_y = sigma_z = tau_yz = tau_zx = tau_xy  = 0.).
     * The 3d compliance should be computed for integration point passed as parameter.
     * @param answer computed 2d plane strain compliance matrix
     * @param form material response form
     * @param gp integration point
     * @param complMtrx3d 3d compliance matrix (in full form)
     */
    void reduceTo1dStressComplMtrx(FloatMatrix & answer, MatResponseForm,
                                   GaussPoint * gp, FloatMatrix & complMtrx3d) const;
    /**
     * Reduces full 3d compliance matrix to 2d beam layer matrix.
     * (2dbeamLayer sigma_y=sigma_z=tau_zy=tau_xy = 0)
     * The 3d compliance should be computed for integration point passed as parameter.
     * @param answer computed 2d plane strain compliance matrix
     * @param form material response form
     * @param gp integration point
     * @param complMtrx3d 3d compliance matrix (in full form)
     */
    void reduceTo2dBeamLayerComplMtrx(FloatMatrix & answer, MatResponseForm,
                                      GaussPoint * gp, FloatMatrix & complMtrx3d) const;
    /**
     * Reduces full 3d compliance matrix to 2d plate layer compliance matrix.
     * (2dplatelayermode assumption sigma_z = 0.)
     * @param answer computed 2d plane strain compliance matrix
     * @param form material response form
     * @param gp integration point
     * @param complMtrx3d 3d compliance matrix (in full form)
     */
    void reduceTo2dPlateLayerComplMtrx(FloatMatrix & answer, MatResponseForm,
                                       GaussPoint * gp, FloatMatrix & complMtrx3d) const;
    /**
     * Reduces full 3d compliance matrix to 3d shell layer compliance matrix.
     * @see StructuralMaterial::reduceTo2dPlateLayerComplMtrx for reference.
     */
    void reduceTo3dShellLayerComplMtrx(FloatMatrix & answer, MatResponseForm,
                                       GaussPoint * gp, FloatMatrix & complMtrx3d) const;
    /**
     * Reduces full 3d compliance matrix to 1d fiber layer compliance matrix.
     * (1dfiber assumption sigma_y = sigma_z = tau_yz = 0.)
     * @param answer computed 2d plane strain compliance matrix
     * @param form material response form
     * @param gp integration point
     * @param complMtrx3d 3d compliance matrix (in full form)
     */
    void reduceTo1dFiberComplMtrx(FloatMatrix & answer, MatResponseForm,
                                  GaussPoint * gp, FloatMatrix & complMtrx3d) const;
    /**
     * Method for computing plane stress stifness matrix of receiver.
     * Default implementation computes 3d stifness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to plane stress stiffness using reduce method described above.
     * Howewer, this reduction is quite time consuming and if it is possible,
     * it is recomended to overload this method and provide direct method for computing
     * particular stifness matrix.
     * @param answer stifness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    virtual void givePlaneStressStiffMtrx(FloatMatrix & answer,
                                          MatResponseForm, MatResponseMode, GaussPoint * gp,
                                          TimeStep * atTime);
    /**
     * Method for computing plane strain stifness matrix of receiver.
     * Default implementation computes 3d stifness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to plane strain stiffness using reduce method described above.
     * Howewer, this reduction is quite time consuming and if it is possible,
     * it is recomended to overload this method and provide direct method for computing
     * particular stifness matrix.
     * Note: as already described, if zero strain component is imposed
     * (Plane strain, ..) this condition must be taken into account in geometrical
     * relations, and corresponding component has to be included in reduced vector.
     * (So plane strain conditions are eps_z = gamma_xz = gamma_yz = 0, but relations
     * for eps_z and sigma_z are included).
     * @param answer stifness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    virtual void givePlaneStrainStiffMtrx(FloatMatrix & answer,
                                          MatResponseForm, MatResponseMode, GaussPoint * gp,
                                          TimeStep * atTime);
    /**
     * Method for computing 1d  stifness matrix of receiver.
     * Default implementation computes 3d stifness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to 1d stiffness using reduce method described above.
     * Howewer, this reduction is quite time consuming and if it is possible,
     * it is recomended to overload this method and provide direct method for computing
     * particular stifness matrix.
     * @param answer stifness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    virtual void give1dStressStiffMtrx(FloatMatrix & answer,
                                       MatResponseForm, MatResponseMode, GaussPoint * gp,
                                       TimeStep * atTime);
    /**
     * Method for computing 2d beam layer stifness matrix of receiver.
     * Default implementation computes 3d stifness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to 2d beam layer stiffness using reduce method described above.
     * Howewer, this reduction is quite time consuming and if it is possible,
     * it is recomended to overload this method and provide direct method for computing
     * particular stifness matrix.
     * @param answer stifness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    virtual void give2dBeamLayerStiffMtrx(FloatMatrix & answer,
                                          MatResponseForm, MatResponseMode, GaussPoint * gp,
                                          TimeStep * atTime);
    /**
     * Method for computing 2d plate layer stifness matrix of receiver.
     * Default implementation computes 3d stifness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to 2d plate layer stiffness using reduce method described above.
     * Howewer, this reduction is quite time consuming and if it is possible,
     * it is recomended to overload this method and provide direct method for computing
     * particular stifness matrix.
     * @param answer stifness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    virtual void give2dPlateLayerStiffMtrx(FloatMatrix & answer,
                                           MatResponseForm, MatResponseMode, GaussPoint * gp,
                                           TimeStep * atTime);
    /**
     * Method for computing 3d shell layer stifness matrix of receiver.
     * Default implementation computes 3d stifness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to 3d shell layer stiffness using reduce method described above.
     * Howewer, this reduction is quite time consuming and if it is possible,
     * it is recomended to overload this method and provide direct method for computing
     * particular stifness matrix.
     * @param answer stifness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    virtual void give3dShellLayerStiffMtrx(FloatMatrix & answer,
                                           MatResponseForm, MatResponseMode, GaussPoint * gp,
                                           TimeStep * atTime);

    /**
     * Method for computing 1d fiber stifness matrix of receiver.
     * Default implementation computes 3d stifness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to 1d fiber stiffness using reduce method described above.
     * Howewer, this reduction is quite time consuming and if it is possible,
     * it is recomended to overload this method and provide direct method for computing
     * particular stifness matrix.
     * @param answer stifness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    virtual void give1dFiberStiffMtrx(FloatMatrix & answer,
                                      MatResponseForm, MatResponseMode, GaussPoint * gp,
                                      TimeStep * atTime);
    // some useful transformations

    /**
     * Tranforms 3d strain vector into another coordinate system.
     * @param answer transformed strain vector
     * @param base Transformation matrix.  There are on each column stored unit vectors of
     * coordinate system (so called base vectors) to which we do transformation. These vectors must
     * be expressed in the same coordinate system as source strainVector.
     * @param strainVector transformed 3d strain
     * @param transpose If transpose == 1 we transpose base matrix before transforming
     */
    void transformStrainVectorTo(FloatArray &answer, const FloatMatrix &base,
                                 const FloatArray &strainVector, int transpose = 0) const;
    /**
     * Tranforms 3d stress vector into another coordinate system.
     * @param answer transformed stress vector
     * @param base Transformation matrix.  There are on each column stored unit vectors of
     * coordinate system (so called base vectors) to which we do transformation. These vectors must
     * be expressed in the same coordinate system as source stressVector.
     * @param strainVector transformed 3d strain
     * @param transpose If transpose == 1 we transpose base matrix before transforming
     */
    void transformStressVectorTo(FloatArray &answer, const FloatMatrix &base,
                                 const FloatArray &stressVector, int transpose = 0) const;
    /**
     * Computes 3d strain vector transformation matrix from standart vector transormation matrix.
     * @param answer transformation matrix for strain vector
     * @param base (3,3) matrix, where on each column are stored unit direction vectors of
     * local coordinate axes to which we do transformation.
     * @param If transpose == 1 we transpose base matrix before transforming
     */
    void giveStrainVectorTranformationMtrx(FloatMatrix &answer, const FloatMatrix &base,
                                           int transpose = 0) const;
    /**
     * Computes 3d stress vector transformation matrix from standart vector transormation matrix.
     * @param answer transformation matrix for stress vector
     * @param base (3,3) matrix, where on each column are stored unit direction vectors of
     * local coordinate axes to which we do transformation.
     * @param If transpose == 1 we transpose base matrix before transforming
     */
    void giveStressVectorTranformationMtrx(FloatMatrix &answer, const FloatMatrix &base,
                                           int transpose = 0) const;
    /**
     * Method for sorting newly computed principal values (pVal) and
     * corresponding principal directions (pDir) to be closed
     * to some (often previous) principal directions (toPDir).
     * pDir and toPDir should have eigen vectors stored in columns and normalized.
     * @param pVal new eigenvalues
     * @param pDir new eigen vectors
     * @param toPDir old eigen vector
     */
    void sortPrincDirAndValCloseTo(FloatArray *pVal, FloatMatrix *pDir, FloatMatrix *toPDir);


    // FRIENDS:

    friend class CrossSection;
    friend class StructuralCrossSection;
    friend class SimpleCrossSection;
    friend class LayeredCrossSection;
    friend class MaxwellChainMaterial;
};


#endif // structuralmaterial_h





