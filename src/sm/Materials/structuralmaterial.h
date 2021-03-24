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

#ifndef structuralmaterial_h
#define structuralmaterial_h

#include "material.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "mathfem.h"
#include "matconst.h"
#include "matstatus.h"
#include "../stressstrainprincmode.h"
#include "valuemodetype.h"
#include <vector>

///@name Input fields for StructuralMaterial
//@{
#define _IFT_StructuralMaterial_referencetemperature "referencetemperature"
#define _IFT_StructuralMaterial_talpha "talpha"
//@}

namespace oofem {
#define STRAIN_STEPS 10.0

class GaussPoint;
///@todo Update documentation
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
 * - Computing real/second PK stress vector (tensor) at integration point for given strain increment and updating its
 *   state (still temporary state, after overall equilibrium is reached).
 * - Updating its state (final state), when equilibrium has been reached.
 * - Returning its material stiffness (and/or flexibility) matrices for given material mode.
 * - Storing/restoring its context to stream.
 * - Returning its material properties.
 *
 * Structural material services should not be called directly by elements. Instead, they always should
 * pass their requests to corresponding cross section model. Cross section performs all necessary integration over
 * its volume and invokes material model services.
 *
 * @author almost everyone
 * @author Jim Brouzoulis
 * @author Mikael Öhman
 */
class OOFEM_EXPORT StructuralMaterial : public Material
{
protected:
    /// Reference temperature (temperature, when material has been built into structure).
    double referenceTemperature = 0.;

public:
    /// Voigt index map
    static std :: array< std :: array< int, 3 >, 3 >vIindex;

    /// Symmetric Voigt index map
    static std :: array< std :: array< int, 3 >, 3 >svIndex;

    static int giveSymVI(int ind1, int ind2) { return svIndex [ ind1 - 1 ] [ ind2 - 1 ]; }
    static int giveVI(int ind1, int ind2) { return vIindex [ ind1 - 1 ] [ ind2 - 1 ]; }

    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    StructuralMaterial(int n, Domain *d);

    bool hasMaterialModeCapability(MaterialMode mode) const override;
    const char *giveClassName() const override { return "StructuralMaterial"; }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    /**
     * Computes the stiffness matrix for giveRealStressVector of receiver in given integration point, respecting its history.
     * The algorithm should use temporary or equilibrium  history variables stored in integration point status
     * to compute and return required result.
     * @param answer Contains result.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveStiffnessMatrix(FloatMatrix &answer,
                                     MatResponseMode mode,
                                     GaussPoint *gp,
                                     TimeStep *tStep);

    /**
     * Computes the real stress vector for given total strain and integration point.
     * The total strain is defined as strain computed directly from displacement field at given time.
     * The stress independent parts (temperature, eigenstrains) are subtracted in constitutive driver.
     * The service should use previously reached equilibrium history variables. Also
     * it should update temporary history variables in status according to newly reached state.
     * The temporary history variables are moved into equilibrium ones after global structure
     * equilibrium has been reached by iteration process.
     * @param answer Stress vector in reduced form. For large deformations it is treated as the second Piola-Kirchoff stress.
     * @param gp Integration point.
     * @param reducedStrain Strain vector in reduced form. For large deformations it is treated as the Green-Lagrange strain.
     * @param tStep Current time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);
    /// Default implementation relies on giveRealStressVector for second Piola-Kirchoff stress
    virtual FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const;
    /// Default implementation relies on giveRealStressVector_3d
    virtual FloatArrayF< 4 >giveRealStressVector_PlaneStrain(const FloatArrayF< 4 > &strain, GaussPoint *gp, TimeStep *tStep) const;
    /// Iteratively calls giveRealStressVector_3d to find the stress controlled equal to zero·
    virtual FloatArray giveRealStressVector_StressControl(const FloatArray &reducedE, const IntArray &strainControl, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatArray giveRealStressVector_ShellStressControl(const FloatArray &reducedE, const IntArray &strainControl, GaussPoint *gp, TimeStep *tStep) const;
    /// Default implementation relies on giveRealStressVector_StressControl
    virtual FloatArrayF< 3 >giveRealStressVector_PlaneStress(const FloatArrayF< 3 > &reducedE, GaussPoint *gp, TimeStep *tStep) const;
    /// Default implementation relies on giveRealStressVector_StressControl
    virtual FloatArrayF< 1 >giveRealStressVector_1d(const FloatArrayF< 1 > &reducedE, GaussPoint *gp, TimeStep *tStep) const;

    /// Default implementation relies on giveRealStressVector_StressControl
    virtual FloatArrayF< 2 >giveRealStressVector_Warping(const FloatArrayF< 2 > &reducedE, GaussPoint *gp, TimeStep *tStep) const;
    /// Default implementation relies on giveRealStressVector_StressControl
    virtual FloatArrayF< 2 >giveRealStressVector_2dBeamLayer(const FloatArrayF< 2 > &reducedE, GaussPoint *gp, TimeStep *tStep) const;
    /// Default implementation relies on giveRealStressVector_StressControl
    virtual FloatArrayF< 5 >giveRealStressVector_PlateLayer(const FloatArrayF< 5 > &reducedE, GaussPoint *gp, TimeStep *tStep) const;
    /// Default implementation relies on giveRealStressVector_StressControl
    virtual FloatArrayF< 3 >giveRealStressVector_Fiber(const FloatArrayF< 3 > &reducedE, GaussPoint *gp, TimeStep *tStep) const;

    /// Default implementation is not provided
    virtual FloatArrayF< 3 >giveRealStressVector_2dPlateSubSoil(const FloatArrayF< 3 > &reducedE, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatArrayF< 6 >giveRealStressVector_3dBeamSubSoil(const FloatArrayF< 6 > &reducedE, GaussPoint *gp, TimeStep *tStep) const;

    /**
     * @name Methods associated with the First PK stress tensor.
     * Computes the first Piola-Kirchhoff stress vector for given total deformation gradient and integration point.
     * The total deformation gradient is computed directly from displacement field at the given time step.
     * The stress independent parts (temperature, eigenstrains) are subtracted in constitutive
     * driver.
     * The service should use previously reached equilibrium history variables. Also
     * it should update temporary history variables in status according to newly reached state.
     * The temporary history variables are moved into equilibrium ones after global structure
     * equilibrium has been reached by iteration process.
     * The First Piola-Kirchhoff stress vector is computed in Total Lagrangian mode
     *
     * @param answer Contains result.
     * @param gp Integration point.
     * @param reducedF Deformation gradient in in reduced form.
     * @param tStep Current time step (most models are able to respond only when tStep is current time step).
     */
    //@{
    /// Default implementation relies on giveRealStressVector for second Piola-Kirchoff stress
    virtual FloatArrayF< 9 >giveFirstPKStressVector_3d(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const;
    /// Default implementation relies on giveFirstPKStressVector_3d
    virtual FloatArrayF< 5 >giveFirstPKStressVector_PlaneStrain(const FloatArrayF< 5 > &vF, GaussPoint *gp, TimeStep *tStep) const;
    /// Iteratively calls giveRealStressVector_3d to find the stress controlled equal to zero·
    virtual FloatArray giveFirstPKStressVector_StressControl(const FloatArray &reducedvF, const IntArray &FControl, GaussPoint *gp, TimeStep *tStep) const;
    /// Default implementation relies on giveFirstPKStressVector_StressControl
    virtual FloatArrayF< 4 >giveFirstPKStressVector_PlaneStress(const FloatArrayF< 4 > &vF, GaussPoint *gp, TimeStep *tStep) const;
    /// Default implementation relies on giveFirstPKStressVector_StressControl
    virtual FloatArrayF< 1 >giveFirstPKStressVector_1d(const FloatArrayF< 1 > &vF, GaussPoint *gp, TimeStep *tStep) const;
    //@}

    /**
     * @name Methods associated with the Cauchy stress tensor.
     * Computes the Cauchy stress vector for given increment of deformation gradient and given integration point.
     * The increment of deformation gradient is computed directly from displacement field at the given time step
     * and it is computed wrt configuration which was reached in the last step.
     * The stress independent parts (temperature, eigenstrains) are subtracted in constitutive
     * driver.
     * The service should use previously reached equilibrium history variables. Also
     * it should update temporary history variables in status according to newly reached state.
     * The temporary history variables are moved into equilibrium ones after global structure
     * equilibrium has been reached by iteration process.
     * The Cauchy stress vector is computed in Updated Lagrangian mode
     *
     * @param answer Contains result.
     * @param gp Integration point.
     * @param reducedF Deformation gradient in in reduced form.
     * @param tStep Current time step (most models are able to respond only when tStep is current time step).
     */
    //@{
    virtual void giveCauchyStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedF, TimeStep *tStep)
    { OOFEM_ERROR("not implemented "); }
    virtual void giveCauchyStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedF, TimeStep *tStep)
    { OOFEM_ERROR("not implemented "); }
    virtual void giveCauchyStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedF, TimeStep *tStep)
    { OOFEM_ERROR("not implemented "); }
    virtual void giveCauchyStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedF, TimeStep *tStep)
    { OOFEM_ERROR("not implemented "); }
    //@}

    /**
     * Prototype for computation of Eshelby stress. No default implementation is provided.
     */
    virtual void giveEshelbyStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedF, TimeStep *tStep);

    static FloatMatrixF< 9, 9 >convert_dSdE_2_dPdF_3D(const FloatMatrixF< 6, 6 > &dSdE, const FloatArrayF< 6 > &S, const FloatArrayF< 9 > &F);
    static FloatMatrixF< 5, 5 >convert_dSdE_2_dPdF_PlaneStrain(const FloatMatrixF< 4, 4 > &dSdE, const FloatArrayF< 4 > &S, const FloatArrayF< 5 > &F);
    static FloatMatrixF< 4, 4 >convert_dSdE_2_dPdF_PlaneStress(const FloatMatrixF< 3, 3 > &dSdE, const FloatArrayF< 3 > &S, const FloatArrayF< 4 > &F);
    static FloatMatrixF< 1, 1 >convert_dSdE_2_dPdF_1D(const FloatMatrixF< 1, 1 > &dSdE, const FloatArrayF< 1 > &S, const FloatArrayF< 1 > &F);

    /**
     * Returns a vector of coefficients of thermal dilatation in direction of each material principal (local) axis.
     * @param answer Vector of thermal dilatation coefficients.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual FloatArrayF< 6 >giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const;
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
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     * @param mode Determines response mode (Total or incremental).
     */
    virtual FloatArray computeStressIndependentStrainVector(GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const;
    FloatArrayF< 6 >computeStressIndependentStrainVector_3d(GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const;
    /// Common functions for convenience
    //@{
    /**
     * Auxiliary member function that computes principal values of stress/strain vector.
     * @param answer Computed principal values.
     * @param s Stress/strain vector which eigenvalues are computed.
     * @param mode Stress strain principal mode..
     */
    static void computePrincipalValues(FloatArray &answer, const FloatArray &s, stressStrainPrincMode mode);
    /**
     * Compute principal values of a stress/strain matrix. Matrix must be symmetric.
     * @param answer Computed principal values.
     * @param s Symmetric matrix.
     * @return Sorted eigenvalues.
     */
    static FloatArrayF< 3 >computePrincipalValues(const FloatMatrixF< 3, 3 > &s);

    /**
     * Compute principal values from invariants.
     * @param I1 First invariant
     * @param I2 Second invariant
     * @param I3 Third invariant
     * @return Sorted list of principal values.
     */
    static FloatArrayF< 3 >computePrincipalValues(double I1, double I2, double I3);
    /**
     * Computes principal values and directions of stress or strain vector.
     * @param answer Computed principal values.
     * @param dir Principal directions (stored column wise).
     * @param s Stress/strain vector.
     * @param mode Stress strain principal mode.
     */
    static void computePrincipalValDir(FloatArray &answer, FloatMatrix &dir, const FloatArray &s, stressStrainPrincMode mode);
    /**
     * Computes principal values and directions of stress or strain matrix (must be symmetric).
     * @param dir Principal directions (stored column wise).
     * @param s Stress/strain matrix.
     * @return Computed principal values.
     */
    static std :: pair< FloatArrayF< 3 >, FloatMatrixF< 3, 3 > >computePrincipalValDir(const FloatMatrixF< 3, 3 > &s);

    /**
     * Computes split of receiver into deviatoric and volumetric part.
     * @param dev Deviatoric part.
     * @param s Input vector
     * @return Volumetric part (diagonal components divided by 3).
     */
    static FloatArrayF< 6 >computeDeviator(const FloatArrayF< 6 > &s);
    static std :: pair< FloatArrayF< 6 >, double >computeDeviatoricVolumetricSplit(const FloatArrayF< 6 > &s);
    static FloatArrayF< 6 >computeDeviatoricVolumetricSum(const FloatArrayF< 6 > &dev, double mean);

    static FloatArrayF< 6 >applyDeviatoricElasticCompliance(const FloatArrayF< 6 > &stress, double EModulus, double nu);
    static FloatArrayF< 6 >applyDeviatoricElasticCompliance(const FloatArrayF< 6 > &stress, double GModulus);

    static FloatArrayF< 6 >applyDeviatoricElasticStiffness(const FloatArrayF< 6 > &strain, double EModulus, double nu);
    static FloatArrayF< 6 >applyDeviatoricElasticStiffness(const FloatArrayF< 6 > &strain, double GModulus);

    static FloatArrayF< 6 >applyElasticStiffness(const FloatArrayF< 6 > &strain, double EModulus, double nu);
    static FloatArrayF< 6 >applyElasticCompliance(const FloatArrayF< 6 > &stress, double EModulus, double nu);

    static double computeStressNorm(const FloatArrayF< 6 > &stress);

    static double computeFirstInvariant(const FloatArrayF< 6 > &s);
    static double computeSecondStressInvariant(const FloatArrayF< 6 > &s);
    static double computeThirdStressInvariant(const FloatArrayF< 6 > &s);

    static double computeFirstCoordinate(const FloatArrayF< 6 > &s);
    static double computeSecondCoordinate(const FloatArrayF< 6 > &s);
    static double computeThirdCoordinate(const FloatArrayF< 6 > &s);
    //@}

    /**
     * Computes full 3d material stiffness matrix at given integration point, time, respecting load history
     * in integration point.
     * @param answer Computed results.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual FloatMatrixF< 6, 6 >give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
    { OOFEM_ERROR("not implemented "); }


    virtual FloatMatrixF< 9, 9 >give3dMaterialStiffnessMatrix_dPdF(MatResponseMode mode,
                                                                   GaussPoint *gp, TimeStep *tStep) const;

    virtual void give3dMaterialStiffnessMatrix_dCde(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep);

    /**
     * Returns a mask of the vector indicies corresponding to components in a general
     * (non-symmetric) second order tensor of some stress/strain/deformation measure that
     * performes work. Thus, components corresponding to imposed zero stress (e.g. plane
     * stress etc.) are not included. On the other hand, if zero strain components are
     * imposed( e.g. plane strain etc.) this condition must be taken into account in
     * geometrical relations. Therefore, these corresponding components are included in
     * the reduced vector. Which compnents to include are given by the particular MaterialMode.
     * Ex: PlaneStress -> [1 2 6 9]
     *
     * @param answer Returned mask.
     * @param mmode Material response mode.
     * @return The number of components in the corresponding full vector.
     */
    static int giveVoigtVectorMask(IntArray &answer, MaterialMode mmode);

    /**
     * The same as giveVoigtVectorMask but returns a mask corresponding to a symmetric
     * second order tensor.
     *
     * Returns a mask of the vector indicies corresponding to components in a symmetric
     * second order tensor of some stress/strain/deformation measure that performes work.
     * Thus, components corresponding to imposed zero stress (e.g. plane stress etc.) are
     * not included. On the other hand, if zero strain components are imposed( e.g. plane
     * strain etc.) this condition must be taken into account in geometrical relations.
     * Therefore, these corresponding components are included in the reduced vector.
     * Which compnents to include are given by the particular MaterialMode.
     * Ex: PlaneStress -> [1 2 6]
     *
     * @param answer Returned mask.
     * @param mmode Material response mode.
     * @return The number of components in the corresponding full vector.
     */
    static int giveVoigtSymVectorMask(IntArray &answer, MaterialMode mmode);
    /**
     * Gives the inverted version of giveVoigtVectorMask.
     * @deprecated This will eventually be removed. The normal mapping given in giveVoigtVectorMask can be used instead.
     */
    static void giveInvertedVoigtVectorMask(IntArray &answer, MaterialMode mmode);

    /**
     * Returns the size of reduced stress/strain vector according to given mode.
     * @param mmode Material response mode.
     */
    static int giveSizeOfVoigtVector(MaterialMode mmode);
    /**
     * Returns the size of symmetric part of a reduced stress/strain vector according to given mode.
     * @param mmode Material response mode.
     */
    static int giveSizeOfVoigtSymVector(MaterialMode mmode);

    /// Converts the reduced symmetric Voigt vector (2nd order tensor) to full form.
    static void giveFullVectorForm(FloatArray &answer, const FloatArray &strainVector,  MaterialMode matMode);
    /// Converts the reduced deformation gradient Voigt vector (2nd order tensor).
    static void giveFullVectorFormF(FloatArray &answer, const FloatArray &strainVector,  MaterialMode matMode);
    /// Converts the reduced unsymmetric Voigt vector (2nd order tensor) to full form.
    static void giveFullSymVectorForm(FloatArray &answer, const FloatArray &vec, MaterialMode matMode);
    /// Converts the full symmetric Voigt vector (2nd order tensor) to reduced form.
    static void giveReducedVectorForm(FloatArray &answer, const FloatArray &vec, MaterialMode matMode);
    /// Converts the full unsymmetric Voigt vector (2nd order tensor) to reduced form.
    static void giveReducedSymVectorForm(FloatArray &answer, const FloatArray &vec, MaterialMode matMode);

    /// Converts the full unsymmetric Voigt matrix (4th order tensor) to reduced form.
    static void giveFullSymMatrixForm(FloatMatrix &answer, const FloatMatrix &red, MaterialMode matMode);
    /// Converts the full symmetric Voigt matrix (4th order tensor) to reduced form.
    static void giveReducedMatrixForm(FloatMatrix &answer, const FloatMatrix &full, MaterialMode matMode);
    /// Converts the full unsymmetric Voigt matrix (4th order tensor) to reduced form.
    static void giveReducedSymMatrixForm(FloatMatrix &answer, const FloatMatrix &full, MaterialMode matMode);

    /**
     * Method for subtracting from reduced space strain vector its stress-independent parts
     * (caused by temperature, shrinkage, creep and possibly by other phenomena).
     * Calls computeStressIndependentStrainVector to obtain stress
     * independent part of strain.
     * @param answer Computed strain vector.
     * @param gp Integration point.
     * @param reducedStrainVector Reduced strain vector.
     * @param tStep Time step.
     * @param mode Determines value mode.
     */
    void giveStressDependentPartOfStrainVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrainVector,
                                               TimeStep *tStep, ValueModeType mode) const;

    int setIPValue(const FloatArray &value, GaussPoint *gp, InternalStateType type) override;
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    /**
     * Method for computing plane stress stiffness matrix of receiver.
     * Default implementation computes 3d stiffness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to plane stress stiffness using reduce method described above.
     * However, this reduction is quite time consuming and if it is possible,
     * it is recommended to overload this method and provide direct method for computing
     * particular stiffness matrix.
     * @param answer Stiffness matrix.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    //@{
    virtual FloatMatrixF< 3, 3 >givePlaneStressStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatMatrixF< 4, 4 >givePlaneStressStiffMtrx_dPdF(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const;

    virtual void givePlaneStressStiffMtrx_dCde(FloatMatrix &answer,
                                               MatResponseMode mmode, GaussPoint *gp,
                                               TimeStep *tStep);
    //@}

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
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    //@{
    virtual FloatMatrixF< 4, 4 >givePlaneStrainStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatMatrixF< 5, 5 >givePlaneStrainStiffMtrx_dPdF(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const;

    virtual void givePlaneStrainStiffMtrx_dCde(FloatMatrix &answer,
                                               MatResponseMode mmode, GaussPoint *gp,
                                               TimeStep *tStep);
    //@}

    /**
     * Method for computing 1d stiffness matrix of receiver.
     * Default implementation computes 3d stiffness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to 1d stiffness using reduce method described above.
     * However, this reduction is quite time consuming and if it is possible,
     * it is recommended to overload this method and provide direct method for computing
     * particular stiffness matrix.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     * @return Stiffness matrix.
     */
    //@{
    virtual FloatMatrixF< 1, 1 >give1dStressStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatMatrixF< 1, 1 >give1dStressStiffMtrx_dPdF(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const;

    virtual void give1dStressStiffMtrx_dCde(FloatMatrix &answer,
                                            MatResponseMode mmode, GaussPoint *gp,
                                            TimeStep *tStep);
    //@}

    /**
     * Method for computing 2d beam layer stiffness matrix of receiver.
     * Default implementation computes 3d stiffness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to 2d beam layer stiffness using reduce method described above.
     * However, this reduction is quite time consuming and if it is possible,
     * it is recommended to overload this method and provide direct method for computing
     * particular stiffness matrix.
     * @param answer Stiffness matrix.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual FloatMatrixF< 2, 2 >give2dBeamLayerStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const;
    /**
     * Method for computing 2d plate layer stiffness matrix of receiver.
     * Default implementation computes 3d stiffness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to 2d plate layer stiffness using reduce method described above.
     * However, this reduction is quite time consuming and if it is possible,
     * it is recommended to overload this method and provide direct method for computing
     * particular stiffness matrix.
     * @param answer Stiffness matrix.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual FloatMatrixF< 5, 5 >givePlateLayerStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const;
    /**
     * Method for computing 1d fiber stiffness matrix of receiver.
     * Default implementation computes 3d stiffness matrix using give3dMaterialStiffnessMatrix and
     * reduces it to 1d fiber stiffness using reduce method described above.
     * However, this reduction is quite time consuming and if it is possible,
     * it is recommended to overload this method and provide direct method for computing
     * particular stiffness matrix.
     * @param answer Stiffness matrix.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual FloatMatrixF< 3, 3 >giveFiberStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const;

    /**
     * Method for computing stiffness matrix of plate subsoil model.
     * Default method is emty; the implementation should be provided by the particular model.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     * @return Stiffness matrix.
     */
    virtual FloatMatrixF< 3, 3 >give2dPlateSubSoilStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const;
    /**
     * Method for computing stiffness matrix of beam3d subsoil model.
     * Default method is emty; the implementation should be provided by the particular model.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     * @return Stiffness matrix.
     */
    virtual FloatMatrixF< 6, 6 >give3dBeamSubSoilStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const;
    /**
     * Transforms 3d strain vector into another coordinate system.
     * @param answer Transformed strain vector
     * @param base Transformation matrix. There are on each column stored unit vectors of
     * coordinate system (so called base vectors) to which we do transformation. These vectors must
     * be expressed in the same coordinate system as source strainVector.
     * @param strainVector 3d strain.
     * @param transpose Determines if we transpose matrix before transforming.
     */
    static FloatArrayF< 6 >transformStrainVectorTo(const FloatMatrixF< 3, 3 > &base,
                                                   const FloatArrayF< 6 > &strain, bool transpose = false);
    /**
     * Transforms 3d stress vector into another coordinate system.
     * @param answer Transformed stress vector.
     * @param base Transformation matrix. There are on each column stored unit vectors of
     * coordinate system (so called base vectors) to which we do transformation. These vectors must
     * be expressed in the same coordinate system as source stressVector.
     * @param stressVector Transformed 3d strain.
     * @param transpose Determines if we transpose matrix before transforming.
     */
    static FloatArrayF< 6 >transformStressVectorTo(const FloatMatrixF< 3, 3 > &base,
                                                   const FloatArrayF< 6 > &stress, bool transpose = false);
    /**
     * Computes equivalent of von Mises stress. Returns 0 if six stress components do not exist on the material.
     * @param currentStress Stress vector given by 6 components.
     */
    static double computeVonMisesStress(const FloatArray &currentStress);
    static double computeVonMisesStress_3D(const FloatArrayF< 6 > &stress);
    static double computeVonMisesStress_PlaneStress(const FloatArrayF< 3 > &stress);

    /**
     * Computes 3d strain vector transformation matrix from standard vector transformation matrix.
     * @param base A (3,3) matrix, where on each column are stored unit direction vectors of
     * local coordinate axes to which we do transformation.
     * @param transpose Determines if we transpose matrix before transforming.
     * @returns Transformation matrix for strain vector.
     */
    static FloatMatrixF< 6, 6 >giveStrainVectorTranformationMtrx(const FloatMatrixF< 3, 3 > &base, bool transpose = false);

    /**
     * Computes 2d strain vector transformation matrix from standard vector transformation matrix.
     * @param base A (2,2) matrix, where on each column are stored unit direction vectors of
     * local coordinate axes to which we do transformation.
     * @param transpose Determines if we transpose matrix before transforming.
     * @returns Transformation matrix for strain vector.
     */
    static FloatMatrixF< 3, 3 >give2DStrainVectorTranformationMtrx(const FloatMatrixF< 2, 2 > &base, bool transpose = false);

    /**
     * Computes 3d stress vector transformation matrix from standard vector transformation matrix.
     * @param base A (3,3) matrix, where on each column are stored unit direction vectors of
     * local coordinate axes to which we do transformation.
     * @param transpose Determines if we transpose matrix before transforming.
     * @returns Transformation matrix for stress vector.
     */
    static FloatMatrixF< 6, 6 >giveStressVectorTranformationMtrx(const FloatMatrixF< 3, 3 > &base, bool transpose = false);
    /**
     * Computes 2d stress vector transformation matrix from standard vector transformation matrix.
     * @param base A (2,2) matrix, where on each column are stored unit direction vectors of
     * local coordinate axes to which we do transformation.
     * @param transpose Determines if we transpose matrix before transforming.
     * @return Transformation matrix for stress vector.
     */
    static FloatMatrixF< 3, 3 >givePlaneStressVectorTranformationMtrx(const FloatMatrixF< 2, 2 > &base, bool transpose = false);
    /**
     * Method for sorting newly computed principal values (pVal) and
     * corresponding principal directions (pDir) to be closed
     * to some (often previous) principal directions (toPDir).
     * pDir and toPDir should have eigenvectors stored in columns and normalized.
     * @param pVal New eigenvalues.
     * @param pDir New eigenvectors.
     * @param toPDir Old eigenvector.
     */
    static void sortPrincDirAndValCloseTo(FloatArray &pVal, FloatMatrix &pDir, const FloatMatrix &toPDir);

    friend class CrossSection;
    friend class StructuralCrossSection;
    friend class SimpleCrossSection;
    friend class LayeredCrossSection;
};
} // end namespace oofem
#endif // structuralmaterial_h
