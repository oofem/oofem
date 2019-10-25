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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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
#include "matconst.h"
#include "matstatus.h"
#include "../stressstrainprincmode.h"
#include "valuemodetype.h"
#include <vector>

///@name Input fields for LatticeStructuralMaterial
//@{
#define _IFT_LatticeStructuralMaterial_referencetemperature "referencetemperature"
#define _IFT_LatticeStructuralMaterial_talpha "talpha"
//@}

namespace oofem {
#define STRAIN_STEPS 10.0

class GaussPoint;
///@todo Update documentation
/**
 * Abstract base class for all lattice "structural" constitutive models. It declares common  services provided
 * by all lattice structural material models. The implementation of these services is partly left on derived classes,
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
class LatticeStructuralMaterial : public Material
{
protected:
    /// Reference temperature (temperature, when material has been built into structure).
    double referenceTemperature;

public:
    /// Voigt index map
    static std :: vector< std :: vector< int > >vIindex;

    /// Symmetric Voigt index map
    static std :: vector< std :: vector< int > >svIndex;

    static int giveSymVI(int ind1, int ind2) { return svIndex [ ind1 - 1 ] [ ind2 - 1 ]; }
    static int giveVI(int ind1, int ind2) { return vIindex [ ind1 - 1 ] [ ind2 - 1 ]; }

    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    LatticeStructuralMaterial(int n, Domain *d);

    bool hasMaterialModeCapability(MaterialMode mode) const override;
    const char *giveClassName() const override { return "LatticeStructuralMaterial"; }

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    virtual void giveRealStressVector_Lattice1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep) {};
    virtual void giveRealStressVector_Lattice2d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep);
    virtual void giveRealStressVector_Lattice3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep);

    /**
     * Method for computing 2d lattice stiffness matrix of receiver.
     * @param answer Stiffness matrix.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void give1dLatticeStiffMtrx(FloatMatrix &answer,
                                        MatResponseMode mmode, GaussPoint *gp,
                                        TimeStep *tStep);

    /**
     * Method for computing 2d lattice stiffness matrix of receiver.
     * @param answer Stiffness matrix.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void give2dLatticeStiffMtrx(FloatMatrix &answer,
                                        MatResponseMode mmode, GaussPoint *gp,
                                        TimeStep *tStep);


    /**
     * Method for computing 3d lattice stiffness matrix of receiver.
     * @param answer Stiffness matrix.
     * @param mmode Material response mode.
     * @param gp Integration point, which load history is used.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void give3dLatticeStiffMtrx(FloatMatrix &answer,
                                        MatResponseMode mmode, GaussPoint *gp,
                                        TimeStep *tStep);

};
} // end namespace oofem
#endif // latticestructuralmaterial_h
