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

#ifndef latticestructuralmaterial_h
#define latticestructuralmaterial_h

#include "../structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matconst.h"
#include "matstatus.h"
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
 * Structural material services should not be called directly by elements. Instead, they always should
 * pass their requests to corresponding cross section model. Cross section performs all necessary integration over
 * its volume and invokes material model services.
 *
 * @author Peter Grassl
 */
class LatticeStructuralMaterial : public StructuralMaterial
{
protected:
    /// Reference temperature (temperature, when material has been built into structure).
    double referenceTemperature = 0.;

public:

    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    LatticeStructuralMaterial(int n, Domain *d);

    void giveStiffnessMatrix(FloatMatrix &answer,
                             MatResponseMode mode,
                             GaussPoint *gp,
                             TimeStep *tStep) override;

    virtual bool hasAnalyticalTangentStiffness() const { return true; }

    bool hasMaterialModeCapability(MaterialMode mode) const override;
    const char *giveClassName() const override { return "LatticeStructuralMaterial"; }

    virtual double giveLatticeStress1d(const double strain, GaussPoint *gp, TimeStep *tStep);

    virtual FloatArrayF< 3 >giveLatticeStress2d(const FloatArrayF< 3 > &strain, GaussPoint *gp, TimeStep *tStep);

    virtual FloatArrayF< 6 >giveLatticeStress3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep);

    virtual FloatArrayF< 6 >giveFrameForces3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep);

    virtual FloatMatrixF< 1, 1 >give1dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatMatrixF< 3, 3 >give2dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatMatrixF< 6, 6 >give3dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const;

    virtual FloatMatrixF< 6, 6 >give3dFrameStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime) override;
};
} // end namespace oofem
#endif // latticestructuralmaterial_h
