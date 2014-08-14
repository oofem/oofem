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

#ifndef simplevitrificationmaterial_h
#define simplevitrificationmaterial_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "floatarray.h"

///@name Input fields for SimpleVitrificationMaterial
//@{
#define _IFT_SimpleVitrificationMaterial_Name "simplevitrificationmaterial"
#define _IFT_SimpleVitrificationMaterial_vitrificationTime "vitrificationtime" ///< Describes the time where the material switches response.
#define _IFT_SimpleVitrificationMaterial_E "e"
#define _IFT_SimpleVitrificationMaterial_nu "nu"
#define _IFT_SimpleVitrificationMaterial_G "g"
#define _IFT_SimpleVitrificationMaterial_alpha "alpha"
#define _IFT_SimpleVitrificationMaterial_E_r "e_r"
#define _IFT_SimpleVitrificationMaterial_nu_r "nu_r"
#define _IFT_SimpleVitrificationMaterial_G_r "g_r"
#define _IFT_SimpleVitrificationMaterial_alpha_r "alpha_r"
//@}

namespace oofem {
/**
 * Model describing the vitrification process of a glass like material.
 * The model is characterized by a linear anisotropic response followed by switch in tangent stiffness as the material vitrifies.
 *
 * @author Tomasz Garstka
 * @author Mikael Öhman
 */
class SimpleVitrificationMaterial : public StructuralMaterial
{
private:
    /// Vitrification time (when equal or larger than this time the material changes response).
    double vitrTime;
    /// Material parameters for the glassy part of the model (after vitrification).
    FloatArray E, nu, G, alpha;
    /// Material parameters for the rubbery part of the model (before vitrification).
    FloatArray E_r, nu_r, G_r, alpha_r;

public:
    /// Constructor.
    SimpleVitrificationMaterial(int n, Domain * d) : StructuralMaterial(n, d) { }
    /// Destructor.
    virtual ~SimpleVitrificationMaterial();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual int checkConsistency();

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                         const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual int hasNonLinearBehaviour() { return true; }
    virtual const char *giveClassName() const { return "SimpleVitrificationMaterial"; }
    virtual const char *giveInputRecordName() const { return _IFT_SimpleVitrificationMaterial_Name; }
};
} // end namespace oofem
#endif // simplevitrificationmaterial_h
