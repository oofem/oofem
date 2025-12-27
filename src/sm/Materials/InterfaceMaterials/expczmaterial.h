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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef expczmaterial_h
#define expczmaterial_h

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for ExpCZMaterial
//@{
#define _IFT_ExpCZMaterial_kn "kn"
#define _IFT_ExpCZMaterial_ks "ks"
#define _IFT_ExpCZMaterial_knc "knc"
#define _IFT_ExpCZMaterial_g1c "g1c"
#define _IFT_ExpCZMaterial_sigfn "sigfn"
#define _IFT_ExpCZMaterial_sigfs "sigfs"
//@}

namespace oofem {

/**
 * This class implements associated status to ExpCZMaterial.
 */
class ExpCZMaterialStatus : public StructuralInterfaceMaterialStatus
{
protected:

public:
    /// Constructor
    ExpCZMaterialStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) override;

    const char *giveClassName() const override { return "ExpCZMaterialStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    //void saveContext(DataStream &stream, ContextMode mode) override;
    //void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * Simple isotropic damage based model for 2d interface elements.
 * In 2d, the interface elements are used to model contact layer between
 * element edges. The generalized strain vector contains two relative displacements
 * (in normal and shear direction). The generalized stress vector contains corresponding
 * tractions in normal and tangent direction.
 *
 * The behaviour of the model is elastic, described by normal and shear stiffness components.
 * Isotropic damage is initiated  when the stress reaches the tensile strength. Damage evolution
 * is governed by normal component of generalized strain vector (normal relative displacement)
 * by an exponential softening law.
 */
class ExpCZMaterial : public StructuralInterfaceMaterial
{
protected:
    /// Material parameters
    double kn0 = 0.;
    double ks0 = 0.;
    double GIc = 0.;
    double GIIc = 0.;
    double sigfn = 0.;
    double sigfs = 0.;

    /// normal jump at damage initiation
    double gn0 = 0.;
    /// shear jump at damage initiations
    double gs0 = 0.;
    double q = 0.;
    double r = 0.;

public:
    /// Constructor
    ExpCZMaterial(int n, Domain * d);

    int checkConsistency() override;

    const char *giveClassName() const override { return "ExpCZMaterial"; }

    void giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep) override;
    void give3dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    void initializeFrom(InputRecord &ir) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new ExpCZMaterialStatus(gp); }
    void printYourself() override;
};
} // end namespace oofem
#endif // isointerfacedamage01_h
