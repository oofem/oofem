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

#ifndef bondceb_h
#define bondceb_h

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for BondCEBMaterial
//@{
#define _IFT_BondCEBMaterial_Name "bondceb"
#define _IFT_BondCEBMaterial_kn "kn"
#define _IFT_BondCEBMaterial_ks "ks"
#define _IFT_BondCEBMaterial_s1 "s1"
#define _IFT_BondCEBMaterial_s2 "s2"
#define _IFT_BondCEBMaterial_s3 "s3"
#define _IFT_BondCEBMaterial_al "al"
#define _IFT_BondCEBMaterial_taumax "taumax"
#define _IFT_BondCEBMaterial_tauf "tauf"
//@}

namespace oofem {

/**
 * This class implements associated status to BondCEBInterfaceMaterial.
 */
class BondCEBMaterialStatus : public StructuralInterfaceMaterialStatus
{
protected:
    /// Cumulative slip.
    double kappa;
    /// Non-equilibrated cumulative slip.
    double tempKappa;

public:
    /// Constructor
    BondCEBMaterialStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor
    virtual ~BondCEBMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    /// Returns the last equilibrated cumulative slip.
    double giveKappa() { return kappa; }
    /// Returns the temporary cumulative slip.
    double giveTempKappa() { return tempKappa; }
    /// Sets the temporary cumulative slip to the given value.
    void setTempKappa(double newKappa) { tempKappa = newKappa; }

    // definition
    virtual const char *giveClassName() const { return "BondCEBMaterialStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
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
class BondCEBMaterial : public StructuralInterfaceMaterial
{
protected:
    /// Normal elastic stiffness.
    double kn;
    /// Shear elastic stiffness.
    double ks;
    /// Shear strength.
    double taumax;
    /// Residual shear stress.
    double tauf;
    /// Characteristic slip values.
    double s0, s1, s2, s3;
    /// Exponent.
    double alpha;

public:
    /// Constructor
    BondCEBMaterial(int n, Domain * d);
    /// Destructor
    virtual ~BondCEBMaterial();

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual bool hasAnalyticalTangentStiffness() const { return true; }

    virtual const char *giveInputRecordName() const { return _IFT_BondCEBMaterial_Name; }
    virtual const char *giveClassName() const { return "BondCEBMaterial"; }

    virtual void giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep);
    virtual void give3dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new BondCEBMaterialStatus(1, domain, gp); }

protected:
    double evaluateBondStress(const double kappa);
 };
} // end namespace oofem
#endif // bondceb_h
