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

#ifndef isointerfacedamage02_h
#define isointerfacedamage02_h

#include "../structuralinterfacematerial.h"
#include "../structuralinterfacematerialstatus.h"

#include <fstream>

///@name Input fields for IsoInterfaceDamageMaterial
//@{
#define _IFT_IsoInterfaceDamageMaterial_2_Name "isointrfdm02"
#define _IFT_IsoInterfaceDamageMaterial_2_tablename "tablename"
#define _IFT_IsoInterfaceDamageMaterial_2_kn "kn"
#define _IFT_IsoInterfaceDamageMaterial_2_ks "ks"
#define _IFT_IsoInterfaceDamageMaterial_2_ft "ft"
#define _IFT_IsoInterfaceDamageMaterial_2_maxOmega "maxomega"
//@}

namespace oofem {

/**
 * This class implements associated Material Status to IsoInterfaceDamageMaterial_2.
 * @author Kristoffer Carlsson
 */
class IsoInterfaceDamageMaterialStatus_2 : public StructuralInterfaceMaterialStatus
{
protected:
    /// Scalar measure of the largest equivalent displacement ever reached in material.
    double kappa;
    /// Non-equilibrated scalar measure of the largest equivalent displacement.
    double tempKappa;
    /// Damage level of material.
    double damage;
    /// Non-equilibrated damage level of material.
    double tempDamage;

public:
    /// Constructor
    IsoInterfaceDamageMaterialStatus_2(int n, Domain * d, GaussPoint * g);
    /// Destructor
    virtual ~IsoInterfaceDamageMaterialStatus_2();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    /// Returns the last equilibrated scalar measure of the largest strain level.
    double giveKappa() { return kappa; }
    /// Returns the temp. scalar measure of the largest strain level.
    double giveTempKappa() { return tempKappa; }
    /// Sets the temp scalar measure of the largest strain level to given value.
    void setTempKappa(double newKappa) { tempKappa = newKappa; }
    /// Returns the last equilibrated damage level.
    double giveDamage() { return damage; }
    /// Returns the temp. damage level.
    double giveTempDamage() { return tempDamage; }
    /// Sets the temp damage level to given value.
    void setTempDamage(double newDamage) { tempDamage = newDamage; }

    // definition
    virtual const char *giveClassName() const { return "IsoInterfaceDamageMaterialStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};


/**
 * Simple isotropic damage based model for 2d and 3d interface elements.
 * In 2d, the interface elements are used to model contact layer between
 * element edges. The generalized strain vector contains two relative displacements
 * (in normal and shear direction). The generalized stress vector contains corresponding
 * tractions in normal and tangent direction.
 *
 * In 3d, the interface elements are used to model contact layer between
 * element surfaces. The generalized strain vector contains two relative displacements
 * (in normal and shear direction). The generalized stress vector contains corresponding
 * tractions in normal and tangent direction
 *
 * The behaviour of the model is elastic, described by normal and shear stiffness components.
 * Isotropic damage is initiated  when the stress reaches the tensile strength. Damage evolution
 * is governed by normal component of generalized strain vector (normal relative displacement)
 * by a table given by a file that relates the normal displacement to the damage. A linear interpolation
 * is made between the values given in the table. If the strain is greater than the largest value
 * in the table the largest damage in the table will be used.
 *
 * Differences between this class and IsoInterfaceDamageMaterial written by:
 * @author Kristoffer Carlsson
 */
class IsoInterfaceDamageMaterial_2 : public StructuralInterfaceMaterial
{
protected:
    /// Elastic properties (normal moduli).
    double kn;
    /// Shear moduli.
    double ks;
    /// Tension strength.
    double ft;
    /// Limit elastic deformation.
    double e0;
    /// Maximum limit on omega. The purpose is elimination of a too compliant material which may cause convergency problems. Set to something like 0.99 if needed.
    double maxOmega;
    /// Name of table file
    std :: string tablename;
    /// Damages read from the second column in the table file
    FloatArray damages;
    /// Strains read from the first column in the table file
    FloatArray strains;

public:
    /// Constructor
    IsoInterfaceDamageMaterial_2(int n, Domain * d);
    /// Destructor
    virtual ~IsoInterfaceDamageMaterial_2();

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual bool hasAnalyticalTangentStiffness() const { return true; }

    virtual const char *giveInputRecordName() const { return _IFT_IsoInterfaceDamageMaterial_2_Name; }
    virtual const char *giveClassName() const { return "IsoInterfaceDamageMaterial"; }

    virtual void giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep);
    virtual void give3dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    /**
     * Computes the equivalent strain measure from given strain vector (full form).
     * @param[out] kappa Return parameter containing the corresponding equivalent strain.
     * @param strain Total strain vector in full form.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);

    /**
     * computes the value of damage parameter omega, based on given value of equivalent strain.
     * @param[out] omega Contains result.
     * @param kappa Equivalent strain measure.
     * @param strain Total strain vector in full form. (unnecessary?)
     * @param gp Integration point.
     */
    virtual void computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new IsoInterfaceDamageMaterialStatus_2(1, domain, gp); }
};
} // end namespace oofem
#endif // isointerfacedamage01_h
