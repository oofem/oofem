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

#ifndef structuralinterfacematerialstatus_h
#define structuralinterfacematerialstatus_h

#include "matstatus.h"
#include "floatarray.h"
#include "matstatmapperint.h"
#include "floatmatrix.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"

namespace oofem {
class GaussPoint;

/**
 * This class implements a structural interface material status information. It is attribute of
 * gaussPoint. This is only an abstract class, for every instance of this material class
 * there should be a specialized derived class, which handles are history variables.
 *
 * This is a base class for all material statuses corresponding to materials derived from
 * StructuralInterfaceMaterial class.
 * It defines the traction, jump (discontinuity), deformation gradient and their temporary counterparts.
 * Functions for accessing these components are defined.
 */
class StructuralInterfaceMaterialStatus : public MaterialStatus, public MaterialStatusMapperInterface
{
protected:
    /// Equilibrated jump (discontinuity)
    FloatArrayF<3> jump;
    /// Equilibrated (engineering) traction vector
    FloatArrayF<3> traction;
    /// Temporary (engineering) traction vector
    FloatArrayF<3> tempTraction;
    /// Temporary jump (discontinuity)
    FloatArrayF<3> tempJump;

    /// Equilibrated first Piola-Kirchhoff traction vector T
    FloatArrayF<3> firstPKTraction;
    /// Temporary first Piola-Kirchhoff traction vector (to find balanced state)
    FloatArrayF<3> tempFirstPKTraction;
    /// Equilibrated deformation gradient in reduced form
    FloatMatrixF<3,3> F;
    /// Temporary deformation gradient in reduced form (to find balanced state)
    FloatMatrixF<3,3> tempF;

    /// Interface normal direction
    FloatArrayF<3> mNormalDir;

    bool mNewlyInserted = true;

    FloatArrayF<2> projectedTraction;

public:
    /// Constructor. Creates new StructuralInterfaceMaterialStatus with number n, belonging to domain d and IntegrationPoint g.
    StructuralInterfaceMaterialStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    /// Returns the const pointer to receiver's jump.
    const FloatArrayF<3> &giveJump() const { return jump; }
    /// Returns the const pointer to receiver's traction vector.
    const FloatArrayF<3> &giveTraction() const { return traction; }
    /// Returns the const pointer to receiver's first Piola-Kirchhoff traction vector.
    const FloatArrayF<3> &giveFirstPKTraction() const { return firstPKTraction; }
    /// Returns the const pointer to receiver's deformation gradient vector.
    const FloatMatrixF<3,3> &giveF() const { return F; }
    /// Returns the const pointer to receiver's temporary jump.
    const FloatArrayF<3> &giveTempJump() const { return tempJump; }
    /// Returns the const pointer to receiver's temporary traction vector.
    const FloatArrayF<3> &giveTempTraction() const { return tempTraction; }
    /// Returns the const pointer to receiver's temporary first Piola-Kirchhoff traction vector.
    const FloatArrayF<3> &giveTempFirstPKTraction() const { return tempFirstPKTraction; }
    /// Returns the const pointer to receiver's temporary deformation gradient vector.
    const FloatMatrixF<3,3> &giveTempF() const { return tempF; }
    /// Returns const reference to normal vector.
    const FloatArrayF<3> &giveNormal() const { return mNormalDir; }
    /// Returns the projected traction.
    const FloatArrayF<2> &giveProjectedTraction() const { return projectedTraction; }
    /// Assigns jump to given vector v.
    void letJumpBe(const FloatArrayF<3> v) { jump = v; }
    /// Assigns traction to given vector v.
    void letTractionBe(const FloatArrayF<3> v) { traction = v; }
    /// Assigns firstPKTraction to given vector v.
    void letFirstPKTractionBe(const FloatArrayF<3> v) { firstPKTraction = v; }
    /// Assigns FVector to given vector v.
    void letFBe(const FloatMatrixF<3,3> v) { F = v; }
    /// Assigns tempTraction to given vector v.
    void letTempTractionBe(const FloatArrayF<3> v) { tempTraction = v; }
    /// Assigns tempJump to given vector v
    void letTempJumpBe(const FloatArrayF<3> v) { tempJump = v; }
    /// Assigns tempFirstPKTraction to given vector v
    void letTempFirstPKTractionBe(const FloatArrayF<3> v) { tempFirstPKTraction = v; }
    /// Assigns tempFVector to given vector v
    void letTempFBe(const FloatMatrixF<3,3> &v) { tempF = v; }
    /// Assigns normal vector
    void letNormalBe(const FloatArrayF<3> &iN) { mNormalDir = iN; }
    /// Assigns projeted traction
    void letProjectedTractionBe(const FloatArrayF<2> &iProjectedTraction) { projectedTraction = iProjectedTraction; }
    ///@TODO Projected tractions are *never* set. What is it supposed to be good for?

    const char *giveClassName() const override { return "StructuralInterfaceMaterialStatus"; }

    /// Functions for MaterialStatusMapperInterface
    void copyStateVariables(const MaterialStatus &iStatus) override;
    void addStateVariables(const MaterialStatus &iStatus) override;

    bool giveNewlyInserted() const { return mNewlyInserted; }
    void setNewlyInserted(bool iNewlyInserted) { mNewlyInserted = iNewlyInserted; }

    virtual double giveDamage() const { return 0.0; }     // no default damage
    virtual double giveTempDamage() const { return 0.0; } // no default damage
};
} // end namespace oofem
#endif // structuralinterfacematerialstatus_h
