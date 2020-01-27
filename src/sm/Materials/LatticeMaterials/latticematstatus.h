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

#ifndef latticematstatus_h
#define latticematstatus_h


#include "../structuralms.h"
#include "randommaterialext.h"
#include "floatarrayf.h"

namespace oofem {
class GaussPoint;
class Dictionary;
class Domain;
class NonlocalMaterialStatusExtension;

/**
 * This class implements a base lattice material status.
 * In this class services are defined that are used by other
 * lattice material statuses.
 */
class LatticeMaterialStatus : public MaterialStatus, public RandomMaterialStatusExtensionInterface
{
protected:

    /// Equilibriated lattice strain
    FloatArrayF< 6 >latticeStrain;

    /// Non-equilibriated lattice strain
    FloatArrayF< 6 >tempLatticeStrain;

    /// Equilibriated lattice stress
    FloatArrayF< 6 >latticeStress;

    /// Non-equilibriated lattice stress
    FloatArrayF< 6 >tempLatticeStress;

    /// Equilibriated reduced lattice strain, which is free of thermal strain
    FloatArrayF< 6 >reducedLatticeStrain;
    /// Non-equilibrated reduced lattice strain, which is free of thermal strain
    FloatArrayF< 6 >tempReducedLatticeStrain;
    /// Equilibriated plastic lattice strain
    FloatArrayF< 6 >plasticLatticeStrain;
    /// Non-equilibrated plastic lattice strain
    FloatArrayF< 6 >tempPlasticLatticeStrain;
    /// Non-equilibrated plastic lattice strain
    FloatArrayF< 6 >oldPlasticLatticeStrain;


    /// Equilibriated damage lattice strain
    FloatArrayF< 6 >damageLatticeStrain;

    /// Non-equilibriated damage lattice strain
    FloatArrayF< 6 >tempDamageLatticeStrain;


    /// Equilibrated normal stress
    double normalLatticeStress = 0.;

    /// Non-equilibrated normal stress
    double tempNormalLatticeStress = 0.;

    /// dissipation
    double dissipation = 0.;

    /// Non-equilibrated increment of dissipation
    double tempDissipation = 0.;

    ///Increment of dissipation
    double deltaDissipation = 0.;

    /// Non-equilibrated increment of dissipation
    double tempDeltaDissipation = 0.;

    /// Characteristic length
    double le = 0.;

    /** the crack_flag indicates if the gp is damaged (cracked):
     *  crack_flag = 0 gp is undamaged
     *  crack_flag = 1 gp is damaged and damage grows
     *  crack_flag = 2 gp is damaged and damage does not grow
     */
    int crackFlag = 0;

    /// Non-equilibrated temp flag.
    int tempCrackFlag = 0;

    /// Non-equilibrated crack width
    double tempCrackWidth = 0.;

    /// Crack width
    double crackWidth = 0.;


    int updateFlag = 0;

public:
    LatticeMaterialStatus(GaussPoint *g);

    const char *giveClassName() const override { return "LatticeMaterialStatus"; }

    void initTempStatus() override;

    void updateYourself(TimeStep *) override;

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    /// Returns lattice strain.
    const FloatArrayF< 6 > &giveLatticeStrain() const { return this->latticeStrain; }
    /// Returns lattice strain.
    const FloatArrayF< 6 > &giveTempLatticeStrain() const { return this->tempLatticeStress; }

    /// Returns reduced lattice strain.
    const FloatArrayF< 6 > &giveReducedLatticeStrain() const { return reducedLatticeStrain; }

    /// Returns temp reduced lattice strain.
    const FloatArrayF< 6 > &giveTempReducedLatticeStrain() const { return tempReducedLatticeStrain; }

    /// Returns plastic lattice strain.
    const FloatArrayF< 6 > &givePlasticLatticeStrain() const { return this->plasticLatticeStrain; }

    /// Returns temp plastic lattice strain.
    const FloatArrayF< 6 > &giveTempPlasticLatticeStrain() const { return this->tempPlasticLatticeStrain; }

    /// Returns plastic lattice strain.
    const FloatArrayF< 6 > &giveOldPlasticLatticeStrain() const { return this->oldPlasticLatticeStrain; }

    /// Returns lattice stress.
    const FloatArrayF< 6 > &giveLatticeStress() const { return this->latticeStress; }
    /// Returns temp lattice stress.
    const FloatArrayF< 6 > &giveTempLatticeStress() const { return this->tempLatticeStress; }

    /// Returns temp damage lattice strain.
    const FloatArrayF< 6 > &giveTempDamageLatticeStrain() const { return this->tempDamageLatticeStrain; }

    /// Assigns the temp value of lattice strain.
    void letTempLatticeStrainBe(const FloatArrayF< 6 > &v) { this->tempLatticeStrain = v; }

    /// Assigns the temp value of lattice strain.
    void letTempReducedLatticeStrainBe(const FloatArrayF< 6 > &v) { this->tempReducedLatticeStrain = v; }

    /// Assigns the temp value of lattice strain.
    void letTempPlasticLatticeStrainBe(const FloatArrayF< 6 > &v) { this->tempPlasticLatticeStrain = v; }

    /// Assigns the temp value of lattice stress.
    void letTempLatticeStressBe(const FloatArrayF< 6 > &v) { this->tempLatticeStress = v; }

    /// Assigns the temp value of damage lattice strain.
    void letTempDamageLatticeStrainBe(const FloatArrayF< 6 > &v) { this->tempDamageLatticeStrain = v; }


    /// Sets the temp normalStress
    void setTempNormalLatticeStress(double val) { this->tempNormalLatticeStress = val; }

    /// Gives the last equilibrated normal stress
    double giveNormalLatticeStress() const { return this->normalLatticeStress; }

    /// Gives the last equilibrated normal stress
    double giveTempNormalLatticeStress() const { return this->tempNormalLatticeStress; }


    ///Sets the temp_crack_flag
    void setTempCrackFlag(int val) { tempCrackFlag = val; }

    ///Sets the temp_crack_width
    void setTempCrackWidth(double val) { tempCrackWidth = val; }

    /**
     * Returns the crack flag
     * @return crack flag
     */
    virtual int giveCrackFlag() const { return this->crackFlag; }

    /**
     * @return crack width
     */
    virtual double giveCrackWidth() const { return this->crackWidth; }


    /// Returns characteristic length stored in receiver
    double giveLe() const { return le; }

    /// Sets characteristic length to given value
    void setLe(double ls) { le = ls; }

    virtual int hasBeenUpdated() const { return this->updateFlag; }

    /**
     * Returns the energy dissipation computed at the GaussPoint of the element.
     * This function is used for the lattice specific vtk export.
     * @return dissipation
     */
    virtual double giveDissipation() const { return dissipation; }
    double giveTempDissipation() const { return tempDissipation; }
    void setTempDissipation(double newDiss) { tempDissipation = newDiss; }

    /**
     * Returns the increment of dissipation computed at the GaussPoint of the element.
     * This function is used for the lattice specific vtk export.
     * @return increment of dissipation
     */
    virtual double giveDeltaDissipation() const { return deltaDissipation; }
    double giveTempDeltaDissipation() const { return tempDeltaDissipation; }
    void setTempDeltaDissipation(double newDiss) { tempDeltaDissipation = newDiss; }

    Interface *giveInterface(InterfaceType) override;

    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;
};
} // end namespace oofem
#endif // matstatus_h
