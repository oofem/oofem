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


#include "structuralms.h"
#include "randommaterialext.h"

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
class LatticeMaterialStatus : public StructuralMaterialStatus, public RandomMaterialStatusExtensionInterface
{
protected:

    /// Equilibrated normal stress
    double normalStress;

    /// Non-equilibrated normal stress
    double tempNormalStress;

    /// Reduced strain, which is temperature free
    FloatArray reducedStrain;

    /// Non-equilibrated reduced strain, which is temperature free
    FloatArray tempReducedStrain;

    /// dissipation
    double dissipation;

    /// Non-equilibrated increment of dissipation
    double tempDissipation;

    ///Increment of dissipation
    double deltaDissipation;

    /// Non-equilibrated increment of dissipation
    double tempDeltaDissipation;

    FloatArray plasticStrain;

    FloatArray tempPlasticStrain;

    FloatArray oldPlasticStrain;

    /// Characteristic length
    double le;

    /** the crack_flag indicates if the gp is damaged (cracked):
     *  crack_flag = 0 gp is undamaged
     *  crack_flag = 1 gp is damaged and damage grows
     *  crack_flag = 2 gp is damaged and damage does not grow
     */
    int crackFlag;

    /// Non-equilibrated temp flag.
    int tempCrackFlag;

    /// Non-equilibrated crack width
    double tempCrackWidth;

    /// Crack width
    double crackWidth;


    int updateFlag;


public:

    /// Constructor
    LatticeMaterialStatus(GaussPoint *g);
    /// Destructor.
    virtual ~LatticeMaterialStatus() { }

    const char *giveClassName() const override { return "LatticeMaterialStatus"; }

    void initTempStatus() override;

    void updateYourself(TimeStep *) override;

    void printOutputAt(FILE *file, TimeStep *tStep) override;


    /*
     * Assign the temp value of reduced strain.
     * @v new temp value of reduced strain
     */
    virtual void  letTempReducedStrainBe(const FloatArray &v) { tempReducedStrain = v; }

    /// Gives the temp value of reduced strain.
    virtual const FloatArray &giveTempReducedStrain() const { return tempReducedStrain; }

    /// Gives the old equilibrated value of reduced strain.
    virtual const FloatArray &giveReducedStrain() const { return reducedStrain; }



    /// Sets the temp normalStress
    virtual void setTempNormalStress(double val) { tempNormalStress = val; }

    /// Gives the last equilibrated normal stress
    virtual double giveNormalStress() { return tempNormalStress; }

    /// Gives the last equilibrated normal stress
    virtual double giveOldNormalStress() { return normalStress; }

    ///Sets the temp_crack_flag
    virtual void setTempCrackFlag(int val) { tempCrackFlag = val; }

    ///Sets the temp_crack_width
    void setTempCrackWidth(double val) { tempCrackWidth = val; }

    virtual FloatArray &givePlasticStrain()
    { return this->plasticStrain; }

    virtual FloatArray &giveTempPlasticStrain()
    { return this->tempPlasticStrain; }

    virtual FloatArray &giveOldPlasticStrain()
    { return this->oldPlasticStrain; }


    /**
     * Returns the crack flag
     * @return crack flag
     */
    virtual int giveCrackFlag() { return this->crackFlag; }

    /**
     * @return crack width
     */
    virtual double giveCrackWidth() { return this->crackWidth; }


    /// Returns characteristic length stored in receiver
    double giveLe()  { return le; }

    /// Sets characteristic length to given value
    void   setLe(double ls) { le = ls; }

    virtual int hasBeenUpdated() { return this->updateFlag; }

    /**
     * Returns the energy dissipation computed at the GaussPoint of the element.
     * This function is used for the lattice specific vtk export.
     * @return dissipation
     */
    virtual double giveDissipation() { return dissipation; }
    double giveTempDissipation() { return tempDissipation; }
    void setTempDissipation(double newDiss) { tempDissipation = newDiss; }

    /**
     * Returns the increment of dissipation computed at the GaussPoint of the element.
     * This function is used for the lattice specific vtk export.
     * @return increment of dissipation
     */
    virtual double giveDeltaDissipation() { return deltaDissipation; }
    double giveTempDeltaDissipation() { return tempDeltaDissipation; }
    void setTempDeltaDissipation(double newDiss) { tempDeltaDissipation = newDiss; }

    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }


    Interface *giveInterface(InterfaceType) override;

    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;
};
} // end namespace oofem
#endif // matstatus_h
