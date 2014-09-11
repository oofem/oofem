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

#ifndef latticematstatus_h
#define latticematstatus_h

#include "../sm/Materials/structuralms.h"

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
class LatticeMaterialStatus : public StructuralMaterialStatus
{
public:
    LatticeMaterialStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor.
    virtual ~LatticeMaterialStatus() { }

    virtual void printOutputAt(FILE *, TimeStep *) { }

    virtual void initTempStatus() { }

    virtual void updateYourself(TimeStep *) { }

    /// Gives the last equilibrated normal stress
    virtual double giveNormalStress() { return 0; }

    /// Gives the last equilibrated normal stress
    virtual double giveOldNormalStress(){return 0;}

    /// Gives the last equilibrated normal stress
    virtual int hasBeenUpdated(){return 0;}


    virtual const char *giveClassName() const { return "LatticeMaterialStatus"; }

    ///Sets the temp_crack_flag
    virtual void setTempCrackFlag(int val) = 0;

    /**
     * Returns the crack flag
     * @return crack flag
     */
    virtual int giveCrackFlag() { return 0; }

    /**
     * @return crack width
     */
    virtual double giveCrackWidth() { return 0; }


    /**
     * @return old crack width
     */
    virtual double giveOldCrackWidth() { return 0; }


    /**
     * Returns the energy dissipation computed at the GaussPoint of the element.
     * This function is used for the lattice specific vtk export.
     * @return dissipation
     */
    virtual double giveDissipation() { return 0; }

    /**
     * Returns the increment of dissipation computed at the GaussPoint of the element.
     * This function is used for the lattice specific vtk export.
     * @return increment of dissipation
     */
    virtual double giveDeltaDissipation() { return 0; }

    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
};
} // end namespace oofem
#endif // matstatus_h
