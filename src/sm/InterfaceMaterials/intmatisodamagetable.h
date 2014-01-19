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

#ifndef intmatisodamagetable_h
#define intmatisodamagetable_h


#include "intmatisodamage.h"
#include <fstream>

///@name Input fields for IntMatIsoDamageTable
//@{
#define _IFT_IntMatIsoDamageTable_Name "intmatisodamagetable"
#define _IFT_IntMatIsoDamageTable_kn "kn"
#define _IFT_IntMatIsoDamageTable_ks "ks"
#define _IFT_IntMatIsoDamageTable_ft "ft"
#define _IFT_IntMatIsoDamageTable_tablename "tablename"
#define _IFT_IntMatIsoDamageTable_maxOmega "maxomega"

//@}

namespace oofem {

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
 * @author Jim Brouzoulis
 */

class IntMatIsoDamageTable : public IntMatIsoDamage
{
protected:
    ///Additional parameters
    /// Name of table file
    std :: string tablename;
    /// Damages read from the second column in the table file
    FloatArray tableDamages;
    /// Jumps read from the first column in the table file
    FloatArray tableJumps;

public:
    /// Constructor
    IntMatIsoDamageTable(int n, Domain *d);
    /// Destructor
    virtual ~IntMatIsoDamageTable();

    virtual const char *giveInputRecordName() const { return _IFT_IntMatIsoDamageTable_Name; }

    /**
     * Computes the value of damage parameter omega, based on given value of equivalent strain.
     * It uses a table (jump vs. damage) and interpolates linearly in between.
     * @param[out] omega Contains result.
     * @param kappa Equivalent strain measure.
     */
    virtual void computeDamageParam(double &omega, double kappa);
    
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

};
} // end namespace oofem
#endif // isointerfacedamage01_h
