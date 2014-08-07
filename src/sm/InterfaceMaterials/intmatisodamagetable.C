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

#include "intmatisodamagetable.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

#include <algorithm>

namespace oofem {
REGISTER_Material(IntMatIsoDamageTable);

IntMatIsoDamageTable :: IntMatIsoDamageTable(int n, Domain *d) : IntMatIsoDamage(n, d){}



IntMatIsoDamageTable :: ~IntMatIsoDamageTable(){}




IRResultType
IntMatIsoDamageTable :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    std :: ifstream is;
    int nbrOfLinesToRead;

    IR_GIVE_FIELD(ir, kn, _IFT_IntMatIsoDamageTable_kn);
    IR_GIVE_FIELD(ir, ks, _IFT_IntMatIsoDamageTable_ks);
    IR_GIVE_FIELD(ir, ft, _IFT_IntMatIsoDamageTable_ft);

    this->e0 = ft / kn;

    //Set limit on the maximum isotropic damage parameter if needed
    IR_GIVE_OPTIONAL_FIELD(ir, maxOmega, _IFT_IntMatIsoDamageTable_maxOmega);
    maxOmega = min(maxOmega, 0.999999);
    maxOmega = max(maxOmega, 0.0);

    // Parse the table file
    IR_GIVE_FIELD(ir, tablename, _IFT_IntMatIsoDamageTable_tablename);

    is.open(tablename.c_str(), std :: ifstream :: in);
    if ( !is.is_open() ) {
        OOFEM_ERROR(" Can't open table file %s.", tablename.c_str());
    }

    // Read first line
    if ( is >> nbrOfLinesToRead ) {
        printf("NumberofLinestoRead: %d\n", nbrOfLinesToRead);
    } else {
        OOFEM_ERROR("IntMatIsoDamageTable :: initializeFrom: Error reading table file, first line should be "
                    "an integer stating how many strain damage pairs that exist in the file.");
    }

    tableDamages.resize(nbrOfLinesToRead + 1);
    tableJumps.resize(nbrOfLinesToRead + 1);

    // Insert a (0,0) pair
    tableJumps(0) = tableDamages(0) = 0.0;

    for ( int i = 0; i < nbrOfLinesToRead; i++ ) {
        if ( !( is >> tableJumps(i + 1) >> tableDamages(i + 1) ) ) {
            OOFEM_ERROR("Error reading table file at line %d, expected a "
                         "strain damage pair.", i + 2);
        }

        if ( ( tableDamages(i + 1) < tableDamages(i) ) || ( tableJumps(i + 1) < tableJumps(i) ) ) {
            OOFEM_ERROR("Error reading table file at line %d, strain "
                         "and damage must be given in an increasing order and be positive.", i + 2);
        }
    }

    // We add e0 to the tableJumps since these should be given as the increase in
    // strain relative to e0.
    tableJumps.add(e0);

    return StructuralInterfaceMaterial :: initializeFrom(ir);
}


void
IntMatIsoDamageTable :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);

    input.setField(this->kn, _IFT_IntMatIsoDamageTable_kn);
    input.setField(this->ks, _IFT_IntMatIsoDamageTable_ks);
    input.setField(this->ft, _IFT_IntMatIsoDamageTable_ft);
    input.setField(this->tablename, _IFT_IntMatIsoDamageTable_tablename);
    input.setField(this->maxOmega, _IFT_IntMatIsoDamageTable_maxOmega);
}



void
IntMatIsoDamageTable :: computeDamageParam(double &omega, double kappa)
{
    if ( kappa > this->e0 ) {
        // Linear interpolation between able values.

        // If out of bounds damage is set to the last given damage value in the table
        if ( kappa >= tableJumps.at( tableJumps.giveSize() ) ) {
            omega = tableDamages.at( tableDamages.giveSize() );
        } else {
            // std::lower_bound uses binary search to find index with value bounding kappa from above
            int index = (int)(std :: lower_bound(tableJumps.givePointer(), tableJumps.givePointer() + tableJumps.giveSize(), kappa) - tableJumps.givePointer());

#if 0
            printf("e0 %lf\n", e0);
            printf( "sizeof %d\n", sizeof( tableJumps.givePointer() ) );
            printf("pointer: %d\n", index);
            printf("value: %lf\n", * index);
            printf( "Index found: %d\n", index - tableJumps.givePointer() );
            printf( "First index: %d\n", tableJumps.givePointer() );
            printf( "Last index: %d\n", tableJumps.givePointer() + tableJumps.giveSize() );
#endif

            // Pointer arithmetic to find the values used in interpolation
            double x0 = tableJumps(index - 1);
            double x1 = tableJumps(index );
            double y0 = tableDamages(index - 1);
            double y1 = tableDamages(index );

            // Interpolation formula
            omega = y0 + ( y1 - y0 ) * ( kappa - x0 ) / ( x1 - x0 );
        }
    } else {
        omega = 0.0;
    }
}




} // end namespace oofem
