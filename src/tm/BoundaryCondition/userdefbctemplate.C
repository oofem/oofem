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

#include "tm/BoundaryCondition/userdefbctemplate.h"
#include "classfactory.h"
#include "node.h"
#include "masterdof.h"
#include "element.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "gausspoint.h"
#include "sparsemtrx.h"
#include "timestep.h"
#include "function.h"
#include "sparselinsystemnm.h"
#include "unknownnumberingscheme.h"
#include "engngm.h"
#include "mathfem.h"
#include "dynamicinputrecord.h"
#include "tm/Elements/transportelement.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace oofem {
REGISTER_BoundaryCondition(UserDefinedBC);

UserDefinedBC :: UserDefinedBC(int n, Domain *d) :
    ActiveBoundaryCondition(n, d)
{
}


void UserDefinedBC :: initializeFrom(InputRecord &ir)
{
    ActiveBoundaryCondition :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, nlocations, _IFT_UserDefinedBC_nlocations);
    locations.resize(nlocations);
    for (int l=1; l<= nlocations; l++) {
        char kwdname[50];
        sprintf(kwdname, "%s%d", _IFT_UserDefinedBC_location, l); // get keyword for lth location
        // initialize l-th location
        IR_GIVE_FIELD(ir, locations[l-1], kwdname);
    }
}


void UserDefinedBC :: giveInputRecord(DynamicInputRecord &input)
{
    //ActiveBoundaryCondition :: giveInputRecord(input);
    input.setField(nlocations, _IFT_UserDefinedBC_nlocations);
    for (int l=1; l<= nlocations; l++) {
        char kwdname[50];
        sprintf(kwdname, "%s%d", _IFT_UserDefinedBC_location, l); // get keyword for lth location
        // initialize l-th location
        input.setField(locations[l-1], kwdname);
    }
}


void UserDefinedBC :: postInitialize()
{
    ActiveBoundaryCondition :: postInitialize();
    
}




void UserDefinedBC :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                                   CharType type, ValueModeType mode,
                                                   const UnknownNumberingScheme &s, 
                                                   FloatArray *eNorm,
                                                   void* lock)
{
    IntArray loc;
    if ( type == ExternalForcesVector ) {
        FloatArray contribution;

        // sample code illustrating accessing solution field(s) at sampling locations
        EngngModel* em = this->domain->giveEngngModel();
        FieldPtr temperatureField = em->giveField(FT_Temperature, tStep);
        for (int l=1; l<= nlocations; l++) {
            FloatArray temperature;
            temperatureField->evaluateAt(temperature, locations[l-1], VM_Total, tStep);
        }

#ifdef _OPENMP
        if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
        answer.assemble(contribution, loc);
#ifdef _OPENMP
        if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif
    } else if ( type == InternalForcesVector ) {
        
    }
}

void UserDefinedBC :: assemble(SparseMtrx &answer, TimeStep *tStep,
                                          CharType type, const UnknownNumberingScheme &r_s, 
                                          const UnknownNumberingScheme &c_s, double scale,
                                          void* lock)
{
    if ( type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == ElasticStiffnessMatrix ) {
        FloatMatrix contribution;
        IntArray loc; 

    #ifdef _OPENMP
        if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
        answer.assemble(loc, loc, contribution); // Contribution assembled
#ifdef _OPENMP
        if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif
        
    } else   {
        printf("Skipping assembly in UserDefinedBC::assemble().\n");
    }
}

void UserDefinedBC :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                                       const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
}

} /* namespace oofem */
