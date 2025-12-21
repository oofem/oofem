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

#include "mpm.h"
#include "interpolationcatalogue.h"
#include "engngm.h"

namespace oofem {

void
Variable::initializeFrom(InputRecord &ir)
{
    // read interpolation type
    std::string name;
    IR_GIVE_FIELD(ir, name, "interpolation");
    // get corresponding interpolation from catalogue
    this->interpolation = interpolationCatalogue.getInterpolationByName(name);
    
    // read variable type
    ir.giveField(this->type,"type");
    // read quantity
    ir.giveField(this->q,"quantity");
    // read variable size
    IR_GIVE_FIELD(ir, this->size, "size");
    // read dofs 
    IR_GIVE_FIELD(ir, this->dofIDs, "dofs");
}    

void
Term::initializeFrom(InputRecord &ir, EngngModel* problem)
{
    // read field and test field ids (names)
    std::string name;
    IR_GIVE_FIELD(ir, name, "variable");
    this->field = problem->giveVariableByName(name);
    IR_GIVE_FIELD(ir, name, "testvariable");
    this->testField = problem->giveVariableByName(name);
    IR_GIVE_FIELD(ir, this->mode, "mmode");
}



} // end namespace oofem

