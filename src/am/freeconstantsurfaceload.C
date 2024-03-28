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

#include "freeconstantsurfaceload.h"
#include "dynamicinputrecord.h"
#include "function.h"
#include "floatarray.h"
#include "intarray.h"
#include "timestep.h"
#include "classfactory.h"
#include "feinterpol3d.h"

#include "connectivitytable.h"
#include <iostream>
#include <map>

#include <algorithm>

namespace oofem {
REGISTER_BoundaryCondition(FreeConstantSurfaceLoad);

FreeConstantSurfaceLoad :: FreeConstantSurfaceLoad(int i, Domain *d) : SurfaceLoad(i, d)
{
    this->loadOffset = 0.0;
}

void
FreeConstantSurfaceLoad :: initializeFrom(InputRecord &ir)
{
    Load :: initializeFrom(ir);

    int dummy;
    IR_GIVE_OPTIONAL_FIELD(ir, dummy, "ndofs");

    int value = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_BoundaryLoad_loadtype);
    lType = ( bcType ) value;

    value = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_BoundaryLoad_cstype);
    coordSystemType = ( CoordSystType ) value;

    IR_GIVE_OPTIONAL_FIELD(ir, propertyDictionary, _IFT_BoundaryLoad_properties);
    IR_GIVE_OPTIONAL_FIELD(ir, propertyTimeFunctDictionary, _IFT_BoundaryLoad_propertyTimeFunctions);

    IR_GIVE_OPTIONAL_FIELD(ir, propertyMultExpr, _IFT_BoundaryLoad_propertyMultExpr);

    //temperOffset = 0.0;
    temperOffset = 273.15;
    IR_GIVE_OPTIONAL_FIELD(ir, temperOffset, _IFT_BoundaryLoad_temperOffset);

    IR_GIVE_OPTIONAL_FIELD(ir, this->loadOffset, _IFT_FreeConstantSurfaceLoad_LoadOffset);
}

void
FreeConstantSurfaceLoad :: giveInputRecord(DynamicInputRecord &input)
{
    SurfaceLoad :: giveInputRecord(input);
    input.setField(this->loadOffset, _IFT_FreeConstantSurfaceLoad_LoadOffset);
}

void
FreeConstantSurfaceLoad :: computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode)
{
    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        OOFEM_ERROR("mode not supported");
    }

    double factor = this->giveTimeFunction()->evaluate(tStep, mode);
    answer.beScaled(factor, componentArray);
}

bool FreeConstantSurfaceLoad :: isImposed(TimeStep *tStep)
{
    if(this->alreadyDisabled) return false;

    //std::cout << "start\n";
    IntArray boundaries = domain->giveSet(set)->giveBoundaryList();

    // !!! only elementboundaries with single element and single surface accepted at the moment !!! 
    Element *el = domain->giveElement(boundaries.at(1));
    IntArray surfNodes = el->giveBoundarySurfaceNodes(boundaries.at(2));
    int nSurfNodes = surfNodes.giveSize();

    int nElements = domain->giveNumberOfElements();
    
    // map of node occurance
    std::map<int, int> connmap;
    // iterate over BC's surface nodes
    for(int k = 1; k <= nSurfNodes; k++) {
        auto nodeId = el->giveDofManagerNumber(surfNodes.at(k));

        // connected element ids
        auto conns = domain->giveConnectivityTable()->giveDofManConnectivityArray(nodeId);
        int connSize = conns->giveSize();
        //std::cout << "\nel=" << el->giveNumber() << ", ";
        //std::cout << connSize << ", ";
        //std::cout << conns->at(1) << "\n";

        // TODO: ITS WRONG WITH INDEPENDENT ELEMENTS?
        if(connSize > 8 || connSize == 0) return true;

        for(int c = 1; c <= connSize; c++) {
            if(el->giveNumber() != conns->at(c)) {
                if (connmap.find(conns->at(c)) == connmap.end()) {
                    connmap[conns->at(c)] = 1;
                } else {
                    connmap[conns->at(c)] = connmap[conns->at(c)]+1;
                }
            }
            //std::cout << conns->at(c) << ", ";
        }
        //giveDofManager(nodeId)->;
    }

    //std::cout << "end\n";

    // if some element is connected to 4 nodes, its connected to the whole face!
    for(auto &el : connmap) {
        if(el.second == 4) {
            this->alreadyDisabled = true;
            return false; 
        }
    }

    return true;
}
} // end namespace oofem
