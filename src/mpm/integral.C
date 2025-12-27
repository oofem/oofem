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
#include "integral.h"
#include "mpm.h"
#include "set.h"
#include "intarray.h"
#include "domain.h"
#include "sparsemtrx.h"
#include "engngm.h"

namespace oofem {
void
Integral::initializeFrom (InputRecord &ir, EngngModel *emodel) {
    int di;
    IR_GIVE_FIELD (ir, di, "domain");
    this->domain = emodel->giveDomain(di);
    int si;
    IR_GIVE_FIELD (ir, si, "set");
    this->setIndex  = si;
    this->set = nullptr; // be resolved in initialize
    int ti;
    IR_GIVE_FIELD (ir, ti, "term");
    this->term = emodel->giveTerm(ti);
    
    IR_GIVE_OPTIONAL_FIELD(ir, factor, "factor");
}
        

} // namespace oofem

