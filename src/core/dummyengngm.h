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

#ifndef dummyengngm_h
#define dummyengngm_h

#include "engngm.h"
#include "inputrecord.h"



#define _IFT_DummyEngngModel_Name "Dummy"

namespace oofem {



/**
 * Dummy Engneering model. Does not solve any problem, but invokes the configured export modules.
 * Usefull for exporting model geometry without solving the problem.
 */
class OOFEM_EXPORT DummyEngngModel : public EngngModel
{
public:

    /**
     * Constructor. Creates Engng model with number i.
     */
    DummyEngngModel(int i, EngngModel * _master = NULL);
    virtual ~DummyEngngModel() {}

    void solveYourselfAt(TimeStep *tStep) override;
    void initializeFrom(InputRecord &ir) override;
    TimeStep * giveNextStep() override;
    const char *giveClassName() const override { return "DummyEngngModel"; }
 

};

} // end namespace oofem
#endif // engngm_h
