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

/**
 *	Interface used to allow export modules to export data
 *	from boundary conditions.
 *
 *  @author Erik Svenning
 *
 *  Created on: Feb 12, 2014
 */

#ifndef BCEXPORTINTERFACE_H_
#define BCEXPORTINTERFACE_H_

namespace oofem {

class PrescribedGradient;
class TimeStep;

class BCExportInterface {
public:
	BCExportInterface();
	virtual ~BCExportInterface();

    /**
     * Boundary condition output
     */
    virtual void outputBoundaryCondition(PrescribedGradient &iBC, TimeStep *tStep) = 0;
};

} // end namespace oofem

#endif /* BCEXPORTINTERFACE_H_ */
