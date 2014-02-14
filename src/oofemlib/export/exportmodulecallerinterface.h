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
#ifndef EXPORTMODULECALLERINTERFACE_H_
#define EXPORTMODULECALLERINTERFACE_H_

namespace oofem {
/**
 *	Interface used to call an export module from a component
 *	implementing the interface. It allows different output for
 *	different component types, e.g. output from a PrescribedGradient
 *	boundary condition can be treated differently than output from a
 *	"regular" Dirichlet boundary condition.
 *
 *	The interface provides a double dispatch mechanism
 *	to call the correct method.
 *
 *  @author Erik Svenning
 *
 *  Created on: Feb 12, 2014
 */

class TimeStep;

class ExportModule;
class BCExportInterface;

class ExportModuleCallerInterface {
public:
	ExportModuleCallerInterface();
	virtual ~ExportModuleCallerInterface();

    virtual void callExportModule(BCExportInterface &iExpMod, TimeStep *tStep) = 0;

};
} // end namespace oofem

#endif /* EXPORTMODULECALLERINTERFACE_H_ */
