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

#ifndef PRIMVARMAPPER_H_
#define PRIMVARMAPPER_H_

#include "valuemodetype.h"

namespace oofem {
class FloatArray;
class Domain;
class TimeStep;

/**
 * Base class for mapping of primary variables between domains.
 * The reason for not using the existing PrimaryUnknownMapper is that I need to handle
 * cases where the number of dofs per node may change (for example XFEM enrichments).
 *
 * Created on: Nov 7, 2013
 * @author Erik Svenning
 */
class PrimaryVariableMapper
{
public:
    PrimaryVariableMapper();
    virtual ~PrimaryVariableMapper();

    virtual void mapPrimaryVariables(FloatArray &oU, Domain &iOldDom, Domain &iNewDom, ValueModeType iMode, TimeStep &iTStep) = 0;
};

/**
 * LSPrimaryVariableMapper: Least-squares primary variable mapper.
 *
 * @author Erik Svenning
 */
class LSPrimaryVariableMapper
{
public:
    LSPrimaryVariableMapper();
    virtual ~LSPrimaryVariableMapper();

    virtual void mapPrimaryVariables(FloatArray &oU, Domain &iOldDom, Domain &iNewDom, ValueModeType iMode, TimeStep &iTStep);
};
} /* namespace oofem */
#endif /* PRIMVARMAPPER_H_ */
