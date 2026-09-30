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

#ifndef mpmassemblers_h
#define mpmassemblers_h

#include "primaryfield.h"
#include "function.h"
#include "dofdistributedprimaryfield.h"
#include "assemblercallback.h"

namespace oofem {

/**
 * Callback class for assembling mid point effective tangents. 
 * @todo Need to parametrize individual contributing terms, ther locations and multilication factors.
 */
class UPLhsAssembler : public MatrixAssembler
{
protected:
    double alpha;
    double deltaT;

public:
    UPLhsAssembler(double alpha, double deltaT);
    void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const override;
};

/**
 * Callback class for assembling residuals
 */
class UPResidualAssembler : public VectorAssembler
{
    protected:
    double alpha;
    double deltaT;
public:
    UPResidualAssembler(double alpha, double deltaT) : VectorAssembler(), alpha(alpha), deltaT(deltaT) {}
    void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const override;
};

/**
 * Callback class for assembling mid point effective tangents. 
 * @todo Need to parametrize individual contributing terms, ther locations and multilication factors.
 */
class TMLhsAssembler : public MatrixAssembler
{
protected:
    double alpha;
    double deltaT;

public:
    TMLhsAssembler(double alpha, double deltaT);
    void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const override;
};

/**
 * Callback class for assembling residuals
 */
class TMResidualAssembler : public VectorAssembler
{
    protected:
    double alpha;
    double deltaT;
public:
    TMResidualAssembler(double alpha, double deltaT) : VectorAssembler(), alpha(alpha), deltaT(deltaT) {}
    void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const override;
};


} // end namespace oofem
#endif // mpmassemblers_h
